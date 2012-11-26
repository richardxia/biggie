package snap.apps

import snap._
import it.unimi.dsi.fastutil.longs.LongList
import it.unimi.dsi.fastutil.longs.LongArrayList
import spark.SparkContext
import SparkContext._

object SimFinder {
  def main(args : Array[String]): Unit = {
    // load arguments
    args match {
      case Array(fastaFile, seedLenStr, readLenStr, numSeedsStr, unionDistStr, outFile, nClusterStr) => {
        val seedLen = seedLenStr.toInt
        val readLen = readLenStr.toInt
        val numSeeds = numSeedsStr.toInt
        val unionDist = unionDistStr.toInt
        val nClusters = nClusterStr.toInt
        
        assert(numSeeds > 1)
        
        // open output file
        val out = new java.io.File(outFile)
        if (!out.exists)
          out.createNewFile
          
        val fw = new java.io.FileWriter(out.getName)
        val bw = new java.io.BufferedWriter(fw)

        // load genome
        val genome = FASTA.read(fastaFile)
        val numBases = genome.totalSize

        // build index (at some point, we want to load it)
        val builder = new IntHashIndexBuilder(seedLen, numBases)
        genome.addToIndex(builder)
        val index = builder.build()

        // initialize variables
        var validPositions = 0
        val valid = Array.fill(numBases.toInt)(false)
        var hits = new LongArrayList
        val unionFind = new UnionFind(numBases)

        // for each substring, find its hits, and join it via union-find to any similar strings
        var pos = 0L
        var i = 0
        var j = 0
        val posStr = new Array[Byte](readLen)
        val hitPosStr = new Array[Byte](readLen)
        
        while (pos < numBases) {
          if (pos % 100000 == 0)
            println("Position: " + pos + "(" + math.round(pos * 100.0) / numBases + "%)")
            
          genome.getSubstring(pos, pos + readLen, posStr)
          
          if (!containsN(posStr)) {
            validPositions += 1
            valid(pos.toInt) = true
            
            // try each seed from the read
            i = 0
            while (i < numSeeds) {
              val seedPos = (i * (readLen - seedLen)) / (numSeeds - 1)
              val seed = DNA.substringToLong(posStr, seedPos, seedPos + seedLen)
              hits.clear()
              index.get(seed, hits, Int.MaxValue)
              val nHits = hits.size()

              j = nHits - 1
              while (j >= 0) {
                val hitPos = hits.getLong(j) - seedPos
                genome.getSubstring(hitPos, hitPos + readLen, hitPosStr)
                
                if (hitPos > pos && unionFind.find(hitPos.toInt) != unionFind.find(pos.toInt)) {
                  if (!containsN(hitPosStr)) {
                    if (HammingDistance.distance(posStr, hitPosStr, unionDist) != -1) {
                      unionFind.union(pos.toInt, hitPos.toInt)
                    }
                  }
                } else if (hitPos <= pos) {
                  j = -1
                }
                      
                j -= 1
              }
              
              i += 1
            }
          }
          
          pos += 1
        }

        // finalize clusters
        unionFind.findClusters(readLen, genome.validSubstring)

        // print cluster stats
        unionFind.printClusterStats(validPositions)

        // print clusters to file
        unionFind.writeClustersToFile(bw, genome, readLen)
        
        // print kernel matrix & cluster members
        unionFind.writeTopNKernelMatricesToFile(outFile, genome, readLen, nClusters)
        unionFind.writeTopNClusterMembersToFile(outFile, genome, readLen, nClusters)
                
        bw.close
      }
      case Array("parallel", master, seedLenStr, readLenStr, numSeedsStr, unionDistStr, outFile, numTasksStr, minClusterSizeStr) => {
        val startRun = System.currentTimeMillis
        
        val seedLen = seedLenStr.toInt
        val readLen = readLenStr.toInt
        val numSeeds = numSeedsStr.toInt
        val unionDist = unionDistStr.toInt
        val numTasks = numTasksStr.toInt
        val minClusterSize = minClusterSizeStr.toInt

        assert(numSeeds > 1)
        //assert(numTasks > 1)
        
        // open output file
        val out = new java.io.File(outFile)
        if (!out.exists)
          out.createNewFile
          
        val fw = new java.io.FileWriter(out.getName)
        val bw = new java.io.BufferedWriter(fw)
        
        // load genome
        val startLoad = System.currentTimeMillis
        val numBases = GenomeLoader.genome.totalSize
        val rangeLength = (numBases / (numTasks * 2)).toInt  // figure out if ceiling is necessary
        println("Genome load took " + math.round(System.currentTimeMillis - startLoad) / 1000.0 + "s")

        val sc = new SparkContext(getSparkDest(master), "SimFinder")

        val startMap = System.currentTimeMillis

        println("Ranges:")
        (0 until (numBases / 2).toInt by rangeLength).foreach(rangeStart => {
          val startPos = 2 * rangeStart
          val endPos = math.min(2 * rangeStart + 2 * rangeLength, numBases) - 1
          println(startPos + ", " + endPos)
        })

        val clusterGroups = sc.parallelize(0 until (numBases / 2).toInt by rangeLength).map(rangeStart => {
          
          val startPos = 2 * rangeStart
          val endPos = math.min(2 * rangeStart + 2 * rangeLength, numBases) - 1

          // build index for given range
          val builder = new IntHashIndexBuilder(seedLen, rangeLength)
          GenomeLoader.genome.addToIndex(builder, startPos, endPos)
          val index = builder.build()

          // initialize variables
          var validPositions = 0
          val midPoint = (numBases / 2).toInt

          //val valid = Array.fill(numBases.toInt)(false)
          val validL = Array.fill(midPoint)(false)
          val validH = Array.fill(midPoint)(false)
          
          var hits = new LongArrayList
          val unionFindL = new UnionFindL(numBases)
          
          // for each substring, find its hits, and join it via union-find to any similar strings
          var pos = 0L
          var i = 0
          var j = 0
          val posStr = new Array[Byte](readLen)
          val hitPosStr = new Array[Byte](readLen)
          
          while (pos < numBases) {
            if (pos % 100000 == 0)
              println("Position: " + pos + "(" + math.round(pos * 100.0) / numBases + "%)")
            
            GenomeLoader.genome.getSubstring(pos, pos + readLen, posStr)
            
            if (!containsN(posStr)) {
              validPositions += 1
              if (pos < midPoint)
                validL(pos.toInt) = true
              else
                validH((pos - midPoint).toInt) = true
            
              // try each seed from the read
              i = 0
              while (i < numSeeds) {
                val seedPos = (i * (readLen - seedLen)) / (numSeeds - 1)
                val seed = DNA.substringToLong(posStr, seedPos, seedPos + seedLen)
                hits.clear()
                index.get(seed, hits, Int.MaxValue)
                val nHits = hits.size()
              
                j = nHits - 1
                while (j >= 0) {
                  val hitPos = hits.getLong(j) - seedPos
                  GenomeLoader.genome.getSubstring(hitPos, hitPos + readLen, hitPosStr)
                
                  if (hitPos > pos && unionFindL.find(hitPos) != unionFindL.find(pos)) {
                    if (!containsN(hitPosStr)) {
                      if (HammingDistance.distance(posStr, hitPosStr, unionDist) != -1) {
                        unionFindL.union(pos, hitPos)
                      }
                    }
                  } else if (hitPos <= pos) {
                    j = -1
                  }
                
                  j -= 1
                }
              
                i += 1
              }
            }
                        
            pos += 1
          }
          
          // return clusters
          //unionFindL.findClusters(readLen, GenomeLoader.genome.validSubstring)  // TODO:  shouldn't use validSubstring (super slow) -- should use !containsN
          unionFindL.findClusters(readLen, containsN)
          unionFindL
        }).collect()

        println("Map took " + math.round(System.currentTimeMillis - startMap) / 1000.0 + "s")

        // debug
        // TODO:  remove this
        println
        println("Intermediate results:")
        val ranges = (0 until (numBases / 2).toInt by rangeLength)
        (0 until ranges.length).foreach(i => {
          val rangeStart = ranges(i)
          val startPos = 2 * rangeStart
          val endPos = math.min(2 * rangeStart + 2 * rangeLength, numBases) - 1
          println(startPos + ", " + endPos)
          clusterGroups(i).printClusterStats(validPositions(readLen))
          println
        })

        // merge clusters
        val startMerge = System.currentTimeMillis
        
        val unionFindLClusters = new UnionFindL(numBases)
        clusterGroups.foreach(g => {
          //val (startPos, clusters) = g
          val clusters = g.nonTrivialMembers

          clusters.keySet.foreach(k => {
            val members = clusters.get(k).get
            val firstMember = members.getLong(0)
            var m = 1
            while (m < members.size) {
              val currentMember = members.getLong(m)
              if (unionFindLClusters.find(firstMember) != unionFindLClusters.find(currentMember))
                unionFindLClusters.union(firstMember, currentMember)
                
              m += 1
            }
          })
        })
        
        println("Merge took " + math.round(System.currentTimeMillis - startMerge) / 1000.0 + "s")
        
        // finalize clusters
        val startFind = System.currentTimeMillis
        //unionFindLClusters.findClusters(readLen, GenomeLoader.genome.validSubstring)
        unionFindLClusters.findClusters(readLen, containsN)
        println("Find took " + math.round(System.currentTimeMillis - startFind) / 1000.0 + "s")

        val startAdmin = System.currentTimeMillis
        // print cluster stats
        unionFindLClusters.printClusterStats(validPositions(readLen))

        // print clusters to file
        unionFindLClusters.writeClustersToFileSnapCompliant(bw, GenomeLoader.genome, readLen, unionDist, minClusterSize)
        println("Admin took " + (System.currentTimeMillis - startAdmin) / 1000.0 + "s")
        
        bw.close
        
        println("Run took " + (System.currentTimeMillis - startRun) / 1000.0 + "s")
      }
      case Array("gridParallel", master, seedLenStr, readLenStr, numSeedsStr, unionDistStr, outFile, gridDimStr, minClusterSizeStr) => {
        val startRun = System.currentTimeMillis
        
        val params = new GridParallelSimFinderParams(master, seedLenStr, readLenStr, numSeedsStr, unionDistStr, outFile, gridDimStr, minClusterSizeStr)
        
        // open output file
        val out = new java.io.File(outFile)
        if (!out.exists)
          out.createNewFile
          
        val fw = new java.io.FileWriter(out.getName)
        val bw = new java.io.BufferedWriter(fw)
                
        val numBases = GenomeLoader.genome.totalSize
        val rangeLength = (numBases / (params.gridDim * 2)).toInt  // division by 2 here makes things confusing :(
        val rangeStarts = (0 until (numBases / 2).toInt by rangeLength)

        println("Ranges:")
        rangeStarts.foreach(rangeStart => {
          val startPos = 2 * rangeStart
          val endPos = math.min(2 * rangeStart + 2 * rangeLength, numBases) - 1
          println(startPos + ", " + endPos)
        })

        val partitionStarts = rangeStarts.map(i => List(i).padTo(rangeStarts.length, i).zip(rangeStarts)).flatten
        //val partitionStarts = getUniquePartitions(numBases, gridDim)  // actually, you want all the partitions

        val sc = new SparkContext(getSparkDest(master), "SimFinder", "/root/spark", Seq("target/scala-2.9.2/snap_2.9.2-0.0.jar"))
        val clusterGroups = sc.parallelize(partitionStarts).map(p => {
          getPartitionClusters(params, p, rangeLength)
        }).collect()
        
        // Merge clusters
        val ufClusters = new UnionFindL(numBases) // must use "L" here b/c we're dealing with whole genome

        clusterGroups.foreach(g => {
          val clusters = g.getNonTrivialMembers

          clusters.keySet.foreach(k => {
            val members = clusters.get(k).get
            val firstMember = members.getLong(0)  // convert to pos; requires knowing whether this was grid or grid diagonal
            var m = 1
            while (m < members.size) {
              val currentMember = members.getLong(m)  // convert to pos; requires knowing whether this was grid or grid diagonal
              if (ufClusters.find(firstMember) != ufClusters.find(currentMember))
                ufClusters.union(firstMember, currentMember)
                
              m += 1
            }
          })
        })

        //ufClusters.findClusters(params.readLen, GenomeLoader.genome.validSubstring)
        ufClusters.findClusters(params.readLen, containsN)
        
        // Print cluster stats
        ufClusters.printClusterStats(validPositions(params.readLen))
        
        // Print clusters to file
        ufClusters.writeClustersToFileSnapCompliant(bw, GenomeLoader.genome, params.readLen, params.unionDist, params.minClusterSize)
        
        bw.close
        
        println("Run took " + (System.currentTimeMillis - startRun) / 1000.0 + "s")
      }
      case Array("gridParallelAccumulator", master, seedLenStr, readLenStr, numSeedsStr, unionDistStr, outFile, gridDimStr, minClusterSizeStr) => {
        val startRun = System.currentTimeMillis
        
        val params = new GridParallelSimFinderParams(master, seedLenStr, readLenStr, numSeedsStr, unionDistStr, outFile, gridDimStr, minClusterSizeStr)
        
        // open output file
        val out = new java.io.File(outFile)
        if (!out.exists)
          out.createNewFile
          
        val fw = new java.io.FileWriter(out.getName)
        val bw = new java.io.BufferedWriter(fw)
                
        val numBases = GenomeLoader.genome.totalSize
        val rangeLength = (numBases / (params.gridDim * 2)).toInt  // division by 2 here makes things confusing :(
        val rangeStarts = (0 until (numBases / 2).toInt by rangeLength)

        println("Ranges:")
        rangeStarts.foreach(rangeStart => {
          val startPos = 2 * rangeStart.toLong
          val endPos = math.min(2 * rangeStart.toLong + 2 * rangeLength.toLong, numBases) - 1
          println(startPos + ", " + endPos)
        })

        val partitionStarts = rangeStarts.map(i => List(i).padTo(rangeStarts.length, i).zip(rangeStarts)).flatten
        //val partitionStarts = getUniquePartitions(numBases, gridDim)  // actually, you want all the partitions
        var ufClusters = new UnionFindL(numBases) // must use "L" here b/c we're dealing with whole genome

	//System.setProperty("spark.master.host", "169.229.49.195") // for gene2
	//System.setProperty("spark.hostname", "169.229.49.195") // for gene2
        val sc = new SparkContext(getSparkDest(master), "SimFinder", "", Seq("target/scala-2.9.2/snap_2.9.2-0.0.jar"))
        val ufAccumulator = sc.accumulator(ufClusters.asInstanceOf[UnionFindAbstract])(UnionFindAP)

        /*val clusterGroups =*/ sc.parallelize(partitionStarts, partitionStarts.size).foreach(p => {
          val uf = getPartitionClusters(params, p, rangeLength)
          ufAccumulator += uf
        })//.collect()

	println("Finished spark part")

        /*
        // Merge clusters
        clusterGroups.foreach(g => {
          val clusters = g.getNonTrivialMembers

          clusters.keySet.foreach(k => {
            val members = clusters.get(k).get
            val firstMember = members.getLong(0)  // convert to pos; requires knowing whether this was grid or grid diagonal
            var m = 1
            while (m < members.size) {
              val currentMember = members.getLong(m)  // convert to pos; requires knowing whether this was grid or grid diagonal
              if (ufClusters.find(firstMember) != ufClusters.find(currentMember))
                ufClusters.union(firstMember, currentMember)
                
              m += 1
            }
          })
        })
        */

	println("About to find clusters on ufClusters...")

        //ufClusters.findClusters(params.readLen, GenomeLoader.genome.validSubstring)
	ufClusters = ufAccumulator.value.asInstanceOf[UnionFindL]
        ufClusters.findClusters(params.readLen, containsN)

	println("Finished finding clusters; about to print stats & to file...")
        
        // Print cluster stats
        ufClusters.printClusterStats(validPositions(params.readLen))
        
        // Print clusters to file
        ufClusters.writeClustersToFileSnapCompliant(bw, GenomeLoader.genome, params.readLen, params.unionDist, params.minClusterSize)
        
        bw.close
        
        println("Run took " + (System.currentTimeMillis - startRun) / 1000.0 + "s")
      }
      case Array("compareAlignmentResults", samFile1, samFile2, outputPrefix) => {
        // input is two sam files
        // assumes reads will be in same order in both files
        // load corresponding reads
        // check:  are they different?  are they correct?
        val samIt1 = SAM.read(samFile1)
        val samIt2 = SAM.read(samFile2)
        
        var onlyFirstCorrect: List[SAMEntry] = Nil
        var onlySecondCorrect: List[SAMEntry] = Nil
        var neitherCorrect: List[SAMEntry] = Nil
        
        while (samIt1.hasNext && samIt2.hasNext) {
          val sam1 = samIt1.next
          val sam2 = samIt2.next
          
          val firstCorrect = Wgsim.isCorrect(sam1)
          val secondCorrect = Wgsim.isCorrect(sam2)
          
          if (firstCorrect && !secondCorrect) onlyFirstCorrect = sam1 :: onlyFirstCorrect
          else if (!firstCorrect && secondCorrect) onlySecondCorrect = sam1 :: onlySecondCorrect
          else if (!firstCorrect && !secondCorrect) neitherCorrect = sam1 :: neitherCorrect
        }
        
        onlyFirstCorrect = onlyFirstCorrect.reverse
        onlySecondCorrect = onlySecondCorrect.reverse
        neitherCorrect = neitherCorrect.reverse
        
        // create fastq writers
        val onlyFirstCorrectWriter = FASTQ.writer(outputPrefix + "_onlyFirstCorrect.fastq")
        val onlySecondCorrectWriter = FASTQ.writer(outputPrefix + "_onlySecondCorrect.fastq")
        val neitherCorrectWriter = FASTQ.writer(outputPrefix + "_neitherCorrect.fastq")

        // write out to files
        onlyFirstCorrect.map(e => new Read(e.readId.getBytes, e.sequence.getBytes, e.quality.getBytes)).foreach(r => onlyFirstCorrectWriter.write(r))
        onlySecondCorrect.map(e => new Read(e.readId.getBytes, e.sequence.getBytes, e.quality.getBytes)).foreach(r => onlySecondCorrectWriter.write(r))
        neitherCorrect.map(e => new Read(e.readId.getBytes, e.sequence.getBytes, e.quality.getBytes)).foreach(r => neitherCorrectWriter.write(r))
        
        onlyFirstCorrectWriter.close
        onlySecondCorrectWriter.close
        neitherCorrectWriter.close
      }
      case _ => println("Incorrect parameters")
    }
  }
  
  def containsN(str: Array[Byte]): Boolean = {
    var i = 0
    while (i < str.length) {
      if (str(i) == 'N')
        return true
      i += 1
    }
    return false
  }
  
  def hammingDistance(s1: String, s2: String, maxDistance: Int): Int = {
    assert(s1.length == s2.length)

    var distance = 0
    
    (0 until s1.length).foreach(i => {
      if (s1(i) != s2(i)) {
        distance += 1
        if (distance > maxDistance)
          return -1
      }
    })
    
    distance
  }

  def validPositions(readLen: Int): Long = {
    var pos = 0L
    var validPositions = 0L
     
    val numBases = GenomeLoader.genome.totalSize
    val posStr = new Array[Byte](readLen)

    while (pos < numBases) {
      GenomeLoader.genome.getSubstring(pos, pos + readLen, posStr)
      
      if (!containsN(posStr))
        validPositions += 1

      pos += 1
    }
    
    validPositions
  }

  def getSparkDest(master: String): String = {
    if (master.startsWith("local"))
      master
    else
      "1@" + master + ":5050"
  }

  /*
  def getUniquePartitions(numBases: Long, gridDim: Int): List[(Long, Long)] = {
    val rangeLength = (numBases / (gridDim * 2)).toInt
    val rangeStarts = (0 until (numBases / 2).toInt by rangeLength).map(2L * _ /* convert to startPos */)

    /*
    println("Ranges:")
    rangeStarts.foreach(rangeStart => {
      val startPos = 2 * rangeStart
      val endPos = min(2 * rangeStart + 2 * rangeLength, numBases) - 1
      println(startPos + ", " + endPos)
    })
    */

    //val partitionStarts = rangeStarts.map(i => List(i).padTo(rangeStarts.length, i).zip(rangeStarts)).flatten
    var partitionStarts: List[(Long, Long)] = Nil
    var i = 0
    var j = 0
    while (i < rangeStarts.length) {
      j = i
      while (j < rangeStarts.length) {
        partitionStarts = ((rangeStarts(i), rangeStarts(j))) :: partitionStarts
        j += 1
      }
      i += 1
    }
    
    partitionStarts.reverse
  }
  */
  
  def getPartitionClusters(params: ParallelSimFinderParams, p: (Int, Int), rangeLength: Int) = {
    val numBases = GenomeLoader.genome.totalSize
    
    val indexRangeStart = p._1  // range start for row, ie, which positions to index
    val scanRangeStart = p._2   // range start for column, ie, which positions to scan over for clustering
    
    val indexStartPos = 2 * indexRangeStart.toLong
    val indexEndPos = math.min(indexStartPos + 2 * rangeLength.toLong, numBases) - 1

    val scanStartPos = 2 * scanRangeStart.toLong
    val scanEndPos = math.min(scanStartPos + 2 * rangeLength.toLong, numBases) - 1
    
    println("indexing " + indexStartPos + " to " + indexEndPos + "; scanning " + scanStartPos + " to " + scanEndPos)

    // build index for given range
    val builder = new IntHashIndexBuilder(params.seedLen, 2 * rangeLength)
    GenomeLoader.genome.addToIndex(builder, indexStartPos, indexEndPos)
    val index = builder.build()
    
    // initialize variables
    var validPositions = 0
    
    // TODO:  keep track of valid positions -- a bit difficult
    // need to do similar thing to what i did with grid, grid diagonal versions of union find (different behavior if range is same or not)
    // TODO:  I think I can skip it, if I just use containsN
    
    var hits = new LongArrayList
    
    val uf = 
    if (indexStartPos == scanStartPos) new UnionFindGridDiagonal((indexStartPos, indexEndPos))
    else new UnionFindGrid((indexStartPos, indexEndPos), (scanStartPos, scanEndPos))
    
    // for each substring, find its hits, and join it via union-find to any similar strings
    var pos = scanStartPos
    var i = 0
    var j = 0
    val posStr = new Array[Byte](params.readLen)
    val hitPosStr = new Array[Byte](params.readLen)
    
    while (pos < scanEndPos) {
      if (pos % 100000 == 0)
        println("Position: " + pos + "(" + math.round((pos - scanStartPos) * 100.0) / (2 * rangeLength) + "%)")
      
      GenomeLoader.genome.getSubstring(pos, pos + params.readLen, posStr)
      
      if (!containsN(posStr)) {
        //validPositions += 1
        // TODO:  update valid positions array
        // TODO:  figure out if I can skip this entirely
        
        // try each seed from the read
        i = 0
        while (i < params.numSeeds) {
          val seedPos = (i * (params.readLen - params.seedLen)) / (params.numSeeds - 1)
          val seed = DNA.substringToLong(posStr, seedPos, seedPos + params.seedLen)
          hits.clear()
          index.get(seed, hits, Int.MaxValue)
          val nHits = hits.size()
          
          j = nHits - 1
          while (j >= 0) {
            val hitPos = hits.getLong(j) - seedPos
            GenomeLoader.genome.getSubstring(hitPos, hitPos + params.readLen, hitPosStr)

      	    // hitPos could be out of range if a seed that occurs early in the range occurs late in the read (makes read extend into prior range)
      	    // so, make sure hitPos is within range
            if (hitPos >= indexStartPos && hitPos > pos && uf.find(hitPos) != uf.find(pos)) {
              if (!containsN(hitPosStr)) {
                if (HammingDistance.distance(posStr, hitPosStr, params.unionDist) != -1) {
                  uf.union(pos, hitPos)
                }
              }
            } else if (hitPos <= pos) {
              j = -1
            }
            
            j -= 1
          }
          
          i += 1
        }
      }
      
      pos += 1
    }
    
    // return clusters
    uf.findClusters(params.readLen, containsN)
    uf
  }
}
