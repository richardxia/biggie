package snap

import scala.io.Source
import it.unimi.dsi.fastutil.longs.LongList
import it.unimi.dsi.fastutil.longs.LongArrayList

class ClusterAnalyzer(clusterFile: String) {
  // Facilitate analysis of clusters produced by SimFinder
  var nonTrivialMembers = scala.collection.mutable.Map[Long, LongArrayList]()
  var clusterIdsBySize = scala.collection.mutable.Map[Int, LongArrayList]()
  var genomeSize = 0L
  var nClusters = 0
  var readLen = 0
  var mergeDist = 0
  var maxClusterSize = 0
  var genome: Genome = null
  
  read(clusterFile)
  
  // read in cluster file
  def read(file: String) = {
    val lines = Source.fromFile(file).getLines
    
    // read in header
    val header = 
    if (lines.hasNext) {
      val l = lines.next
      val entries = l.split(" ")
      assert(entries.length == 4)
      genomeSize = entries(0).toLong
      nClusters = entries(1).toInt
      readLen = entries(2).toInt
      mergeDist = entries(3).toInt
      
      l
    } else ""

    // read in clusters
    for (cluster <- lines) {
      // parse cluster line
      // first # is cluster size
      // subsequent #s are positions belonging to the cluster
      val entries = cluster.split(" ")
      assert(entries.length >= 3)
      val clusterSize = entries(0).toInt
      val clusterId = entries(1).toLong
      val clusterMembers = new LongArrayList

      if (clusterSize > maxClusterSize)
        maxClusterSize = clusterSize

      var i = 1
      while (i < entries.size) {
        clusterMembers.add(entries(i).toLong)
        i += 1
      }
      
      // store cluster
      nonTrivialMembers += ((clusterId, clusterMembers))
      
      // store cluster id by size
      clusterIdsBySize.get(clusterSize) match {
        case None => {
          val l = new LongArrayList
          l.add(clusterId)
          clusterIdsBySize += ((clusterSize, l))
        }
        case Some(l) => l.add(clusterId)
      }
    }
  }
  
  def loadGenome(fastaFile: String) = {
    genome = FASTA.read(fastaFile)
  }
  
  // analyze clusters
  def nClustersWithSize(size: Int): Int = {
    clusterIdsBySize.get(size) match {
      case None => 0
      case Some(l) => l.size
    }
  }

  def nClustersWithSize(min: Int, max: Int): Int = {
    (min to max).map(nClustersWithSize(_)).reduce(_ + _)
  }

  def nClustersWithSizeAtLeast(size: Int): Int = {
    if (size <= maxClusterSize)
      nClustersWithSize(size, maxClusterSize)
    else
      0
  }
  
  def nPosInClustersWithSize(size: Int): Int = {
    clusterIdsBySize.get(size) match {
      case None => 0
      case Some(l) => l.size /* how many clusters with given size */ * size
    }
  }
  
  def nPosInClustersWithSize(min: Int, max: Int): Int = {
    (min to max).map(nPosInClustersWithSize(_)).reduce(_ + _)
  }
  
  def nPosInClustersWithSizeAtLeast(size: Int): Int = {
    if (size <= maxClusterSize)
      nPosInClustersWithSize(size, maxClusterSize)
    else
      0
  }

  def topNClusterIds(n: Int): List[Long] = {
    val effectiveN = math.min(n, nClusters)

    // returns ids of N largest clusters
    // get sorted list of cluster sizes
    val sizes = clusterIdsBySize.keySet.toList.sortWith(_ > _)
    
    var topIds: List[Long] = Nil
    
    var i = 0
    while (topIds.size < effectiveN) {
      val ids = clusterIdsBySize.get(sizes(i)).get
      val nToAdd = 
      if (ids.size < effectiveN - topIds.size)
        ids.size
      else
        effectiveN - topIds.size
      (0 until nToAdd).foreach(j => {
        topIds = ids.getLong(j) :: topIds
      })
      
      i += 1
    }

    topIds.reverse
  }
  
  def getConsensusStr(clusterMembers: Array[Long])/*: String*/ = {
    val counts = Array.fill(readLen)(Array.fill(4)(0))
    val s = new Array[Byte](readLen)
    
    // for each member in cluster, get genome substring
    (0 until clusterMembers.size).foreach(i => {
      genome.getSubstring(clusterMembers(i), clusterMembers(i) + readLen, s)
      (0 until s.length).foreach(j => {
        counts(j)(DNA.BASE_TO_CODE(s(j))) += 1
      })
    })
    
    // get consensus string
    val consensus = Array.fill(readLen)(' '.toByte)
    (0 until readLen).foreach(i => {
      consensus(i) = DNA.CODE_TO_BASE(counts(i).indexOf(counts(i).max)).toByte
    })
    //new String(consensus)
    consensus
  }
  
  def getHistogramOfDiffsToConsensus(clusterMembers: Array[Long]): Array[Int] /* may want to make a histogram class at some point */ = {
    assert(genome != null)

    val consensus = getConsensusStr(clusterMembers)
    val histogram = Array.fill(readLen)(0)
    val member = new Array[Byte](readLen)
    
    (0 until clusterMembers.size).foreach(i => {
      // get ith member
      genome.getSubstring(clusterMembers(i), clusterMembers(i) + readLen, member)
      
      // get edit distance between ith member & consensus
      val editDist = HammingDistance.distance(consensus, member, readLen)
      
      // update histogram
      histogram(editDist) += 1
    })
    
    histogram
  }
  
  def writeKernelMatrixToFile(clusterId: Long, filenamePrefix: String) = {
    val members = nonTrivialMembers.get(clusterId)
    members match {
      case None => println("Invalid cluster ID")
      case Some(l) => {
        val fw = new java.io.FileWriter(filenamePrefix + clusterId + ".matrix")
        val bw = new java.io.BufferedWriter(fw)

        val buf1 = new Array[Byte](readLen)
        val buf2 = new Array[Byte](readLen)
        var j = 0

        // Compute edit distance for each pair in cluster
        while (j < l.size) {
          genome.getSubstring(l.getLong(j), l.getLong(j) + readLen, buf1)
          var k = 0
          while (k < l.size) {
            genome.getSubstring(l.getLong(k), l.getLong(k) + readLen, buf2)
            val dist = HammingDistance.distance(buf1, buf2, readLen)
            bw.write(" " + dist)
        
            k += 1
          }
          bw.newLine
      
          j += 1
        }
        
        bw.close
      }
    }
  }
  
  def getClusterDiameter(clusterMembers: Array[Long]): Int = {
    // cluster diameter is defined as max pairwise distance
    var diameter = 0
    val buf1 = new Array[Byte](readLen)
    val buf2 = new Array[Byte](readLen)

    (0 until clusterMembers.size).foreach(i => {
      (i until clusterMembers.size).foreach(j => {
        // get distance between i & j
        val pos1 = clusterMembers(i)
        val pos2 = clusterMembers(j)

        genome.getSubstring(pos1, pos1 + readLen, buf1)
        genome.getSubstring(pos2, pos2 + readLen, buf2)
        
        val dist = HammingDistance.distance(buf1, buf2, readLen)
        
        if (dist > diameter)
          diameter = dist
      })
    })
    
    diameter
  }
  
  def getAveragePairwiseDistance(clusterMembers: Array[Long]): Double = {
    var sum = 0
    val buf1 = new Array[Byte](readLen)
    val buf2 = new Array[Byte](readLen)

    (0 until clusterMembers.size).foreach(i => {
      (i + 1 until clusterMembers.size).foreach(j => {
        // get distance between i & j
        val pos1 = clusterMembers(i)
        val pos2 = clusterMembers(j)

        genome.getSubstring(pos1, pos1 + readLen, buf1)
        genome.getSubstring(pos2, pos2 + readLen, buf2)
        
        val dist = HammingDistance.distance(buf1, buf2, readLen)
        sum += dist
      })
    })
    
    sum.toDouble / (clusterMembers.size * (clusterMembers.size - 1) / 2)
  }

  def getClusterMembers(clusterId: Long): Array[Long] = {
    val members = nonTrivialMembers.get(clusterId)
    members match {
      case None => {
        throw new IllegalArgumentException("No cluster with given ID")
        null
      }
      case Some(l) => l.toLongArray(Array[Long]())
    }
  }

  def getInterClusterDistance(clusterMembers1: Array[Long], clusterMembers2: Array[Long]): Int = {
    val consensus1 = getConsensusStr(clusterMembers1)
    val consensus2 = getConsensusStr(clusterMembers2)
    HammingDistance.distance(consensus1, consensus2, readLen)
  }
  
  def getHistogramOfInterClusterDistances: Array[Int] = {
    val histogram = Array.fill(readLen + 1)(0)
    val consensusByClusterId = scala.collection.mutable.Map[Long, Array[Byte]]()

    // find & store the consensus of each cluster
    val clusterIds = nonTrivialMembers.keySet.toArray
    var i = 0
    while (i < clusterIds.size) {
      val id = clusterIds(i)
      consensusByClusterId += ((id, getConsensusStr(getClusterMembers(id))))
      
      i += 1
    }

    // find histogram of inter-cluster distances
    i = 0
    while (i < clusterIds.size) {
      var j = i + 1
      while (j < clusterIds.size) {
        val dist = HammingDistance.distance(consensusByClusterId.get(clusterIds(i)).get, consensusByClusterId.get(clusterIds(j)).get, readLen)
        histogram(dist) += 1
        
        j += 1
      }
      i += 1
    }

    histogram
  }
  
  def getHistogramOfClusterDiameters: Array[Int] = {
    val histogram = Array.fill(readLen + 1)(0)
    val clusterIds = nonTrivialMembers.keySet.toArray

    var i = 0
    while (i < clusterIds.size) {
      val id = clusterIds(i)
      val dist = getClusterDiameter(getClusterMembers(id))
      histogram(dist) += 1
      
      i += 1
    }
    
    histogram
  }
  
  def getMinimumTotalDistance: Long = {
    // total of the sum of distances of each str to its centroid, & the sum of the distances of the centroids to the global centroid
    // means I have to know the global centroid
    var sum = 0
    
    // get global centroid
    val counts = Array.fill(readLen)(Array.fill(4)(0))
    val s = new Array[Byte](readLen)
    
    println("Finding global consensus...")
    assert(genome != null)
    for (p <- genome.pieces) {
      println("Scanning " + p.name)
      val data = p.data
      var lastN = -1 // last index where we saw an 'N'
      var key = 0L
      val mask = (1L << (readLen * 2)) - 1
      var i = 0
      val printInterval = 1000000
      var nextPrint = printInterval
      while (i < data.length) {
        val base = data(i)
        key = ((key << 2) | DNA.BASE_TO_CODE(base)) & mask
        if (base == 'N')
          lastN = i
        if (lastN <= i - readLen) {
          // find consensus
          genome.getSubstring(i, i + readLen, s)
          (0 until s.length).foreach(j => {
            counts(j)(DNA.BASE_TO_CODE(s(j))) += 1
          })
        }
        if (i == nextPrint) {
          println("Position: %d/%d".format(i, data.length))
          nextPrint += printInterval
        }
        i += 1
      }
    }
    
    val globalConsensus = Array.fill(readLen)(' '.toByte)
    (0 until readLen).foreach(i => {
      globalConsensus(i) = DNA.CODE_TO_BASE(counts(i).indexOf(counts(i).max)).toByte
    })
        
    // for each cluster
    val clusterIds = nonTrivialMembers.keySet.toArray
    var i = 0
    while (i < clusterIds.size) {
      // get consensus
      val id = clusterIds(i)
      val members = getClusterMembers(id)
      val consensus = getConsensusStr(members)
      
      // for each member
      var j = 0
      while (j < members.size) {
        // get distance to consensus
        genome.getSubstring(members(j), members(j) + readLen, s)
        val distFromMemberToConsensus = HammingDistance.distance(s, consensus, readLen)
        sum += distFromMemberToConsensus

        // get distance from consensus to global centroid
        val distFromConsensusToGlobalConsensus = HammingDistance.distance(consensus, globalConsensus, readLen)
        sum += distFromConsensusToGlobalConsensus
        
        j += 1
      }
      
      i += 1
    }
    
    // get overall minimum total distance
    sum
  }
  
  def writeClusterSizesToFile(filename: String) = {
    val fw = new java.io.FileWriter(filename)
    val bw = new java.io.BufferedWriter(fw)
    
    clusterIdsBySize.keySet.foreach(size => {
      (0 until clusterIdsBySize.get(size).get.size).foreach(i => {
        bw.write(size.toString)
        bw.newLine
      })
    })
    
    bw.close
  }
  
  def subCluster(clusterId: Long, kmeansFilename: String, n: Int) = {
    // read in the kmeans file
    val kmeans = Source.fromFile(kmeansFilename).getLines.map(_.toInt).toArray
    
    val members = nonTrivialMembers.get(clusterId)
    members match {
      case None => println("Invalid cluster ID")
      case Some(l) => {
        assert(kmeans.size == l.size)
        
        val clusterMembers = l.toLongArray(Array[Long]())
        for (k: Int <- kmeans.toSet) {
          // filter cluster members that are in subcluster i
          val subClusterMembers = (0 until clusterMembers.length).flatMap(i => 
            if (kmeans(i) == k)
              List(clusterMembers(i)) 
            else
              Nil
          ).toArray
          
          // print info about subcluster
          println("Consensus for cluster " + k + ": ")
          println(getConsensusStr(subClusterMembers))
          
          // print a random sample of n subcluster members
          val shuffledMembers = scala.util.Random.shuffle(subClusterMembers.toList)
          println("Members:")
          shuffledMembers.take(n).foreach(pos => {
            println(genome.substring(pos, pos + readLen))
          })
          
          // print cluster size
          println("Cluster size: " + subClusterMembers.size)
          
          // print cluster diameter
          println("Cluster diameter: " + getClusterDiameter(subClusterMembers))
          
          // print average pairwise distance
          println("Average pairwise distance: " + getAveragePairwiseDistance(subClusterMembers))
          
          // print histogram of diffs to consensus
          println("Histogram of diffs to consensus: " + getHistogramOfDiffsToConsensus(subClusterMembers).mkString(", "))
          
          for (k2: Int <- kmeans.toSet.diff(Set(k))) {
            val subClusterMembers2 = (0 until clusterMembers.length).flatMap(i => 
              if (kmeans(i) == k2)
                List(clusterMembers(i))
              else
                Nil
            ).toArray
            
            val interDist = getInterClusterDistance(subClusterMembers, subClusterMembers2)
            println("Distance between clusters " + k + " & " + k2 + ": " + interDist)
            println
          }
        }
      }
    }
  }
  
  def printStats = {
    (5 to 100 by 5).foreach(i => {
      println("# clusters with size at least " + i + ": " + nClustersWithSizeAtLeast(i))
    })
  }

  // does not change array in place
  // might want to make it change array in place to avoid creating an extra array if this is too slow
  def mutateDNA(str: Array[Byte], mutationRate: Double): Array[Byte] = {
    // determine how many mutations to create
    val mutatedStr = str.clone
    val numMutations = (mutationRate * readLen).toInt
    val posToMutate = scala.util.Random.shuffle(0 until readLen).take(numMutations)
    
    var i = 0
    while (i < numMutations) {
      // substitute RC for old base
      val oldChar = str(posToMutate(i))
      mutatedStr(posToMutate(i)) = DNA.COMPLEMENT(oldChar)
      
      i += 1
    }
    mutatedStr
  }

  def generateReadsFromCluster(clusterId: Long, numReads: Int, mutationRate: Double, filename: String) = {
    assert(genome != null)

    val writer = FASTQ.writer(filename)
    val members = getClusterMembers(clusterId)
    val numMembers = members.size
    val qualityStr = ("2" * readLen).getBytes
    
    var i = 0
    val str = new Array[Byte](readLen)
    
    while (i < numReads) {
      // pull out a random cluster member
      val pos = members(scala.util.Random.nextInt(numMembers))
      
      // get the genome substring
      genome.getSubstring(pos, pos + readLen, str)
      
      // mutate it
      val mutatedStr = mutateDNA(str, mutationRate)
      
      // write to file
      val idStr = Wgsim.getWgsimId("", pos, pos + readLen).getBytes
      writer.write(new Read(idStr, mutatedStr, qualityStr))
      
      i += 1
    }
    
    writer.close
  }

  def generateFileWithSubclusters(clusterId: Long, kmeansFilename: String, outputFilename: String) = {
    val kmeans = Source.fromFile(kmeansFilename).getLines.map(_.toInt).toArray
    val members = getClusterMembers(clusterId)
    
    val fw = new java.io.FileWriter(outputFilename)
    val bw = new java.io.BufferedWriter(fw)

    bw.write(List(genomeSize, kmeans.max, readLen, mergeDist).mkString(" "))
    bw.newLine
    
    for (k: Int <- kmeans.toSet) {
      // filter cluster members that are in subcluster i
      val subClusterMembers = (0 until members.length).flatMap(i => 
        if (kmeans(i) == k)
          List(members(i)) 
        else
          Nil
      ).toArray
      
      bw.write(subClusterMembers.size + " " + subClusterMembers.mkString(" "))
      bw.newLine
    }
    
    bw.close
  }
}