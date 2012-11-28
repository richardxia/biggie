package snap.apps

import snap._

import spark.SparkContext
import SparkContext._
import scala.math._
import it.unimi.dsi.fastutil.longs.LongArrayList
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import scala.collection.mutable.Set

object BruteForce {
  def main(args: Array[String]) : Unit = {
    args match {
      // Prerequisites
      case Array("index", fname, seedLen) => {
        val idxStart = System.currentTimeMillis
        val builder = new IntHashIndexBuilder(seedLen.toInt, GenomeLoader.genome.totalSize)
        GenomeLoader.genome.addToIndex(builder)
        val index = builder.build()
        printf("Indexing took %.3fs%n", (System.currentTimeMillis - idxStart) / 1000.0)

        val fos = new java.io.FileOutputStream(fname)
        val oos = new java.io.ObjectOutputStream(fos)
        oos.writeObject(index)
        oos.close
      }
      case Array("snapAlign", fastqFile, numSeeds, seedLen, seedsToTry, maxDist, confDiff, maxHits, outputFname) => {
        val reads = FASTQ.read(fastqFile, true)

        println("Indexing genome...")
        val idxStart = System.currentTimeMillis
        val builder = new IntHashIndexBuilder(seedLen.toInt, GenomeLoader.genome.totalSize)
        GenomeLoader.genome.addToIndex(builder)
        val index = builder.build()
        printf("Indexing took %.3fs%n", (System.currentTimeMillis - idxStart) / 1000.0)

        println("Aligning reads with SNAP...")
        val snapResults = BruteForceAligner.snapAlign(index, reads, numSeeds.toInt, seedsToTry.toInt, maxDist.toInt, confDiff.toInt, maxHits.toInt)

        println("Saving align results to file...")
        val fos = new java.io.FileOutputStream(outputFname)
        val oos = new java.io.ObjectOutputStream(fos)
        oos.writeObject(snapResults)
        oos.close
      }
      case Array("snapAlignPaired", fastqFile1, fastqFile2, numSeeds, seedLen, seedsToTry, maxDist, confDiff, maxHits, outputFname1, outputFname2) => {
        println("Indexing genome...")
        val idxStart = System.currentTimeMillis
        val builder = new IntHashIndexBuilder(seedLen.toInt, GenomeLoader.genome.totalSize)
        GenomeLoader.genome.addToIndex(builder)
        val index = builder.build()
        printf("Indexing took %.3fs%n", (System.currentTimeMillis - idxStart) / 1000.0)

        var reads = FASTQ.read(fastqFile1, true)
        println("Aligning reads with SNAP...")
        var snapResults = BruteForceAligner.snapAlign(index, reads, numSeeds.toInt, seedsToTry.toInt, maxDist.toInt, confDiff.toInt, maxHits.toInt)

        println("Saving align results to file...")
        var fos = new java.io.FileOutputStream(outputFname1)
        var oos = new java.io.ObjectOutputStream(fos)
        oos.writeObject(snapResults)
        oos.close

        reads = FASTQ.read(fastqFile2, true)
        println("Aligning reads with SNAP...")
        snapResults = BruteForceAligner.snapAlign(index, reads, numSeeds.toInt, seedsToTry.toInt, maxDist.toInt, confDiff.toInt, maxHits.toInt)

        println("Saving align results to file...")
        fos = new java.io.FileOutputStream(outputFname2)
        oos = new java.io.ObjectOutputStream(fos)
        oos.writeObject(snapResults)
        oos.close

      }
      case Array("snapSamToAlignResult", snapSam) => {
        val snapSamEntries = AlignerJudge.loadSam(snapSam).map(s => new RichSamEntry(s))
        
        // load snap results (sam format)
        val snapSamEntries1 = snapSamEntries.filter(e => e.ID.split("/")(1).startsWith("1"))
        val snapSamEntries2 = snapSamEntries.filter(e => e.ID.split("/")(1).startsWith("2"))
        
        // create align results from sam entries
        val snapResults1 = snapSamEntries1.map(e => samEntryToAlignResult(e)).toArray
        val snapResults2 = snapSamEntries2.map(e => samEntryToAlignResult(e)).toArray
        
        // write out to file
        var fos = new java.io.FileOutputStream(snapSam + "1.snap")
        var oos = new java.io.ObjectOutputStream(fos)
        oos.writeObject(snapResults1)
        oos.close
        
        fos = new java.io.FileOutputStream(snapSam + "2.snap")
        oos = new java.io.ObjectOutputStream(fos)
        oos.writeObject(snapResults2)
        oos.close
      }
      case Array("checkSnapEditDistances", snapSam) => {
        val snapSamEntries = AlignerJudge.loadSam(snapSam).map(s => new RichSamEntry(s))
        
        // load snap results (sam format)
        val snapSamEntries1 = snapSamEntries.filter(e => e.ID.split("/")(1).startsWith("1"))
        val snapSamEntries2 = snapSamEntries.filter(e => e.ID.split("/")(1).startsWith("2"))

        (0 until snapSamEntries1.length).foreach(i => {
          val dist1 = snapSamEntries1(i).readEditDist
          val dist2 = snapSamEntries2(i).readEditDist
          val pairDist = dist1 + dist2
          if (pairDist > 20)
            println(i + ": " + pairDist)
        })
      }
      case Array("removeDuplicateReads", fastqFile) => {
        // load reads
        val reads = FASTQ.read(fastqFile, true /* keep quality */)
        
        val idsSeen = scala.collection.mutable.Set[String]()
        var uniqueReads = List[Read]()
        
        reads.foreach(r => {
          // get id
          val id = r.idStr
          
          // remove erroneous suffix
          val i = id.indexOf("/")
          val truncatedId = id.substring(0, i + 2 /* to include slash & char after slash */)
          
          // see if this id has been seen yet
          // if yes, skip
          // if not, add this read to the list
          if (!idsSeen.contains(truncatedId)) {
            idsSeen.add(truncatedId)
            val newR = new Read(truncatedId.getBytes, r.data, r.quality)
            uniqueReads = newR :: uniqueReads
          }
        })
        uniqueReads = uniqueReads.reverse
        
        // print unique reads to file (with new fname)
        val w = FASTQ.writer(fastqFile + ".unique")
        uniqueReads.foreach(r => {
          w.write(r)
        })
        
        w.close
      }
      case Array("removeDuplicateReadsFromBfaResults", accessKey, secretKey, bfaBucket) => {
        val sc = new SparkContext("local[4]", "bruteForceAligner")
        val readsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)](s3Path(accessKey, secretKey, None, bfaBucket)).collect
        
        val idsSeen1 = scala.collection.mutable.Set[String]()
        var uniqueResults = List[((Read, Read), AlignResult)]()
        
        readsWithBfaResults.foreach(resTuple => {
          val (readPair, bfaRes) = resTuple
    	    val (read1, read2) = readPair
    	    
          // get id of 1st read
          val id1 = read1.idStr
          
          // remove erroneous suffix
          val i1 = id1.indexOf("/")
          val truncatedId1 = id1.substring(0, i1 + 2 /* to include slash & char after slash */)
          
          // see if this id has been seen yet
          // if yes, skip
          // if not, add this tuple to the list (with both reads' ids fixed)
          if (!idsSeen1.contains(truncatedId1)) {
            idsSeen1.add(truncatedId1)
            val newR1 = new Read(truncatedId1.getBytes, read1.data, read1.quality)
            
            val id2 = read2.idStr
            val i2 = id2.indexOf("/")
            val truncatedId2 = id2.substring(0, i2 + 2)
            val newR2 = new Read(truncatedId2.getBytes, read2.data, read2.quality)
            
            uniqueResults = ((newR1, newR2), bfaRes) :: uniqueResults
          }
          uniqueResults = uniqueResults.reverse
        })
        
        // save to s3
        val path = s3Path(accessKey, secretKey, None, bfaBucket + "Unique")
        sc.parallelize(uniqueResults).saveAsObjectFile(path)
        
        println("# unique results: " + uniqueResults.length)
        println("# total results: " + readsWithBfaResults.length)
      }

      // Brute Force
      case Array("BFA", master, fastqFile, numSeeds, seedLen, seedsToTry, maxDist, confDiff, maxHits) => {
        println("Loading reads...")
        val reads = FASTQ.read(fastqFile, true)
        println("Indexing genome...")
        val idxStart = System.currentTimeMillis
        val builder = new IntHashIndexBuilder(seedLen.toInt, GenomeLoader.genome.totalSize)
        GenomeLoader.genome.addToIndex(builder)
        val index = builder.build()
        printf("Indexing took %.3fs%n", (System.currentTimeMillis - idxStart) / 1000.0)

        println("Aligning reads with SNAP...")
        val snapResults = BruteForceAligner.snapAlign(index, reads, numSeeds.toInt, seedsToTry.toInt, maxDist.toInt, confDiff.toInt, maxHits.toInt)

        val bfaResults = new Array[AlignResult](reads.size)

        println("Aligning reads with BFA...")
        val alignStart = System.currentTimeMillis

        var i = 0
        while (i < reads.size) {
          println("Aligning read " + i + "...")

          // if snap returned single hit or multi hit, use snap's edit distance instead of maxDist
          val (initDist, snapRC) = 
            snapResults(i) match {
              case RichSingleHit(pos, isRC, editDistance) => (editDistance, isRC)
              case RichMultipleHits(bestPos, bestScore, isRC, secondBestPos, secondBestScore) => (bestScore, isRC) // maybe secondBestScore?
              case _ => (maxDist.toInt, false)
            }

          val start = System.currentTimeMillis
          bfaResults(i) = BruteForceAligner.align(reads(i).data, initDist, maxDist.toInt, confDiff.toInt, snapRC)
          printf("Aligning read took %.3fs%n", (System.currentTimeMillis - start) / 1000.0)

          i += 1
        }
        printf("Aligning took %.3fs%n", (System.currentTimeMillis - alignStart) / 1000.0)

        // score bfa alignments against ground truth from simulator
        BruteForceAligner.scoreBFA(reads, bfaResults, maxDist.toInt, Some(snapResults))
      }
      case Array("BFA", fastqFile, maxDist, confDiff, snapResultsFname) => {
        val reads = FASTQ.read(fastqFile, true)

        // requires that you've already aligned reads via SNAP
        // do so via "snapAlign" option
        val fileIn = new java.io.FileInputStream(snapResultsFname)
        val in = new java.io.ObjectInputStream(fileIn)
        val snapResults = in.readObject.asInstanceOf[Array[AlignResult]]

        //val readsWithSnapResults = reads.zip(snapResults) // for scoring
        val bfaResults = new Array[AlignResult](reads.size)

        println("Aligning reads with BFA...")
        val alignStart = System.currentTimeMillis

        var i = 0
        while (i < reads.size) {
          println("Aligning read " + i + "...")

          // if snap returned single hit or multi hit, use snap's edit distance instead of maxDist
          val (initDist, snapRC) = 
            snapResults(i) match {
              case RichSingleHit(pos, isRC, editDistance) => (editDistance, isRC)
              case RichMultipleHits(bestPos, bestScore, isRC, secondBestPos, secondBestScore) => (bestScore, isRC) // maybe secondBestScore?
              case _ => (maxDist.toInt, false)
            }

          val start = System.currentTimeMillis
          bfaResults(i) = BruteForceAligner.align(reads(i).data, initDist, maxDist.toInt, confDiff.toInt, snapRC)
          printf("Aligning read took %.3fs%n", (System.currentTimeMillis - start) / 1000.0)

          i += 1
        }
        printf("Aligning took %.3fs%n", (System.currentTimeMillis - alignStart) / 1000.0)

      }
      case Array("parallelBFA", master, fastqFile, numSeeds, seedLen, seedsToTry, maxDist, confDiff, maxHits) => {
        val reads = FASTQ.read(fastqFile, true)
        println("Indexing genome...")
        val idxStart = System.currentTimeMillis
        val builder = new IntHashIndexBuilder(seedLen.toInt, GenomeLoader.genome.totalSize)
        GenomeLoader.genome.addToIndex(builder)
        val index = builder.build()
        printf("Indexing took %.3fs%n", (System.currentTimeMillis - idxStart) / 1000.0)

        println("Aligning reads with SNAP...")
        val snapResults = BruteForceAligner.snapAlign(index, reads, numSeeds.toInt, seedsToTry.toInt, maxDist.toInt, confDiff.toInt, maxHits.toInt)

        val readsWithSnapResults = reads.zip(snapResults)

        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")

        val alignStart = System.currentTimeMillis
        BruteForceAligner.parallelAlign(sc, master, readsWithSnapResults, maxDist.toInt, confDiff.toInt)
        val bfaTimeInSec = (System.currentTimeMillis - alignStart) / 1000.0
        println("Aligning with BFA (using Spark) took " + bfaTimeInSec + "s for " + reads.length + " reads.")
        println("Speed: " + reads.length/bfaTimeInSec + " reads/s.")

        val readsWithBfaResults = sc.objectFile[(Read, AlignResult)]("hdfs://" + master + ":9000/bfaResults").collect
        BruteForceAligner.parallelScoreBFA(readsWithBfaResults, maxDist.toInt)

        System.exit(0)
      }
      case Array("parallelBFA", master, accessKey, secretKey, fastqFile, maxDist, confDiff, snapResultsFname, chooseOutputDest, outputPrefix) => {
        val reads = FASTQ.read(fastqFile, true)

        // requires that you've already aligned reads via SNAP
        // do so via "snapAlign" option
        val fileIn = new java.io.FileInputStream(snapResultsFname)
        val in = new java.io.ObjectInputStream(fileIn)
        val snapResults = in.readObject.asInstanceOf[Array[AlignResult]]

        // make sure I've gotten the right snap, read files
        if (reads.length != snapResults.length)
          throw new IllegalArgumentException("Mismatch between # reads (" + reads.length + ") and # snap align results (" + snapResults.length + ").")

        val readsWithSnapResults = reads.zip(snapResults)
        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")
        var outputDest = ""

        chooseOutputDest match {
          case "hdfs" => outputDest = "hdfs://" + master + ":9000/"
          case "s3" => outputDest = "s3n://" + accessKey + ":" + secretKey + "@hashingGenome/"
          case _ => {
            println("Incorrect output destination.  Should be hdfs or s3.")
            System.exit(0)
          }
        }

        val alignStart = System.currentTimeMillis
        BruteForceAligner.parallelAlign(sc, outputDest, readsWithSnapResults, maxDist.toInt, confDiff.toInt, outputPrefix)
        val bfaTimeInSec = (System.currentTimeMillis - alignStart) / 1000.0
        println("Aligning with BFA (using Spark) took " + bfaTimeInSec + "s for " + reads.length + " reads.")
        println("Speed: " + reads.length/bfaTimeInSec + " reads/s")

        val readsWithBfaResults = sc.objectFile[(Read, AlignResult)](outputDest + outputPrefix).collect

        // write out reads with wgsim-like id strings
        BruteForceAligner.writeToFastq(readsWithBfaResults, outputPrefix + ".fq")

        System.exit(0)
      }
      case Array("parallelPairedBFA", master, accessKey, secretKey, fastqFile1, fastqFile2, snapFile1, snapFile2, maxDist, confDiff, separationDist, chooseOutputDest, outputPrefix) => {
        // load raw reads
        println("Loading raw reads...")
        val reads1 = FASTQ.read(fastqFile1, true)
        val reads2 = FASTQ.read(fastqFile2, true)

        // load SNAP results
        println("Loading SNAP results...")
        var fileIn = new java.io.FileInputStream(snapFile1)
        var in = new java.io.ObjectInputStream(fileIn)
        val snapResults1 = in.readObject.asInstanceOf[Array[AlignResult]]
        val readsWithSnapResults1 = reads1.zip(snapResults1)

        fileIn = new java.io.FileInputStream(snapFile2)
        in = new java.io.ObjectInputStream(fileIn)
        val snapResults2 = in.readObject.asInstanceOf[Array[AlignResult]]
        val readsWithSnapResults2 = reads2.zip(snapResults2)

        val readsWithSnapResults = readsWithSnapResults1.zip(readsWithSnapResults2)

        // set output destination
        var outputDest = ""
        chooseOutputDest match {
          case "hdfs" => outputDest = "hdfs://" + master + ":9000/"
          case "s3" => outputDest = "s3n://" + accessKey + ":" + secretKey + "@hashingGenome/"
          case _ => {
            println("Incorrect output destination.  Should be hdfs or s3.")
            System.exit(0)
          }
        }

        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")

        val alignStart = System.currentTimeMillis
        BruteForceAligner.parallelPairedAlign(sc, outputDest, outputPrefix, readsWithSnapResults, maxDist.toInt, confDiff.toInt, separationDist.toInt)

        val bfaTimeInSec = (System.currentTimeMillis - alignStart) / 1000.0
        println("Aligning with BFA (using Spark) took " + bfaTimeInSec + "s for " + reads1.length + " reads.")
        println("Speed: " + reads1.length/bfaTimeInSec + " reads/s")

        // print out reads (with wgsim-style string)
        val readsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)](outputDest + outputPrefix).collect
        BruteForceAligner.writePairToFastq(readsWithBfaResults, outputPrefix)

        System.exit(0)
      }
      case Array("simplePairedBFA", master, sparkMaster, accessKey, secretKey, fastq1, fastq2, snap1, snap2, maxDist, confDiff, separationDist, outputDest, outputBucket) => {
        val params = new PairedBruteForceParams(maxDist.toInt, confDiff.toInt, separationDist.toInt)
        val c = new PairedBruteForceContext(master, sparkMaster, accessKey, secretKey, fastq1, fastq2, snap1, snap2, outputDest, outputBucket, params)
        val sc = new SparkContext(c.sparkDest, "simpleBruteForceAligner")

        simpleBfaAndSaveOutput(sc, c, params)
      }

      // Postprocessing
      case Array("bfaResultsToFastq", master, accessKey, secretKey, bfaResultsFname, outputFilename) => {
        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")
        val readsWithBfaResults = sc.objectFile[(Read, AlignResult)]("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + bfaResultsFname).collect

        // write out reads with wgsim-like id strings
        BruteForceAligner.writeToFastq(readsWithBfaResults, outputFilename)
      }
      case Array("pairedBfaResultsToFastq", accessKey, secretKey, bfaResultsBucket, outputPrefix, confDiff) => {
        val sc = new SparkContext("local[8]", "bruteForceAligner")
        val readsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)](s3Path(accessKey, secretKey, None, bfaResultsBucket)).collect
        BruteForceAligner.writePairToFastq(readsWithBfaResults, outputPrefix, confDiff.toInt)
      }
      case Array("scoreBFA", master, accessKey, secretKey, outputPrefix, maxDist) => {
        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")

        val readsWithBfaResults = sc.objectFile[(Read, AlignResult)]("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + outputPrefix).collect
        BruteForceAligner.parallelScoreBFA(readsWithBfaResults, maxDist.toInt)

        System.exit(0)
      }
      case Array("scorePairedBFA", master, bfaResultsFname, maxDist, snapFile1, snapFile2) => {
        // load reads & bfa results
        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")
        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)]("hdfs://" + master + ":9000/" + bfaResultsFname).collect
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip

        // load snap results
        println("Loading SNAP results...")
        var fileIn = new java.io.FileInputStream(snapFile1)
        var in = new java.io.ObjectInputStream(fileIn)
        val snapResults1 = in.readObject.asInstanceOf[Array[AlignResult]]

        fileIn = new java.io.FileInputStream(snapFile2)
        in = new java.io.ObjectInputStream(fileIn)
        val snapResults2 = in.readObject.asInstanceOf[Array[AlignResult]]

        // score bfa
        BruteForceAligner.scorePairedBFA(reads1, reads2, bfaResults, maxDist.toInt, Some(snapResults1), Some(snapResults2))
      }
      case Array("scorePairedBFA", master, bfaResultsFname, maxDist) => {
        // load reads & bfa results
        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")
        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)]("hdfs://" + master + ":9000/" + bfaResultsFname).collect
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip

        // score bfa
        BruteForceAligner.scorePairedBFA(reads1, reads2, bfaResults, maxDist.toInt, None, None)
      }
      case Array("getPercentAligned", fname) => {
        // read file
        val fReader = new java.io.FileReader(fname)
        val bReader = new java.io.BufferedReader(fReader)

        var line = bReader.readLine
        var numConfident = 0
        var numReads = 0

        while (line != null) {
          // skip lines in header
          if (!line.startsWith("@")) {
            numReads += 1

            val fields = line.split("\t")

            // pull out quality
            val quality = fields(4).toInt // 5th field (http://samtools.sourceforge.net/SAM1.pdf)

            if (quality > 10 && quality != 255)
              numConfident += 1
          }

          line = bReader.readLine
        }

        bReader.close

        println("# Reads: " + numReads)
        println("# Confident: " + numConfident)
        println("% aligned: " + (numConfident.toDouble / numReads) * 100)
      }
      case Array("compareAlignerResultsToBfaResults", accessKey, secretKey, alignerResultsSam, bfaBucket, maxDist, confDiff, minSeparation, maxSeparation, confidenceThreshold) => {
        // recall that the input to the aligner is just the reads i found a single hit for with the bfa
        
        // load aligner results
        val alignerEntries = AlignerJudge.loadSam(alignerResultsSam).map(samLine => new RichSamEntry(samLine)).sortWith(_.illuminaID < _.illuminaID)
        assert(alignerEntries.length % 2 == 0)

	      //val sortedAlignerEntries = alignerEntries.sortWith(_.illuminaID < _.illuminaID)

        // split into results for first, second sets of reads
        var alignerEntries1 = List[RichSamEntry]()
        var alignerEntries2 = List[RichSamEntry]()
        (0 until alignerEntries.length).foreach(i => {
          if (i % 2 == 0) alignerEntries1 = alignerEntries(i) :: alignerEntries1
          else alignerEntries2 = alignerEntries(i) :: alignerEntries2
        })
        alignerEntries1 = alignerEntries1.reverse
        alignerEntries2 = alignerEntries2.reverse

        assert(alignerEntries1.length == alignerEntries2.length) // make sure I split up the align results correctly
        val alignerResults1 = alignerEntries1.map(e => samEntryToAlignResult(e))
        val alignerResults2 = alignerEntries2.map(e => samEntryToAlignResult(e))
        
        // load reads & bfa results
        val readPairsWithBfaResults = AlignerJudge.loadReadPairsAndBfaResults(None, accessKey, secretKey, None, bfaBucket)
        /*
        val singleHitReadPairsWithBfaResults = readPairsWithBfaResults.filter(t => t._2.isInstanceOf[RichPairSingleHit] || {
          if (t._2.isInstanceOf[RichPairMultiHit]) {
            val RichPairMultiHit(bestPos1, bestPos2, bestScore, secondBestPos1, secondBestPos2, secondBestScore, isRC) = t._2.asInstanceOf[RichPairMultiHit]
            bestScore < secondBestScore
          } else false
        })
        */
        val singleHitReadPairsWithBfaResults = readPairsWithBfaResults.flatMap(t => {
          t._2 match {
            case RichPairSingleHit(_, _, _, _) => List(t)
            case RichPairMultiHit(bestPos1, bestPos2, bestScore, secondBestPos1, secondBestPos2, secondBestScore, isRC) => {
              //if (bestScore < secondBestScore)
              if (bestScore + confDiff.toInt <= secondBestScore)
                List((t._1, RichPairSingleHit(bestPos1, bestPos2, bestScore, isRC)))
              else
                Nil
            }
            case _ => Nil
          }
        })

      	val sortedHitReadPairsWithBfaResults = singleHitReadPairsWithBfaResults.sortWith((s, t) => {
      	    val (sReadPair, sBfaRes) = s
      	    val (tReadPair, tBfaRes) = t
      	    val (sRead1, sRead2) = sReadPair
      	    val (tRead1, tRead2) = tReadPair
      	    sRead1.idStr < tRead1.idStr
      	})

        /*
      	(0 until 5).foreach(i => {
      	   println(i)
      	   println(alignerEntries(2 * i).fields.mkString("\t"))
      	   println(alignerEntries(2 * i + 1).fields.mkString("\t"))
      	   println(sortedHitReadPairsWithBfaResults(i)._1._1.idStr)
      	   println(sortedHitReadPairsWithBfaResults(i)._1._2.idStr)
      	   println(sortedHitReadPairsWithBfaResults(i)._2)
      	})
	      */

        //val (readPairs, bfaResults) = singleHitReadPairsWithBfaResults.unzip
	      val (readPairs, bfaResults) = sortedHitReadPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip

        assert(alignerResults1.length == bfaResults.length) // make sure I filtered out NotFound/MultiHit correctly

        // for each pair of results
        // for each read in the pair
        // check whether the aligner came up with the same answer as the brute force
        // if not, check whether i need to verify this pair 
        // (sometimes want to rule out b/c of mapping reads to different chromosomes, or out-of-range separation distance, whether too small or too large)
        // if need to verify, check whether the aligner or the brute force made the mistake
        var nTotal = alignerEntries.length
        var nConfident = 0
        var nError = 0

      	var error1 = false
      	var error2 = false        
        (0 until bfaResults.length).foreach(i => {
          if (includePair(bfaResults(i).asInstanceOf[RichPairSingleHit], minSeparation.toInt, maxSeparation.toInt)) {
            // read 1
            val read1 = reads1(i)
            //if (alignerResults1(i).isInstanceOf[RichSingleHit])
            if (alignerEntries1(i).confident(confidenceThreshold.toInt)) {
                nConfident += 1

                //if (!AlignmentVerifier.checkForMatch(bfaResults(i).asInstanceOf[RichPairSingleHit], alignerResults1(i).asInstanceOf[RichSingleHit], true, maxDist.toInt)) {
                if (!AlignmentVerifier.checkForMatch(bfaResults(i).asInstanceOf[RichPairSingleHit], alignerResults1(i).asInstanceOf[RichSingleHit], true, maxDist.toInt, alignerEntries1(i).numPrefixBasesClipped)) {
                  nError += 1
          	      error1 = true
          	    }
              }
              
            // read 2get
            val read2 = reads2(i)

            //if (alignerResults2(i).isInstanceOf[RichSingleHit])
            if (alignerEntries2(i).confident(confidenceThreshold.toInt)) {
              nConfident += 1

              //if (!AlignmentVerifier.checkForMatch(bfaResults(i).asInstanceOf[RichPairSingleHit], alignerResults2(i).asInstanceOf[RichSingleHit], false, maxDist.toInt)) {
              if (!AlignmentVerifier.checkForMatch(bfaResults(i).asInstanceOf[RichPairSingleHit], alignerResults2(i).asInstanceOf[RichSingleHit], false, maxDist.toInt, alignerEntries2(i).numPrefixBasesClipped)) {  
                nError += 1
        	      error2 = true
        	    }
            }

      	    if (error1 || error2) {
      	      // error info
      	      println(i)
      	      println("read1:")
      	      println(read1.dataStr)
      	      println(read1.qualityStr)
              println(alignerEntries1(i).cigarStr)
      	      println("read2:")
      	      println(read2.dataStr)
      	      println(read2.qualityStr)
              println(alignerEntries2(i).cigarStr)
      	      println(bfaResults(i))
      	      println(alignerResults1(i))
      	      if (error1) println("Error at " + i + ".1")
      	      println(alignerResults2(i))
      	      if (error2) println("Error at " + i + ".2")
      	      println
      	      
      	      // verify
      	      if (AlignmentVerifier.needToVerify(alignerResults1(i), alignerResults2(i), maxSeparation.toInt)) {
      	        AlignmentVerifier.verifyAlignerError(read1.data, read2.data, bfaResults(i).asInstanceOf[RichPairSingleHit], 
      	          RichPairSingleHitHelper.getPairHit(alignerResults1(i).asInstanceOf[RichSingleHit], alignerResults2(i).asInstanceOf[RichSingleHit]), maxDist.toInt)
      	      }

      	      // reset
      	      error1 = false
      	      error2 = false
      	    }
      	  }
        })
        
        println("% Confident: " + (nConfident.toDouble / nTotal) * 100)
      	println("# Error: " + nError)
      	println("# Confident: " + nConfident)
      	println("# Total: " + nTotal)
        println("% Error: " + (nError.toDouble / nConfident) * 100)
      }
      case Array("findMissingReads", alignerResultsSam, bfaResultsFastq) => {
        // load aligner results
        val alignerSamEntries = AlignerJudge.loadSam(alignerResultsSam).map(samLine => new RichSamEntry(samLine))
        val alignerIds = alignerSamEntries.map(e => e.ID)
        
        // load bfa results
        val reads = FASTQ.read(bfaResultsFastq)
        val readIds = reads.map(read => read.idStr)
        
        // for each read in reads1
        // get id
        // see if aligner results has an entry for that id
        reads.foreach(read => {
          if (!alignerIds.contains(read.idStr.substring(1)))
            println(read.idStr)
        })
      }
      case Array("wgsimScore", alignerResultsSam, maxDist) => {
        val alignerEntries = AlignerJudge.loadSam(alignerResultsSam).map(samLine => new RichSamEntry(samLine))

        val total = alignerEntries.length
	      var nConfident = alignerEntries.count(e => e.confident())
        var nError = 0

        alignerEntries.foreach(e => {
          if (e.confident() && !Wgsim.isCorrect(e.ID, e.readAbsPos, maxDist.toInt, false)) nError += 1
        })
        
        println("Error rate: " + (nError.toDouble / nConfident * 100.0) + "%")
	      println("% Confident: " + (nConfident.toDouble / total * 100.0) + "%")
      }
      case Array("wgsimScoreAndVerify", alignerResultsSam, maxDist, minSeparation, maxSeparation) => {
        val alignerEntries = AlignerJudge.loadSam(alignerResultsSam).map(samLine => new RichSamEntry(samLine))

        var alignerEntries1 = List[RichSamEntry]()
        var alignerEntries2 = List[RichSamEntry]()
        (0 until alignerEntries.length).foreach(i => {
          if (i % 2 == 0) alignerEntries1 = alignerEntries(i) :: alignerEntries1
          else alignerEntries2 = alignerEntries(i) :: alignerEntries2
        })
        alignerEntries1 = alignerEntries1.reverse
        alignerEntries2 = alignerEntries2.reverse

        assert(alignerEntries1.length == alignerEntries2.length) // make sure I split up the align results correctly
        val alignerResults1 = alignerEntries1.map(e => samEntryToAlignResult(e))
        val alignerResults2 = alignerEntries2.map(e => samEntryToAlignResult(e))

        val total = alignerEntries.length
	      var nConfident = alignerEntries.count(e => e.confident())
        var nError = 0

        //alignerEntries.foreach(e => {
        (0 until alignerEntries1.length).foreach(i => {
          // check corresponding entries
          if (alignerEntries1(i).confident() && !Wgsim.isCorrect(alignerEntries1(i).ID, alignerEntries1(i).readAbsPos, maxDist.toInt, false)) nError += 1
          if (alignerEntries2(i).confident() && !Wgsim.isCorrect(alignerEntries2(i).ID, alignerEntries2(i).readAbsPos, maxDist.toInt, false)) nError += 1
          
          // verify
          if (AlignmentVerifier.needToVerify(alignerResults1(i), alignerResults2(i), maxSeparation.toInt)) {
            // pull out the bfa positions
            val idStr = alignerEntries1(i).ID
            val wgsimId = idStr.split(":")(0)
            val bfaPos1 = GenomeLoader.genome.getAbsPos(Wgsim.getChr(wgsimId), Wgsim.getLowEnd(wgsimId))
            val bfaPos2 = GenomeLoader.genome.getAbsPos(Wgsim.getChr(wgsimId), Wgsim.getHighEnd(wgsimId))
            val bfaRes = RichPairSingleHit(bfaPos1, bfaPos2, 0, false /* last 2 fields will be ignored */)
            
  	        AlignmentVerifier.verifyAlignerError(alignerEntries1(i).read.data, alignerEntries2(i).read.data, bfaRes, 
  	          RichPairSingleHitHelper.getPairHit(alignerResults1(i).asInstanceOf[RichSingleHit], alignerResults2(i).asInstanceOf[RichSingleHit]), maxDist.toInt)
  	      }
        })
        
        println("Error rate: " + (nError.toDouble / nConfident * 100.0) + "%")
	      println("% Confident: " + (nConfident.toDouble / total * 100.0) + "%")
      }
      case Array("findIntersection", accessKey, secretKey) => {
        val sc = new SparkContext("local[4]", "bruteForceAligner")
        val readsWithBfaResults1 = sc.objectFile[((Read, Read), AlignResult)](s3Path(accessKey, secretKey, None, "simpleBfa500sample")).collect
        val readsWithBfaResults2 = sc.objectFile[((Read, Read), AlignResult)](s3Path(accessKey, secretKey, None, "simpleBfa800sample")).collect

        val (readPairs1, bfaResults1) = readsWithBfaResults1.unzip
        val (readPairs2, bfaResults2) = readsWithBfaResults2.unzip
        
        // find intersect size
        val readIds1 = readPairs1.map(t => t._1.idStr).toSet
        val readIds2 = readPairs2.map(t => t._1.idStr).toSet
        
        println("Intersect size: " + readIds1.intersect(readIds2).size)
      }
      case Array("getUniqueReads", set1A, set1B, set2A, set2B, outA, outB) => {
        // Note:  |set1| < |set2|

        // load reads from set1
        var set1AReads = FASTQ.read(set1A, true)
        var set1BReads = FASTQ.read(set1B, true)
        
        // load reads from set2
        var set2AReads = FASTQ.read(set2A, true)
        var set2BReads = FASTQ.read(set2B, true)
        
        // print out statistics about intersection
        val set1ReadIds = set1AReads.map(_.idStr).toSet
        val set2ReadIds = set2AReads.map(_.idStr).toSet
        println("Intersect size: " + set1ReadIds.intersect(set2ReadIds).size)
        
        // print out reads in set2 not in set1 to file
        val writerA = new FASTQ.Writer(outA)
        val writerB = new FASTQ.Writer(outB)
        
        (0 until set2AReads.length).foreach(i => {
          if (!set1ReadIds.contains(set2AReads(i).idStr)) {
            writerA.write(set2AReads(i))
            writerB.write(set2BReads(i))
          }
        })
        
        writerA.close
        writerB.close
      }
      case Array("mergeBfaResults", accessKey, secretKey, bfaBucket1, bfaBucket2, outputBucket) => {
        val baseBucket = "hashingGenome"
        val s3Path = "s3n://" + accessKey + ":" + secretKey + "@" + baseBucket + "/"
        val sc = new SparkContext("local[4]", "bruteForceAligner")

        // load results from 1st bucket
        val readsWithBfaResults1 = sc.objectFile[((Read, Read), AlignResult)](s3Path + bfaBucket1).collect
        
        // load results from 2nd bucket
        val readsWithBfaResults2 = sc.objectFile[((Read, Read), AlignResult)](s3Path + bfaBucket2).collect
        
        // merge results
        val mergedReadsWithBfaResults = readsWithBfaResults1 ++ readsWithBfaResults2
        
        // save results to s3
        sc.parallelize(mergedReadsWithBfaResults).saveAsObjectFile(s3Path + outputBucket)
      }
      case Array("countUnique", fname) => {
        // load reads
        val reads = FASTQ.read(fname, true /* keep quality */)
        
        // map to ID
        val ids = reads.map(r => r.idStr)
        
        // strip suffix
        val idPrefixes = ids.map(id => {
          val entries = id.split("/")
          entries(0)
        })
        
        // count unique
        println("# total:  " + reads.length)
        println("# unique: " + idPrefixes.toSet.size)
      }

      // Tests
      case Array("testOutput", master, outputPrefix) => {
        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")

        val readsWithBfaResults = sc.objectFile[(Read, AlignResult)]("hdfs://" + master + ":9000/" + outputPrefix).collect

        // write out reads with wgsim-like id strings
        BruteForceAligner.writeToFastq(readsWithBfaResults, outputPrefix + ".fq")

        System.exit(0)
      }
      case Array("testPairedBFA", accessKey, secretKey, snapFname1, snapFname2, bfaResultsFname, readNum) => {
        // load reads
        val sc = new SparkContext("local[8]", "bruteForceAligner")
        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)]("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + bfaResultsFname).collect
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip

        // load snap results
        var fileIn = new java.io.FileInputStream(snapFname1)
        var in = new java.io.ObjectInputStream(fileIn)
        val snapResults1 = in.readObject.asInstanceOf[Array[AlignResult]]

        fileIn = new java.io.FileInputStream(snapFname2)
        in = new java.io.ObjectInputStream(fileIn)
        val snapResults2 = in.readObject.asInstanceOf[Array[AlignResult]]

        // align with paired BFA
        val read1 = reads1(readNum.toInt)
        val read2 = reads2(readNum.toInt)

        val bfaRes = BruteForceAligner.pairedAlign(read1.data, read2.data, 10, 10, 1, 1000, false)
        println("res: " + bfaRes)
      }
      case Array("testPairedBFAOnErrors", master, accessKey, secretKey, bfaResultsFname, snapFname1, snapFname2, outputFname) => {
        // load reads
        val sc = new SparkContext("local[8]", "bruteForceAligner")
        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)]("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + bfaResultsFname).collect
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip

        // load snap results
        var fileIn = new java.io.FileInputStream(snapFname1)
        var in = new java.io.ObjectInputStream(fileIn)
        val snapResults1 = in.readObject.asInstanceOf[Array[AlignResult]]

        fileIn = new java.io.FileInputStream(snapFname2)
        in = new java.io.ObjectInputStream(fileIn)
        val snapResults2 = in.readObject.asInstanceOf[Array[AlignResult]]        

        // zip reads & snap results
        val readsWithSnapResults1 = reads1.zip(snapResults1)
        val readsWithSnapResults2 = reads2.zip(snapResults2)
        val readsWithSnapResults = readsWithSnapResults1.zip(readsWithSnapResults2)

        // pull out reads that bfa made errors on
        val errorReads = List(9, 21, 55, 58, 68, 87, 95)

        val readSubset = new Array[((Read, AlignResult), (Read, AlignResult))](errorReads.size)
        var numReads = 0
        errorReads.foreach(r => {
          readSubset(numReads) = readsWithSnapResults(r)
          numReads += 1
        })

        // align read subset
        val alignStart = System.currentTimeMillis
        BruteForceAligner.parallelPairedAlign(sc, "s3n://" + accessKey + ":" + secretKey + "@hashingGenome/", outputFname, readSubset, 10, 1, 1000)

        // if it works, then print
        val readsWithBfaResults = sc.objectFile[(Read, AlignResult)]("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + outputFname).collect
        readsWithBfaResults.foreach(println)

        System.exit(0)
      }
      case Array("testVerifyPairedHit", snapFname1, snapFname2) => {
        // load snap results
        var fileIn = new java.io.FileInputStream(snapFname1)
        var in = new java.io.ObjectInputStream(fileIn)
        val snapResults1 = in.readObject.asInstanceOf[Array[AlignResult]]

        fileIn = new java.io.FileInputStream(snapFname2)
        in = new java.io.ObjectInputStream(fileIn)
        val snapResults2 = in.readObject.asInstanceOf[Array[AlignResult]]        

        val errorReads = List(9, 21, 55, 58, 68, 87, 95)

        errorReads.foreach(e => {
          println(e + ":  " + BruteForceAligner.verifyPairedHit(snapResults1(e), snapResults2(e), 1000))
          println(snapResults1(e))
          println(snapResults2(e))
          println
        })
      }
      case Array("testWritePairToFastq", accessKey, secretKey, bfaResultsFname) => {
        // load reads
        val sc = new SparkContext("local[8]", "bruteForceAligner")
        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)]("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + bfaResultsFname).collect

        BruteForceAligner.writePairToFastq(readPairsWithBfaResults, bfaResultsFname)
      }
      case Array("testRichManyHits") => {
        val res = RichManyHits(5, 10)
        (1 to 10).foreach(i => {
          val pos = (i*1000).toLong
          val editDist = (new scala.util.Random).nextInt(15)
          println("Pos: " + pos + ", Edit dist: " + editDist)
          res.addNewMatch(pos, editDist)
        })
        println(res)
      }
      case Array("testSimplePairedBfa", accessKey, secretKey, fastq1, fastq2, snap1, snap2) => {
        val params = new PairedBruteForceParams(15, 1, 600)
        val c = new PairedBruteForceContext("localHost", "local", accessKey, secretKey, fastq1, fastq2, snap1, snap2, "s3", "testSimpleBfa", params)
        val sc = new SparkContext(c.sparkDest, "simpleBruteForceAligner")

        simpleBfaAndSaveOutput(sc, c, params)
        checkSimpleBfaResults(sc, c, params)
      }
      case Array("testSimplePairedBfaSubset", accessKey, secretKey, fastq1, fastq2, snap1, snap2, numPairs) => {
        val params = new PairedBruteForceParams(15, 1, 600)
        val c = new PairedBruteForceContext("localHost", "local", accessKey, secretKey, fastq1, fastq2, snap1, snap2, "s3", "testSimpleBfaSubset", params)
        val sc = new SparkContext(c.sparkDest, "simpleBruteForceAligner")

        simpleBfaSubsetInputSaveOutput(numPairs.toInt, sc, c, params)
        /*
        // get subset of input
        val indices = getRandomSubset(numPairs.toInt, c.readsWithSnapResults.length)
        val readsWithSnapResultsSubset = new Array[((Read, AlignResult), (Read, AlignResult))](numPairs.toInt)

        (0 until indices.length).foreach(i => {
          readsWithSnapResultsSubset(i) = c.readsWithSnapResults(indices(i))
        })

        // run simple bfa
        // save results
        sc.parallelize(readsWithSnapResultsSubset).map(pair => {
          val ((read1, snapRes1), (read2, snapRes2)) = pair

          // set the upper bound for edit distance
          val (snapUpperBound, snapFirstReadIsRC) = params.getSnapSummary(snapRes1, snapRes2)

          // align the reads
          ((read1, read2), SimpleBruteForceAligner.pairedAlign(read1.data, read2.data, params, snapUpperBound, snapFirstReadIsRC))
        })
        .saveAsObjectFile(c.outputPath)
        */
        
        checkSimpleBfaResults(sc, c, params)
      }
      case Array("testSimplePairedBfaRealReads", accessKey, secretKey, snapSam, bwaSam, bfaBucket, maxDist) => {
        // load corresponding snap, bwa, & bfa results for real reads
        // snap
        val snapSamEntries = AlignerJudge.loadSam(snapSam).map(s => new RichSamEntry(s))
        val snapSamEntries1 = snapSamEntries.filter(e => e.ID.split("/")(1).startsWith("1"))
        val snapSamEntries2 = snapSamEntries.filter(e => e.ID.split("/")(1).startsWith("2"))
        val snapResults1 = snapSamEntries1.map(e => samEntryToAlignResult(e)).toArray
        val snapResults2 = snapSamEntries2.map(e => samEntryToAlignResult(e)).toArray

        // bwa
        val bwaSamEntries = AlignerJudge.loadSam(bwaSam).map(s => new RichSamEntry(s))
        var bwaSamEntries1 = List[RichSamEntry]()
        var bwaSamEntries2 = List[RichSamEntry]()
        (0 until bwaSamEntries.length).foreach(i => {
          if (i % 2 == 0) bwaSamEntries1 = bwaSamEntries(i) :: bwaSamEntries1
          else bwaSamEntries2 = bwaSamEntries(i) :: bwaSamEntries2
        })
        bwaSamEntries1 = bwaSamEntries1.reverse
        bwaSamEntries2 = bwaSamEntries2.reverse
        
      	(0 until bwaSamEntries.length).foreach(i => {
      	   println(i)
      	   println(bwaSamEntries(i).fields.toList.mkString(", "))
      	   samEntryToAlignResult(bwaSamEntries(i))
      	})

        val bwaResults1 = bwaSamEntries1.map(e => samEntryToAlignResult(e)).toArray
        val bwaResults2 = bwaSamEntries2.map(e => samEntryToAlignResult(e)).toArray
        
        // bfa
        val readPairsWithBfaResults = AlignerJudge.loadReadPairsAndBfaResults(None, accessKey, secretKey, None, bfaBucket)
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip

        assert(snapResults1.length == snapResults2.length)
        println("# snap pairs: " + snapResults1.length)
        assert(bwaResults1.length == bwaResults2.length)
        println("# bwa pairs:  " + bwaResults1.length)
        println("# bfa pairs:  " + bfaResults.length)

        assert(snapResults1.length == bwaResults1.length && snapResults1.length == bfaResults.length)

        // for each read, check for agreement (b/t snap & bwa, snap & bfa, bwa & bfa)
        var mismatchCount = 0

        var mismatchOnRead1 = false
        var mismatchOnRead2 = false

        (0 until snapResults1.length).foreach(i => {
          mismatchOnRead1 = false
          mismatchOnRead2 = false
          /*
          println("Read 1:")
          println(snapResults1(i))
          println(bwaResults1(i))
          println("Read 2:")
          println(snapResults2(i))
          println(bwaResults2(i))
          println(bfaResults(i))
          println
          */

          /*
          if (snapResults1(i).isInstanceOf[RichSingleHit] && bwaResults1(i).isInstanceOf[RichSingleHit] && bfaResults(i).isInstanceOf[RichPairSingleHit]) {
            assert(abs(snapResults1(i).asInstanceOf[RichSingleHit].position - bwaResults1(i).asInstanceOf[RichSingleHit].position) <= maxDist.toInt)
            assert(abs(snapResults1(i).asInstanceOf[RichSingleHit].position - bfaResults(i).asInstanceOf[RichPairSingleHit].pos1 + 1 /* 0-indexing => 1-indexing */) <= maxDist.toInt)
          }
          
          if (snapResults2(i).isInstanceOf[RichSingleHit] && bwaResults2(i).isInstanceOf[RichSingleHit] && bfaResults(i).isInstanceOf[RichPairSingleHit]) {
            assert(abs(snapResults2(i).asInstanceOf[RichSingleHit].position - bwaResults2(i).asInstanceOf[RichSingleHit].position) <= maxDist.toInt)
            assert(abs(snapResults2(i).asInstanceOf[RichSingleHit].position - bfaResults(i).asInstanceOf[RichPairSingleHit].pos2 + 1 /* 0-indexing => 1-indexing */) <= maxDist.toInt)
          }
          */

	        /*          
          // check for mismatches
          // if any one of the results is a single hit and the others aren't, that's a mismatch.
          // if they're all single hits and they don't all report the same position, that's a mismatch.
          // if they're all not found, don't do anything.
          if (needToCheckForMismatch(snapResults1(i), bwaResults1(i), bfaResults(i))) {
            if (mismatch(snapResults1(i).asInstanceOf[RichSingleHit], bwaResults1(i).asInstanceOf[RichSingleHit], getFirstHit(bfaResults(i).asInstanceOf[RichPairSingleHit]), maxDist.toInt)) {
              mismatchOnRead1 = true
              mismatchCount += 1
              println("Mismatch at pair " + i + ", read 1")
              println("SNAP: " + snapResults1(i))
              println("BWA:  " + bwaResults1(i))
              println("BFA:  " + bfaResults(i))
              println
            }
          }
          
          if (needToCheckForMismatch(snapResults2(i), bwaResults2(i), bfaResults(i))) {
            if (mismatch(snapResults2(i).asInstanceOf[RichSingleHit], bwaResults2(i).asInstanceOf[RichSingleHit], getSecondHit(bfaResults(i).asInstanceOf[RichPairSingleHit]), maxDist.toInt)) {
              mismatchOnRead2 = true
              mismatchCount += 1
              println("Mismatch at pair " + i + ", read 2")
              println("SNAP: " + snapResults2(i))
              println("BWA:  " + bwaResults2(i))
              println("BFA:  " + bfaResults(i))
              println
            }
          }
	        */

	        /*          
          if (mismatchOnRead1 || mismatchOnRead2) {
            // verify alignments at all 3 sets of positions
            // get snap score
            if (snapResults1(i).isInstanceOf[RichSingleHit] && snapResults2(i).isInstanceOf[RichSingleHit]) {
              val snapScore = AlignmentVerifier.getBestScore(reads1(i).data, reads2(i).data, snapResults1(i).asInstanceOf[RichSingleHit].position, snapResults2(i).asInstanceOf[RichSingleHit].position, maxDist.toInt)
              println("Snap score: " + snapScore)
            }
            
            // get bwa score
            if (bwaResults1(i).isInstanceOf[RichSingleHit] && bwaResults2(i).isInstanceOf[RichSingleHit]) {
              val bwaScore = AlignmentVerifier.getBestScore(reads1(i).data, reads2(i).data, bwaResults1(i).asInstanceOf[RichSingleHit].position, bwaResults2(i).asInstanceOf[RichSingleHit].position, maxDist.toInt)
              println("Bwa score: " + bwaScore)
            }
            
            // get bfa score
            if (bfaResults(i).isInstanceOf[RichPairSingleHit]) {
              val res = bfaResults(i).asInstanceOf[RichPairSingleHit]
              val bfaScore = AlignmentVerifier.getBestScore(reads1(i).data, reads2(i).data, res.pos1, res.pos2, maxDist.toInt)
              println("Bfa score: " + bfaScore)
            }
          }
	        */

          //if (i == 8) {
      	  println("Snap:")
      	  println(snapResults1(i))
      	  println(snapResults2(i))
      	  if (snapResults1(i).isInstanceOf[RichSingleHit] && snapResults2(i).isInstanceOf[RichSingleHit]) {
            val snapScore = AlignmentVerifier.getBestScore(reads1(i).data, reads2(i).data, snapResults1(i).asInstanceOf[RichSingleHit].position-1, snapResults2(i).asInstanceOf[RichSingleHit].position-1, maxDist.toInt)
            println("Snap score" + i + ": " + snapScore)
          }
	        println

      	  println("BWA:")
      	  println(bwaResults1(i))
      	  println(bwaResults2(i))
      	  if (bwaResults1(i).isInstanceOf[RichSingleHit] && bwaResults2(i).isInstanceOf[RichSingleHit]) {
              val bwaScore = AlignmentVerifier.getBestScore(reads1(i).data, reads2(i).data, bwaResults1(i).asInstanceOf[RichSingleHit].position-1, bwaResults2(i).asInstanceOf[RichSingleHit].position-1, maxDist.toInt)
              println("Bwa score: " + bwaScore)
          }
	        println

      	  println("BFA:")
      	  println(bfaResults(i))
        })
      }
      case Array("testSimpleBfaOnOldBfaMultiHits", accessKey, secretKey, bfaBucket) => {
        // load bfa results
        val readPairsWithBfaResults = AlignerJudge.loadReadPairsAndBfaResults(None, accessKey, secretKey, None, bfaBucket)
        
        // pull out any multi hit results where both best score & second best score are 0
        val readPairsWithResultsOfInterest = 
        readPairsWithBfaResults.flatMap(t => {
          t._2 match {
            case RichPairMultiHit(_, _, bestScore, _, _, secondBestScore, _) => {
              if (bestScore == 0 && secondBestScore == 0)
                List(t)
              else
                Nil
            }
            case _ => Nil
          }
        })
        
        // take first 4 from list
        assert(readPairsWithResultsOfInterest.length >= 4)
        val readPairsWithResultsToTest = readPairsWithResultsOfInterest.take(4)
        
        // call brute force
        val sc = new SparkContext("local[4]", "simpleBruteForceAligner")
        val params = new PairedBruteForceParams(20, 1, 1000)
        
        sc.parallelize(readPairsWithResultsToTest).map(pair => {
          val (readPair, bfaResult) = pair
          val (read1, read2) = readPair

          // set the upper bound for edit distance
          val (snapUpperBound, snapFirstReadIsRC) = (5, bfaResult.asInstanceOf[RichPairMultiHit].isRC)

          // align the reads
          ((read1, read2), SimpleBruteForceAligner.pairedAlign(read1.data, read2.data, params, snapUpperBound, snapFirstReadIsRC))
        })
        .collect
        .foreach(println)
      }

      // Debugging
      case Array("examinePairedBFAResults", master, accessKey, secretKey, bfaResultsFname) => {
        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")
        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)]("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + bfaResultsFname).collect
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip

        var resNum = 0
        bfaResults.foreach(res => {
          println("res# :" + resNum)
          res match {
            case RichPairSingleHit(_, _, _, _) => println(res)
            case _ =>
          }

          println
          resNum += 1
        })

        System.exit(0)
      }
      case Array("examinePairedBFAErrors", master, bfaResultsFname, maxDist, confDiff) => {
        // load reads
        val sc = new SparkContext("1@" + master + ":5050", "bruteForceAligner")
        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)]("hdfs://" + master + ":9000/" + bfaResultsFname).collect
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip

        // test them against the locations that wgsim reported
        val readNums = List(0, 55, 87, 91)
        val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r

        readNums.foreach(i => {
          println("Read#: " + i)
          val read1 = reads1(i)
          val read2 = reads2(i)

          println("idStr: " + read1.idStr)

          read1.idStr match {
            case ID_REGEX(piece, start, end) => {
              val lowEnd = min(start.toLong, end.toLong)
              val highEnd = max(start.toLong, end.toLong)

              println("lowEnd: " + lowEnd)
              println("highEnd: " + highEnd)

              val res1 = BruteForceAligner.alignInRegion(read1.data, piece, lowEnd, highEnd, maxDist.toInt, maxDist.toInt, confDiff.toInt, false)
              val res2 = BruteForceAligner.alignInRegion(read2.data, piece, lowEnd, highEnd, maxDist.toInt, maxDist.toInt, confDiff.toInt, false)

              println("wgsim says: ")
              println(res1)
              println(res2)
              println
            }
            case _ =>
          }

          println("bfa says: ")
          println(bfaResults(i))
          println
        })
      }
      case Array("debugBFAErrors", accessKey, secretKey, bfaResultsFname) => {
        // load reads
        val outputDest = "s3n://" + accessKey + ":" + secretKey + "@hashingGenome/"
        val sc = new SparkContext("local[8]", "bruteForceAligner")
        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)](outputDest + bfaResultsFname).collect
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip

        // look at bfa errors
        val readNums = List(1167, 534, 1972, 969, 816, 1579, 1853, 96, 1783, 454, 1833, 1906, 1876, 681, 1860, 1847, 1179, 1728, 1383)

        val maxDist = 10
        val confDiff = 1
        val separationDist = 1000
        val lv = new LandauVishkin(2 * maxDist + confDiff - 1)
        val outputPrefix = "bfaErrors"

        sc.parallelize(readNums).map(n => {
          val read1 = reads1(n)
          val read2 = reads2(n)

          // set the upper bound for edit distance
          var snapUpperBound = 12

          // align the reads
          ((read1, read2), BruteForceAligner.pairedAlign(read1.data, read2.data, snapUpperBound, maxDist, confDiff, separationDist, false))
        })
        .saveAsObjectFile(outputDest + outputPrefix)

        // look at the results
        val errorPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)](outputDest + outputPrefix).collect
        errorPairsWithBfaResults.foreach(pair => {
          val (readPair, result) = pair
          val (read1, read2) = readPair

          println(read1.dataStr)
          println(read2.dataStr)
          println(result)
          println
        })
      }
      case Array("analyzeBFAOutput", accessKey, secretKey, bfaResultsFname) => {
        val outputDest = "s3n://" + accessKey + ":" + secretKey + "@hashingGenome/"
        val sc = new SparkContext("local[8]", "bruteForceAligner")

        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)](outputDest + bfaResultsFname).collect
        val (readPairs, results) = readPairsWithBfaResults.unzip

        var single = 0
        var multi = 0
        var none = 0

        results.foreach(r => {
          r match {
            case RichPairSingleHit(_, _, _, _) => single += 1
            case RichPairMultiHit(_, _, bestScore, _, _, secondBestScore, _) => {
              if (bestScore < secondBestScore)
                single += 1
              else
                multi += 1
            }
            case NotFound => none += 1
            case _ =>
          }
        })

        println("# reads:     " + results.length)
        println("# single:    " + single + " (" + round(single.toDouble / results.length * 100.0) + ")")
        println("# multi:     " + multi + " (" + round(multi.toDouble / results.length * 100.0) + ")")
        println("# not found: " + none + " (" + round(none.toDouble / results.length * 100.0) + ")")
      }
      case Array("analyzeSnapResults", snapFname) => {
        val fileIn = new java.io.FileInputStream(snapFname)
        val in = new java.io.ObjectInputStream(fileIn)
        val snapResults = in.readObject.asInstanceOf[Array[AlignResult]]

        var max = 0

        snapResults.foreach(r => {
          r match {
            case RichSingleHit(_, _, editDist) =>
              if (editDist > max)
                max = editDist
            case _ =>
          }
        })

        println("Max: " + max)
      }
      case Array("lookAtSnapSimResults", snapFname1, snapFname2) => {
        var fileIn = new java.io.FileInputStream(snapFname1)
        var in = new java.io.ObjectInputStream(fileIn)
        val snapResults1 = in.readObject.asInstanceOf[Array[AlignResult]]

        fileIn = new java.io.FileInputStream(snapFname2)
        in = new java.io.ObjectInputStream(fileIn)
        val snapResults2 = in.readObject.asInstanceOf[Array[AlignResult]]

        val pairResults = snapResults1.zip(snapResults2)

        val separationDist = 1000

        var resNum = 0

        pairResults.foreach(r => {
          val (res1, res2) = r

          if (resNum == 95) {
            println(res1)
            println(res2)
          }

          resNum += 1
        })
      }
      case Array("examineErrorsBWAandSNAP", accessKey, secretKey, bwaFname, snapFname, bfaFname) => {
        // load errors (sam files)
        val bwaSource = scala.io.Source.fromFile(bwaFname)
        val bwaLines = bwaSource.getLines

        val snapSource = scala.io.Source.fromFile(snapFname)
        val snapLines = snapSource.getLines

        // get pairs
        var bwaSet = scala.collection.mutable.Set[(String, (String, Long))]()
        bwaLines.foreach(l => {
          val entries = l.split("\t")
          val id = entries(0).split("/")(0)
          val piece = entries(2)
          val pos = entries(3).toLong
          bwaSet.add((id, (piece, pos)))
        })

        val bwaMap = bwaSet.groupBy(t => {
          val entries = t._1.split("/")
          entries(0)
        })

        var snapSet = Set[(String, (String, Long))]()
        snapLines.foreach(l => {
          if (!l.startsWith("@")) {
            val entries = l.split("\t")
            val id = entries(0).split("/")(0)
            val piece = entries(2)
            val pos = entries(3).toLong
            snapSet.add((id, (piece, pos)))
          }
        })

        val snapMap = snapSet.groupBy(t => {
          val entries = t._1.split("/")
          entries(0)
        })


        val allErrors = bwaMap.keySet ++ snapMap.keySet
        val errorsInCommon = bwaMap.keySet.intersect(snapMap.keySet)

        println("# error pairs: " + allErrors.size)
        println("# bwa error pairs: " + bwaMap.keySet.size)
        println("# snap error pairs: " + snapMap.keySet.size)
        println("# error pairs in common: " + errorsInCommon.size)
        println

        var sameAnsCt = 0
        errorsInCommon.foreach(e => {
          val bwaRes = bwaMap.get(e)
          val snapRes = snapMap.get(e)
          if (bwaRes == snapRes) 
            sameAnsCt += 1

          println("bwa: " + bwaMap.get(e))
          println("snap: " + snapMap.get(e))
          /*
          else {
            println("bwa: " + bwaMap.get(e))
            println("snap: " + snapMap.get(e))
            println
          }
          */
        })

        println("errors in common:")
        println("agreement: " + sameAnsCt + " (" + math.round(sameAnsCt.toDouble/errorsInCommon.size * 100) + "%)")

        bwaSource.close
        snapSource.close

        // where bwa & snap agree, investigate
        // load bfa results on errors
        val sc = new SparkContext("local[8]", "bruteForceAligner")
        val readPairsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)]("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + bfaFname).collect
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip

        var readNums = List[Int]()
        var errorNum = 0

        errorsInCommon.foreach(e => {
          val bwaRes = bwaMap.get(e)
          val snapRes = snapMap.get(e)
          //if (bwaRes == snapRes) {
            // first, what does bfa say (pos, edit distance)?
            // next, check the read against the position that bwa & snap say
            // what is the edit distance?  is it < the edit distance bfa found?

            // pull out corresponding read from bfa
            val read = reads1.filter(r => {
              val rID = r.idStr
              val prefix = rID.split("/")(0)

              // get out id part of e
              val i = e.indexOf(":")
              val ePrefix = e.substring(i + 1, e.length)

              prefix == ePrefix
            }).head

            // pull out corresponding bfa result
            val readNum = reads1.indexOf(read)
            readNums = readNum :: readNums
            val bfaResult = bfaResults(readNum)
            errorNum += 1

            //if (readNum == 1579) {
              // see what they say
              println("readNum: " + readNum)
              println("bfa says:  " + bfaResult)
              println("bwa says:  " + bwaMap.get(e))
              println("snap says: " + snapMap.get(e))

              // check score reported by snap & bwa
              // get pos reported by snap & bwa
              (1 to 2).foreach(i => {
                val posList = 
                if (i == 1) bwaRes.get.map(t => t._2).toList
                else snapRes.get.map(t => t._2).toList

                  if (posList.length == 2) {
                    val (pieceName1, pos1) = posList(0)
                    val (pieceName2, pos2) = posList(1)

                    /*
                    if (piece1 != piece2) {
                      println("Pieces not equal")
                      System.exit(0)
                    }
                    */

                    val piece1 = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(pieceName1))
                    val piece2 = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(pieceName2))

                    val absPos1 = piece1.startIndex + pos1 - 1
                    val absPos2 = piece2.startIndex + pos2 - 1

                    // get reads
                    val read1 = reads1(readNum)
                    val read2 = reads2(readNum)
                    val readLen = read1.data.length
                    val rc1 = new Array[Byte](readLen)
                    DNA.getReverseComplement(read1.data, rc1)
                    val rc2 = new Array[Byte](readLen)
                    DNA.getReverseComplement(read2.data, rc2)

                    // check edit dist
                    val maxDist = 30
                    val lv = new LandauVishkin(maxDist)

                    val ref1 = new Array[Byte](readLen)
                    GenomeLoader.genome.getSubstring(absPos1, absPos1 + readLen, ref1)

                    val ref2 = new Array[Byte](readLen)
                    GenomeLoader.genome.getSubstring(absPos2, absPos2 + readLen, ref2)

                    // pos1
                    val score1f_1 = lv.distance(ref1, readLen, read1.data, readLen, maxDist)
                    val score2f_1 = lv.distance(ref1, readLen, read2.data, readLen, maxDist)

                    val score1r_1 = lv.distance(ref1, readLen, rc1, readLen, maxDist)
                    val score2r_1 = lv.distance(ref1, readLen, rc2, readLen, maxDist)

                    // pos
                    val score1f_2 = lv.distance(ref2, readLen, read1.data, readLen, maxDist)
                    val score2f_2 = lv.distance(ref2, readLen, read2.data, readLen, maxDist)

                    val score1r_2 = lv.distance(ref2, readLen, rc1, readLen, maxDist)
                    val score2r_2 = lv.distance(ref2, readLen, rc2, readLen, maxDist)

                    val individScores = List(score1f_1, score2r_2, score2f_1, score1r_2, score1f_2, score2r_1, score2f_2, score1r_1)
                    val scores = List(score1f_1 + score2r_2, score2f_1 + score1r_2, score1f_2 + score2r_1, score2f_2 + score1r_1)
                    println("manual:")
                    println(individScores)
                    println(scores)
                    println("pos1: " + absPos1 + " or " + pieceName1 + "(" + pos1 + ")")
                    println("pos2: " + absPos2 + " or " + pieceName2 + "(" + pos2 + ")")

                    /*
                    // print out reads & rc
                    println("read1 (f): " + new String(read1.data))
                    println("read1 (r): " + new String(rc1))
                    println("read2 (f): " + new String(read2.data))
                    println("read2 (r): " + new String(rc2))

                    // also print out reference
                    println("ref at pos1: " + new String(ref1))
                    println("ref at pos2: " + new String(ref2))
                    */

                    // also get edit dist at pos from id
                    // HACK
                    /*
                    val p = "chr18"
                    val pos1A = 108797
                    val pos2A = 108968

                    val pieceA = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(p))

                    val absPos1A = pieceA.startIndex + pos1A
                    val absPos2A = pieceA.startIndex + pos2A

                    val ref1A = new Array[Byte](readLen)
                    GenomeLoader.genome.getSubstring(absPos1A, absPos1A + readLen, ref1A)

                    val ref2A = new Array[Byte](readLen)
                    GenomeLoader.genome.getSubstring(absPos2A, absPos2A + readLen, ref2A)

                    // pos1
                    val score1f_1A = lv.distance(ref1A, readLen, read1.data, readLen, 20)
                    val score2f_1A = lv.distance(ref1A, readLen, read2.data, readLen, 20)

                    val score1r_1A = lv.distance(ref1A, readLen, rc1, readLen, 20)
                    val score2r_1A = lv.distance(ref1A, readLen, rc2, readLen, 20)

                    // pos
                    val score1f_2A = lv.distance(ref2A, readLen, read1.data, readLen, 20)
                    val score2f_2A = lv.distance(ref2A, readLen, read2.data, readLen, 20)

                    val score1r_2A = lv.distance(ref2A, readLen, rc1, readLen, 20)
                    val score2r_2A = lv.distance(ref2A, readLen, rc2, readLen, 20)

                    val individScoresA = List(score1f_1A, score2r_2A, score2f_1A, score1r_2A, score1f_2A, score2r_1A, score2f_2A, score1r_1A)
                    val scoresA = List(score1f_1A + score2r_2A, score2f_1A + score1r_2A, score1f_2A + score2r_1A, score2f_2A + score1r_1A)
                    println("manual:")
                    println(individScoresA)
                    println(scoresA)
                    println("pos1: " + absPos1A + " or " + pieceA + "(" + pos1A + ")")
                    println("pos2: " + absPos2A + " or " + pieceA + "(" + pos2A + ")")
                    */
                  }
                //}

                println

              })
          //}
        })        

        //println(readNums)
      }
      case Array("examineErrorsFromSam", fname) => {
        val source = scala.io.Source.fromFile(fname)
        val lines = source.getLines

        val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r

        lines.foreach(l => {
          if (!l.startsWith("@")) {
            println(l)
            val entries = l.split("\t")

            // get id => bfa ans
            val idEntries = entries(0).split("_")
            var bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(idEntries(0)))
            var bfaPos1 = idEntries(1).toLong
            var bfaPos2 = idEntries(2).toLong - 100

            // get error ans
            val errorPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(entries(2)))
            val errorPos = entries(3).toLong - 1

            // get read & cf to bfa ans, error ans
            val read = entries(9).getBytes

            // cf to bfa ans
            val readLen = read.length
            val rc = new Array[Byte](readLen)
            DNA.getReverseComplement(read, rc)

            val posList = List(bfaPiece.startIndex + bfaPos1, bfaPiece.startIndex + bfaPos2, errorPiece.startIndex + errorPos)
            val maxDist = 30
            val lv = new LandauVishkin(maxDist)

            // compare read & its rc to each pos
            posList.foreach(pos => {
              val ref = new Array[Byte](readLen)
              GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)

              // pos1
              val scoreF = lv.distance(ref, readLen, read, readLen, maxDist)
              val scoreR = lv.distance(ref, readLen, rc, readLen, maxDist)

              println("Pos: " + pos)
              println("ScoreF: " + scoreF)
              println("ScoreR: " + scoreR)
              println
            })
          }
        })
        
        source.close
      }
      case Array("scoreSam", fname) => {
        val source = scala.io.Source.fromFile(fname)
        val lines = source.getLines

        val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r

        var nCorrect = 0
        var nTotal = 0
        var nConfident = 0
        var nClipped = 0

        lines.foreach(l => {
          if (!l.startsWith("@")) {
            val entries = l.split("\t")

            // get id => bfa ans
            val idEntries = entries(0).split(":")(0).split("_")
            
            var bfaPiece: GenomePiece = null
            var bfaPos1 = 0L
            var bfaPos2 = 0L
            
            if (idEntries.length == 6) {
              bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(idEntries(0) + "_" + idEntries(1) + "_" + idEntries(2)))
              bfaPos1 = idEntries(3).toLong
              bfaPos2 = idEntries(4).toLong
            } else if (idEntries.length == 5) {
              bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(idEntries(0) + "_" + idEntries(1)))
              bfaPos1 = idEntries(2).toLong
              bfaPos2 = idEntries(3).toLong
            } else if (idEntries.length == 4) {
              bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(idEntries(0)))
              bfaPos1 = idEntries(1).toLong
              bfaPos2 = idEntries(2).toLong
            }

            // get sam ans
            if (entries(2) != "*") {
              val samPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(entries(2)))
              val samPos = entries(3).toLong - 1

              nConfident += 1
              val gap = 25
              if (bfaPiece == samPiece && (abs(bfaPos1 - samPos) < gap || abs(bfaPos2 - samPos) < gap)) {
                nCorrect += 1
              } else {
                // check for clipping
                val qualityStr = entries(entries.length - 1)
                // count how many #'s the quality string starts with
                val i = qualityStr.indexWhere(c => c != '#')
                if (i > 0) {
                  val clippedPos = samPos - i
                  if (bfaPiece == samPiece && (abs(bfaPos1 - clippedPos) < gap || abs(bfaPos2 - clippedPos) < gap))
                    nClipped += 1
                  else
                    println(l)
                } else println(l)
              } 
            }

            nTotal += 1
          }
        })
        
        println("% correct: " + (nCorrect.toDouble / nConfident) * 100)
        println("% aligned: " + (nConfident.toDouble / nTotal) * 100)
        println("% clipped: " + (nClipped.toDouble / nTotal) * 100)
        println("% correct (adj): " + (nCorrect.toDouble + nClipped.toDouble) / nConfident.toDouble * 100)
        println("# correct: " + nCorrect + ", # clipped: " + nClipped + ", # confident: " + nConfident)
      }
      case Array("debugError") => {
        /*
        var bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId("chr17"))
        var bfaPos1 = 22254553.toLong
        var bfaPos2 = 22254752.toLong

        // get error ans
        val errorPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId("chr17"))
        val errorPos1 = 22252175.toLong
        val errorPos2 = 22252373.toLong

        // get read & cf to bfa ans, error ans
        val read1 = "CTCAGAAACTTCTCTGTGATGATTGCATTCAACTCACAGAGTTGAACCCTCCTATGGATAGAGCAGTGTTGAAACTCTCTTTTTGTGGAATCTGCAAGTG".getBytes
        val rc2 = "TTCAACTCCCAGAGTTTCACGTTGCTTTTCATAGAGTAGTTCTGAAACATGCTTTTCGTAGTGTCTGCAAGTGGACATTTGGAGCGCTTTCAGGCCTGTG".getBytes
        */
        
        /*
        var bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId("chr9"))
        val bfaPos1 = 67745396.toLong
        val bfaPos2 = 67745598.toLong
        
        val errorPiece = bfaPiece
        val errorPos1 = 45039265.toLong
        val errorPos2 = 45039064.toLong

        val read1 = "CCCAGGAGGTGGAGGTTGCATTGAGCTGAAATCATGCCATTGCACTCCAGCCTTGGCAACAAGAGCAAAACTCCATCCCAAAAACAAATAAACAACAACA".getBytes
        val rc2 = "GCCGGGTGCAGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGTCGGGGCGGGCGGGGCGCAAGGGCAAGAGTTTTGGACCAGCCTGGCCAACATGGC".getBytes
        */
        
        /*
        var bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId("chr1"))
        val bfaPos1 = 121485085.toLong
        val bfaPos2 = 121485267.toLong
        
        val errorPiece = bfaPiece
        val errorPos1 = 121484067.toLong
        val errorPos2 = 121484248.toLong
        
        val read1 = "GAAGAATTCTCAGAATCTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATT".getBytes
        val rc2 = "CAGAAACTTCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAACTCTGCAAGTGGA".getBytes
        */

        /*
        var bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(""))
        val bfaPos1 = .toLong
        val bfaPos2 = .toLong
        
        val errorPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(""))
        val errorPos1 = .toLong
        val errorPos2 = .toLong
        
        val read1 = "".getBytes
        val rc2 = "".getBytes
        */
        
        var bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId("chr1"))
        val bfaPos1 = 74236056.toLong
        val bfaPos2 = (74236268 - 1).toLong
        
        val errorPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId("chr1"))
        val errorPos1 = (74236057 - 1).toLong
        val errorPos2 = (74236293 - 1).toLong
        
        val read1 = "GGAGACCCCACCCCAGGCTGCAAGGGGGCACTGCTGGGGCTGCATGCTCCACGGAGCTGGCAAGATCTGGGAGCATGTGGCAGCCCCACCCTCCCAGGTG".getBytes
        val rc2 = "GTGCTCACATTCCAGTCGCCTGCTGCCTCAACACCCTTTGGAATTTAGGTGACAGTGAGTGTGGGTGGACAAGCCAAGGGTGGGCTGAGGGCATCTCGGC".getBytes

        val posList = List(bfaPiece.startIndex + bfaPos1, bfaPiece.startIndex + bfaPos2, errorPiece.startIndex + errorPos1, errorPiece.startIndex + errorPos2)
        val maxDist = 30
        val lv = new LandauVishkin(maxDist)
        
        // compare read & its rc to each pos
        val readLen = 100
        val ref = new Array[Byte](readLen)

        posList.foreach(pos => {
          GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)

          // pos1
          val scoreF = lv.distance(ref, readLen, read1, readLen, maxDist)
          val scoreR = lv.distance(ref, readLen, rc2, readLen, maxDist)

          println("Pos: " + pos)
          println("ScoreF: " + scoreF)
          println("ScoreR: " + scoreR)
          println
        })
      }
      case Array("verifyErrorsFromSam", fname) => {
        val source = scala.io.Source.fromFile(fname)
        val lines = source.getLines

        val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r

        var nCorrect = 0
        var nTotal = 0
        var nConfident = 0
        var nClipped = 0
        var nError = 0

        var error = false

        lines.foreach(l => {
          if (!l.startsWith("@")) {
            error = false
            
            val entries = l.split("\t")

            // get id => bfa ans
            val idEntries = entries(0).split(":")(0).split("_")
            
            var bfaPiece: GenomePiece = null
            var bfaPos1 = 0L
            var bfaPos2 = 0L
            
            if (idEntries.length == 6) {
              bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(idEntries(0) + "_" + idEntries(1) + "_" + idEntries(2)))
              bfaPos1 = idEntries(3).toLong
              bfaPos2 = idEntries(4).toLong - 1
            } else if (idEntries.length == 5) {
              bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(idEntries(0) + "_" + idEntries(1)))
              bfaPos1 = idEntries(2).toLong
              bfaPos2 = idEntries(3).toLong - 1
            } else if (idEntries.length == 4) {
              bfaPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(idEntries(0)))
              bfaPos1 = idEntries(1).toLong
              bfaPos2 = idEntries(2).toLong - 1
            }

            // get sam ans
            var samPiece: GenomePiece = null
            var samPos = 0L
            
            if (entries(2) != "*") {
              samPiece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(entries(2)))
              samPos = entries(3).toLong - 1

              nConfident += 1
              val gap = 25  // may want to change this
              if (bfaPiece == samPiece && (abs(bfaPos1 - samPos) < gap || abs(bfaPos2 - samPos) < gap)) {
                nCorrect += 1
              } else {
                // check for clipping
                val qualityStr = entries(entries.length - 1)
                // count how many #'s the quality string starts with
                val i = qualityStr.indexWhere(c => c != '#')
                if (i > 0) {
                  val clippedPos = samPos - i
                  if (bfaPiece == samPiece && (abs(bfaPos1 - clippedPos) < gap || abs(bfaPos2 - clippedPos) < gap))
                    nClipped += 1
                  else
                    error = true
                } else error = true
              } 
            }

            nTotal += 1
            
            // verify errors
            if (error) {
              // only count the error if bfa's results aren't too far/too close
              if(
                (bfaPos2 - bfaPos1) < 500                                 // filter out cases where bfa's positions are too far apart
                && (bfaPos2 - bfaPos1) > 100                              // filter out cases where bfa's positions are too far close
              ) 
                nError += 1
                
              
              val pairPieceName = entries(6)
              val pairPos = entries(7).toLong - 1
              if (
                (pairPieceName == samPiece.name || pairPieceName == "=")  // filter out cases where aligner being tested put the pair on different chromosomes
                && pairPos != samPos                                      // filter out cases where aligner mapped both pos to same location
                && abs(pairPos - samPos) < 1000                           // filter out cases where aligner mapped pair to too-distant locations
                && (bfaPos2 - bfaPos1) < 500                              // filter out cases where bfa's positions are too far apart
                && (bfaPos2 - bfaPos1) > 100                              // filter out cases where bfa's positions are too far close
              ) { 
                println(l)
                // want to verify b/c this read & its pair were placed on same chromosome
                val posList = List(bfaPiece.startIndex + bfaPos1, bfaPiece.startIndex + bfaPos2, samPiece.startIndex + samPos)
                val posNames = List("brute force (1)", "brute force (2)", "sam")
                val maxDist = 30
                val lv = new LandauVishkin(maxDist)

                // compare read & its rc to each pos
                val read = entries(9).getBytes
                val readLen = read.length
                var rc = new Array[Byte](readLen)
                DNA.getReverseComplement(read, rc)
                
                val ref = new Array[Byte](readLen)

                (0 until posList.length).foreach(i => {
                  val pos = posList(i)
                  GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)

                  // pos1
                  val scoreF = lv.distance(ref, readLen, read, readLen, maxDist)
                  val scoreR = lv.distance(ref, readLen, rc, readLen, maxDist)

                  println(posNames(i))

                  val (myPiece, myOffset) = GenomeLoader.genome.getLocation(pos)
                  
                  println("Pos: " + myPiece + "(" + myOffset + ")")
                  println("ScoreF: " + scoreF)
                  println("ScoreR: " + scoreR)
                  println
                })
              }
              
              
            }
          }
        })
        
        println("% error: " + (nError.toDouble / nConfident) * 100)
      }
      case Array("lookIntoBwaIssue") => {
        // bwa result
        //case class RichPairSingleHit(pos1: Long, pos2: Long, score: Int, isRC: Boolean) extends AlignResult {
        val resBwa = RichPairSingleHit(GenomeLoader.genome.getAbsPos("chr5", 178928349), GenomeLoader.genome.getAbsPos("chr5", 178928171), 1, true)
        
        // bfa result
        val resBfa1 = RichPairSingleHit(GenomeLoader.genome.getAbsPos("chr5", 178876229), GenomeLoader.genome.getAbsPos("chr5", 178876407), 1, false)
        val resBfa2 = RichPairSingleHit(GenomeLoader.genome.getAbsPos("chr5", 178928347), GenomeLoader.genome.getAbsPos("chr5", 178928170), 3, true /*actually I don't know*/)
        
        // verify each
        val read1 = "AAAGCCAGTCACTGAGACAAGGAGCATTGCCAGAGAAGAAGCTTTAATCCCGTGCTGCAGCCCAGGAGATGGGAGATCAGTCTCAGCCCATCCCCTGACC".getBytes
        val read2 = "AGGCCTCCCAGATAGAATAAAAAAAAAAACTCCCAAGACCACCACACCCTCCAGCCTGGGAATTGCCTAACCCACCACCTGCTTCCTGGTGACCGACTCT".getBytes
        
        //def getBestScore(read1: Array[Byte], read2: Array[Byte], pos1: Long, pos2: Long, maxDist: Int): Int = {
        // bwa
        AlignmentVerifier.getBestScore(read1, read2, resBwa.pos1 - 1, resBwa.pos2 - 1, 20)  // 1-indexed
        println
        AlignmentVerifier.getBestScore(read1, read2, resBfa1.pos1, resBfa1.pos2, 20)  // 0-indexed
        println
        AlignmentVerifier.getBestScore(read1, read2, resBfa2.pos1, resBfa2.pos2, 20)  // 0-indexed
        println
        
        // redo alignment & see what happens
        val params = new PairedBruteForceParams (20, 3, 1000)
        SimpleBruteForceAligner.pairedAlign(read1, read2, params, 4, true, true)
      }

      // Misc (WIP)
      case Array("buildGraph", master, readLen, seedLen, seedsToTry, maxDist, maxHits) => {
        // index genome
        println("Indexing genome...")
        val idxStart = System.currentTimeMillis
        val builder = new IntHashIndexBuilder(seedLen.toInt, GenomeLoader.genome.totalSize)
        GenomeLoader.genome.addToIndex(builder)
        val index = builder.build()
        printf("Indexing took %.3fs%n", (System.currentTimeMillis - idxStart) / 1000.0)

        // align "reads" from genome (each 100-base substring in genome)
        // foreach read
        // want to find all locations in the genome that it is no more than <maxDist> edits (only substitutions, not indels) away from
        // modified snap
        val sc = new SparkContext("1@" + master + ":5050", "buildGraph")
        val genomeSize = GenomeLoader.genome.totalSize
        sc.parallelize(0 until (genomeSize/2).toInt).flatMap(pos => {
          // try seedsToTry seeds from the read to get candidate locations
          // if any (or <some min #>) of the seeds violates maxHits, fall back to bfa (check each location)
          // -- gets around problem that index only includes some of the hits for popular seeds
          // else try candidates only
          // keep track of how many you had to do bfa for (should be a small fraction)

          val pos1 = 2 * pos.toLong
          val l1 = BruteForceAligner.multiAlign(index, pos1, readLen.toInt, seedLen.toInt, seedsToTry.toInt, maxDist.toInt, maxHits.toInt)

          val pos2 = 2 * pos.toLong + 1
          val l2 = BruteForceAligner.multiAlign(index, pos2, readLen.toInt, seedLen.toInt, seedsToTry.toInt, maxDist.toInt, maxHits.toInt)

          l1 ++ l2
        })
        .saveAsObjectFile("hdfs://" + master + ":9000/edges")
      }

      case _ => {
        println("Incorrect parameters.")
      }
    }
  }
  
  def s3Path(accessKey: String, secretKey: String, baseBucketName: Option[String], subBucketName: String): String = {
    val baseBucket = 
    baseBucketName match {
      case Some(name) => name
      case None => "hashingGenome"
    }
    
    "s3n://" + accessKey + ":" + secretKey + "@" + baseBucket + "/" + subBucketName
  }

  def simpleBfaAndSaveOutput(sc: SparkContext, c: PairedBruteForceContext, params: PairedBruteForceParams) = {
    sc.parallelize(c.readsWithSnapResults).map(pair => {
      val ((read1, snapRes1), (read2, snapRes2)) = pair

      // set the upper bound for edit distance
      val (snapUpperBound, snapFirstReadIsRC) = params.getSnapSummary(snapRes1, snapRes2)

      // align the reads
      ((read1, read2), SimpleBruteForceAligner.pairedAlign(read1.data, read2.data, params, snapUpperBound, snapFirstReadIsRC))
    })
    .saveAsObjectFile(c.outputPath)
  }

  def getRandomSubset(subsetSize: Int, totalSize: Int): List[Int] = {
    assert(subsetSize <= totalSize)
    
    val which = scala.collection.mutable.Set[Int]()
    while (which.size < subsetSize) {
      which += scala.util.Random.nextInt(totalSize)
    }
    
    which.toList.sortWith(_ < _)
  }

  def simpleBfaSubsetInputSaveOutput(numPairs: Int, sc: SparkContext, c: PairedBruteForceContext, params: PairedBruteForceParams) = {
    // get subset of input
    val indices = getRandomSubset(numPairs, c.readsWithSnapResults.length)
    val readsWithSnapResultsSubset = new Array[((Read, AlignResult), (Read, AlignResult))](numPairs)
    
    (0 until indices.length).foreach(i => {
      readsWithSnapResultsSubset(i) = c.readsWithSnapResults(indices(i))
    })
    
    // run simple bfa
    // save results
    sc.parallelize(readsWithSnapResultsSubset).map(pair => {
      val ((read1, snapRes1), (read2, snapRes2)) = pair

      // set the upper bound for edit distance
      val (snapUpperBound, snapFirstReadIsRC) = params.getSnapSummary(snapRes1, snapRes2)

      // align the reads
      ((read1, read2), SimpleBruteForceAligner.pairedAlign(read1.data, read2.data, params, snapUpperBound, snapFirstReadIsRC))
    })
    .saveAsObjectFile(c.outputPath)
  }

  def checkSimpleBfaResults(sc: SparkContext, c: PairedBruteForceContext, params: PairedBruteForceParams) = {
    // check results
    val readsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)](c.outputPath).collect
    readsWithBfaResults.foreach(t => {
      val ((read1, read2), res) = t
      println(read1.idStr)
      println(res)
      res match {
        case RichPairSingleHit(pos1, pos2, score, isRC) => {
          assert(Wgsim.isCorrect(read1.idStr, pos1, params.maxDist))
          assert(Wgsim.isCorrect(read2.idStr, pos2, params.maxDist))
          println("Align at brute force pair: " + AlignmentVerifier.getBestScore(read1.data, read2.data, pos1, pos2, params.maxDist))
          println("Align at wgsim pair: " + AlignmentVerifier.getBestScore(read1.data, read2.data, Wgsim.getLowEnd(read1.idStr), Wgsim.getHighEnd(read1.idStr), params.maxDist))
        }
        case RichPairMultiHit(bestPos1, bestPos2, bestScore, secondBestPos1, secondBestPos2, secondBestScore, isRC) => {
          println("Align at bestPos pair: " + AlignmentVerifier.getBestScore(read1.data, read2.data, bestPos1, bestPos2, params.maxDist))
          println("Align at secondBestPos pair: " + AlignmentVerifier.getBestScore(read1.data, read2.data, secondBestPos1, secondBestPos2, params.maxDist))
          println("Align at wgsim pair: " + AlignmentVerifier.getBestScore(read1.data, read2.data, Wgsim.getLowEnd(read1.idStr), Wgsim.getHighEnd(read1.idStr), params.maxDist))
        }
        case _ =>
      }
      println
    })
  }

  def samEntryToAlignResult(e: RichSamEntry): AlignResult = {
    // check to see whether this read was mapped or not
    if (e.readIsUnmapped) NotFound
    else RichSingleHit(e.readAbsPos, e.readIsRC, e.readEditDist)
  }
  
  def needToCheckForMismatch(res1: AlignResult, res2: AlignResult, res3: AlignResult) = {
    res1.isInstanceOf[RichSingleHit] && res2.isInstanceOf[RichSingleHit] && res3.isInstanceOf[RichPairSingleHit]
  }
  
  def mismatch(res1: RichSingleHit, res2: RichSingleHit, res3: RichSingleHit, maxDist: Int) = {
    !((abs(res1.position - res2.position) <= maxDist) && (abs(res1.position - res3.position) <= maxDist) && (abs(res2.position - res3.position) <= maxDist))
  }

  def includePair(bfaRes: RichPairSingleHit, minSeparation: Int, maxSeparation: Int): Boolean = {
    val sep = abs(bfaRes.pos1 - bfaRes.pos2)
    (sep >= minSeparation && sep <= maxSeparation)
  }
}
