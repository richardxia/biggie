package snap.apps

import snap._
import scala.math.{max, min}

object SimulatedReadTester {
  def main(args : Array[String]) : Unit = {
    args match {
      case Array(faFile, fqFile, numSeeds, seedLen, seedsToTry, maxDist, confDiff, maxHits) =>
        run(faFile, fqFile, parseSeq(numSeeds), seedLen.toInt, parseSeq(seedsToTry),
            parseSeq(maxDist), parseSeq(confDiff), parseSeq(maxHits))
      case _ =>
        System.err.println(
            "Usage: SimulatedReadTester faFile fqFile numSeeds seedLen seedsToTry maxDist confDiff maxHits")
        System.exit(1)
    }
  }

  // Parse a command-line argument representing a sequence of ints (e.g. 1-10 or 2,3,4)
  def parseSeq(str: String): Seq[Int] = {
    if (str.indexOf(',') != -1) {        // comma-separated list
      str.split(",").map(_.toInt)
    } else if (str.indexOf("-") != -1) { // range
      str.split("-") match {
        case Array(start, end) => start.toInt to end.toInt
        case _ => throw new IllegalArgumentException("Invalid range format: " + str)
      }
    } else { // single number
      Seq(str.toInt)
    }
  }

  def run(faFile: String, fqFile: String, numSeeds: Seq[Int], seedLen: Int, seedsToTry: Seq[Int],
          maxDist: Seq[Int], confDiff: Seq[Int], maxHits: Seq[Int])
  {
    println("Loading genome...")
    val genome = FASTA.read(faFile)
    println("Loading reads...")
    val reads = FASTQ.read(fqFile, false)
    println("Indexing genome...")
    val idxStart = System.currentTimeMillis
    val builder = new IntHashIndexBuilder(seedLen, genome.totalSize)
    genome.addToIndex(builder)
    val index = builder.build()
    printf("Indexing took %.3fs%n", (System.currentTimeMillis - idxStart) / 1000.0)
    val results = new Array[AlignResult](reads.size)
    for (ns <- numSeeds; stt <- seedsToTry; md <- maxDist; cd <- confDiff; mh <- maxHits) {
      System.gc()
      println()
      test(genome, index, reads, ns, stt, md, cd, mh)
    }
  }

  def test(genome: Genome, index: Index, reads: IndexedSeq[Read], numSeeds: Int, seedsToTry: Int,
           maxDist: Int, confDiff: Int, maxHits: Int)
  {
    printf("Trying seedLen=%d numSeeds=%d seedsToTry=%d maxDist=%d confDiff=%d maxHits=%d%n",
      index.seedLen, numSeeds, seedsToTry, maxDist, confDiff, maxHits)
    //println("Aligning reads...")
    val alnStart = System.currentTimeMillis
    val aligner = numSeeds match {
      case 1 => new SingleSeedAligner(genome, index, seedsToTry, maxDist, confDiff, maxHits)
      case 2 => new TwoSeedAligner(genome, index, seedsToTry, maxDist, confDiff, maxHits)
      case _ => throw new IllegalArgumentException("Invalid numSeeds: " + numSeeds)
    }
    val results = new Array[AlignResult](reads.size)
    var i = 0
    while (i < reads.size) {
      if (i % 1000000 == 0 && i > 0)
        printf("Progress: %d/%d%n", i, reads.size)
      results(i) = aligner.align(reads(i).data)
      i += 1
    }
    val alnTime = (System.currentTimeMillis - alnStart) / 1000.0
    printf("Aligning took %.3fs (%.0f reads/s)%n", alnTime, reads.size / alnTime)
    //printf("Checking alignments...%n")
    var good = 0       // Total reads that are not nonsense (too many N's)
    var confident = 0  // Good reads for which we report a single hit
    var correct = 0    // Confident hits that we placed in the right position
    var multiHits = 0  // Good reads for which we report multiple hits
    var unmatched = 0  // Good reads which we didn't match
    val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r
    for (i <- 0 until reads.size) {
      if (i % 1000000 == 0 && i > 0)
        printf("Progress: %d/%d%n", i, reads.size) 
      val read = reads(i)
      if (read.data.count(_ == 'N') < read.data.length / 2) {
        good += 1
        results(i) match {
          case SingleHit(pos, rc) =>
            confident += 1
            // Get the true position of the read and compare it
            read.idStr match {
	          case ID_REGEX(piece, start, end) =>
	            val lowEnd = min(start.toLong, end.toLong)
	            val highEnd = max(start.toLong, end.toLong)
	            val (myPiece, myOffset) = genome.getLocation(pos)
	            if (myPiece == piece && myOffset >= lowEnd - maxDist && myOffset <= highEnd + maxDist)
	              correct += 1
	          case _ =>
	            println("WARNING: Read ID " + read.idStr + " is not of the expected form! Are you using wgsim?")
	        }
          case MultipleHits =>
            multiHits += 1
          case NotFound =>
            unmatched += 1
        }
      }
    }
    val errors = confident - correct
    val matched = good - unmatched
    printf("Good reads: %d/%d (%.2f%%)%n", good, reads.size, good * 100.0 / reads.size)
    printf("Total matched: %d/%d (%.2f%%)%n", matched, good, matched  * 100.0/ good)
    printf("Confident matches: %d/%d (%.2f%%)%n", confident, good, confident * 100.0 / good)
    printf("Alignment errors: %d/%d (%.2f%%)%n", errors, confident, errors * 100.0 / confident) 
    aligner match {
      case s: SingleSeedAligner => s.printStats()
      case _ =>
    }
  }
}
