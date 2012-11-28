package snap.apps

import snap._

/**
 * Checks whether a set of reads in a FASTQ file align to a FASTA file,
 * without figuring out an exact position for each read. This is useful
 * for filtering out, say, human reads from a sample.
 */
object Recognizer {
  def main(args: Array[String]) : Unit = {
    args match {
      case Array(fastaFile, fastqFile, numSeeds, seedLen, maxTries, maxDist) =>
        run(fastaFile, fastqFile, numSeeds.toInt, seedLen.toInt, maxTries.toInt, maxDist.toInt)
      case _ =>
        System.err.println(
            "Usage: Recognizer <fastaFile> <fastqFile> <numSeeds> <seedLen> <maxTries> <maxDist>")
        System.exit(1)
    }
  }
  
  def run(fastaFile: String, fastqFile: String, numSeeds: Int, seedLen: Int, maxTries: Int, maxDist: Int) {
    println("Loading genome...")
    val genome = FASTA.read(fastaFile)
    println("Loading reads...")
    val reads = FASTQ.read(fastqFile)
    println("Indexing genome...")
    val idxStart = System.currentTimeMillis
    val builder = new IntHashIndexBuilder(seedLen, genome.totalSize)
    genome.addToIndex(builder)
    val index = builder.build()
    printf("Indexing took %.3fs%n", (System.currentTimeMillis - idxStart) / 1000.0)
    System.gc()
    println("Aligning reads...")
    val alnStart = System.currentTimeMillis
    val aligner = numSeeds match {
      case 1 => new SingleSeedAligner(genome, index, maxTries, maxDist, 1, 40, false)
      case 2 => new TwoSeedAligner(genome, index, maxTries, maxDist, 1, 40, false)
      case _ => throw new IllegalArgumentException("Invalid numSeeds: " + numSeeds)
    }
    val matched = FASTQ.writer("matched.fq")
    val unmatched = FASTQ.writer("unmatched.fq")
    var numMatched = 0
    var numUnmatched = 0
    for (i <- 0 until reads.size) {
      if (i % 1000000 == 0 && i > 0)
        printf("Progress: %d/%d%n", i, reads.size) 
      val read = reads(i)
      if (aligner.align(read.data).found) {
        matched.write(read)
        numMatched += 1
      } else {
        unmatched.write(read)
        numUnmatched += 1
      }
    }
    matched.close()
    unmatched.close()
    val alnTime = (System.currentTimeMillis - alnStart) / 1000.0
    printf("Aligning took %.3fs (%.0f reads/s)%n", alnTime, reads.size / alnTime)
    val total = numMatched + numUnmatched
    printf("Total matched: %d/%d%n", numMatched, total)
    printf("Total unmatched: %d/%d%n", numUnmatched, total)
  }
}
