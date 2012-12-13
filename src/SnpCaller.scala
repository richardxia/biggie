package biggie

import java.io.File

import scala.io.Source
import scala.math.{max, min}
import scala.collection.JavaConversions._

import net.sf.samtools.SAMFileReader
import net.sf.picard.reference.ReferenceSequence
import net.sf.picard.reference.ReferenceSequenceFile
import net.sf.picard.reference.IndexedFastaSequenceFile

class SnpCaller(bamFile: SAMFileReader, ref: ReferenceSequence, refSeq: String, region: Range, weirdness: Array[Float]) {
  // [region.start, region.end] with respect to reference, 1-indexed

  // Range of coverages for which to call bases, both in total and per direction.
  val TOTAL_COVERAGE_RANGE = 22 to 100
  val DIR_COVERAGE_RANGE = 10 to 50

  // What fraction of coverage the second-most-common base should have for us to report it;
  // you'd think this should be something like .4, but PCR imbalance messes it up.
  val SECOND_THRESHOLD = 0.2
  val SECOND_DIRECTIONAL_THRESHOLD = 0.01  // Ditto but per direction

  val WEIRDNESS_THRESHOLD = 0.1

  val baseCount = Array.ofDim[Int](2, 4, region.size + 100)
  val coverage = Array.ofDim[Int](2, region.size + 100)
  val snps = new Array[SNP](region.size + 100)

  val refBases = ref.getBases()

  def run() {
    val reads = bamFile.queryOverlapping(refSeq, region.start, region.end)
    for (read <- reads) {
      if (!read.getReadUnmappedFlag()) {
        val dir = if (read.getReadNegativeStrandFlag()) 1 else 0
        var posInRef = read.getAlignmentStart()
        var posInRead = 0
        //println("read " + read.getReadName() + " " + read.getAlignmentStart())
        // TODO: stop processing after end of region
        for (cigar <- read.getCigar().getCigarElements()) {
          val op = cigar.getOperator().toString()(0)
          val count = cigar.getLength()
          //println("CIGAR " + op + " " + count)
          op match {
            case 'M' =>
              for (i <- 0 until count) {
                // TODO: Filter based on Phred score
                val base = DNA.BASE_TO_CODE(read.getReadBases()(posInRead).asInstanceOf[Char])
                if (region contains posInRef) {
                  val regionPos = posInRef - region.start
                  //println("baseCount(" + dir + ")(" + base + ")(" + regionPos + ") += 1")
                  //println("coverage(" + dir + ")(" + regionPos + ") += 1")
                  baseCount(dir)(base)(regionPos) += 1
                  coverage(dir)(regionPos) += 1
                }
                posInRead += 1
                posInRef += 1
              }

            case 'I' =>
              posInRead += count

            case 'D' =>
              posInRef += count

            case 'S' =>
              posInRead += count

            case other =>
              System.err.println("Weird CIGAR: " + other)
              throw new RuntimeException("Weird CIGAR: " + other)
              //println("Unhandled CIGAR element: " + other)
          }
        }
      }
    }

    // Need to close BAMFileIterator before starting next one
    reads.close()

    if (weirdness != null) {
      for (pos <- region) {
        if (weirdness(pos) < WEIRDNESS_THRESHOLD)
          call(pos)
      }
    } else {
      for (pos <- region) {
        call(pos)
      }
    }
    snps
  }

  // pos 0 is first base of region which is 1-indexed
  private def call(pos: Int) {
    val regionPos = pos - region.start
    val totalCoverage = coverage(0)(regionPos) + coverage(1)(regionPos)
    if (!TOTAL_COVERAGE_RANGE.contains(totalCoverage) ||
        !DIR_COVERAGE_RANGE.contains(coverage(0)(regionPos)) ||
        !DIR_COVERAGE_RANGE.contains(coverage(1)(regionPos)) ) {
        //println("Skipping " + pos + "; TOTAL_COVERAGE: " + totalCoverage + "; DIR_COVERAGE: " + coverage(0)(regionPos) + ", " + coverage(1)(regionPos))
        return
    }
    val baseFracs = (0 until 4).map(base => {
        val frac = if (baseCount(0)(base)(regionPos) == 0 || baseCount(1)(base)(regionPos) == 0) {
          0.0
        } else {
          (baseCount(0)(base)(regionPos) + baseCount(1)(base)(regionPos)).toDouble / totalCoverage
        }
        (base, frac)
      }).sortWith(_._2 > _._2)
    val (best, bestFrac) = baseFracs(0)
    val (second, secondFrac) = baseFracs(1)

    //println(pos + " " + totalCoverage)
    //println(baseCount(0)(0)(regionPos) + " " + baseCount(0)(1)(regionPos) + " " + baseCount(0)(2)(regionPos) + " " + baseCount(0)(3)(regionPos))
    //println(baseCount(1)(0)(regionPos) + " " + baseCount(1)(1)(regionPos) + " " + baseCount(1)(2)(regionPos) + " " + baseCount(1)(3)(regionPos))
    //println(DNA.CODE_TO_BASE(best) + ": " + bestFrac + "; " + DNA.CODE_TO_BASE(second) + ": " + secondFrac)

    // Just call best for now
    val base1 = if (bestFrac == 0.0) 'N' else DNA.CODE_TO_BASE(best)
    val base2 = if (bestFrac == 0.0) 'N' else DNA.CODE_TO_BASE(second)
    // Now we need to decide whether to call one or both; let's use some arbitrary thresholds
    val secondMinDirFrac = min(baseCount(0)(second)(regionPos).toDouble / coverage(0)(regionPos),
                               baseCount(1)(second)(regionPos).toDouble / coverage(1)(regionPos))
    if (secondFrac >= SECOND_THRESHOLD && secondMinDirFrac >= SECOND_DIRECTIONAL_THRESHOLD) {
      callTwo(pos, base1, base2)
    } else {
      callOne(pos, base1)
    }
  }

  private def callOne(pos: Int, base: Char) {
    val regionPos = pos - region.start
    val refBase: Char = Character.toUpperCase(refBases(pos-1).toChar)
    if (base != refBase) {
      val snp = new SNP(refBase, base, base)
      snps(regionPos) = snp
      println("Calling SNP " + snp + " at " + (pos-1))
    } else {
      snps(regionPos) = SNP.MATCHING(base)
    }
  }

  private def callTwo(pos: Int, base1: Char, base2: Char) {
    val regionPos = pos - region.start
    val snp = {
      val refBase: Char = Character.toUpperCase(refBases(pos-1).toChar)
      if (base1 == refBase) {
        new SNP(refBase, base1, base2)
      } else if (base2 == refBase) {
        new SNP(refBase, base2, base1)
      } else {
        new SNP(refBase, min(base1, base2).toChar, max(base1, base2).toChar)
      }
    }
    snps(regionPos) = snp
    //diffList += ((pos, snp))
    println("Calling SNP " + snp + " at " + (pos-1))
  }
}

object SnpCaller {
  def readRef(refFileName: String): ReferenceSequenceFile = {
    new IndexedFastaSequenceFile(new File(refFileName))
  }

  def readWeirdness(weirdnessFileName: String, size: Int) = {
    val t1 = System.currentTimeMillis()
    val weirdness = new Array[Float](size + 100)
    for (line <- Source.fromFile(weirdnessFileName).getLines) {
      val Array(junk, pos, junk2, score) = line.split(' ')
      weirdness(pos.toInt) = score.toFloat
    }
    val t2 = System.currentTimeMillis()
    println("Weirdness loading: " + (t2 - t1))
    weirdness
  }

  def readRegions(regionsFileName: String) = {
    Source.fromFile(regionsFileName).getLines.map( line => {
      val Array(refSeq, start, end) = line.split('\t')
      (refSeq, start.toInt to end.toInt)
    })
  }

  def main(args: Array[String]) {
    if (args.size != 3 && args.size != 4) {
      println("Usage: SnpCaller alignments.bam reference.fa regions.txt [weirdness.txt]")
      println("We expect an alignments.bam.bai file which is located in the same directory as alignments.bam")
      println("regions.txt should be a text file with tab-delimited reference sequence name, start, and end positions on each line")
    }

    val bamFileName = args(0)
    val refFileName = args(1)
    val regionsFileName = args(2)

    val ref = readRef(refFileName)
    // This doesn't work if you have multiple sequences/chromosomes
    //val weirdness = if (args.size == 4) readWeirdness(arg(3), ref.size) else null
    val weirdness = null
    val regions = readRegions(regionsFileName)

    val bamFile = new SAMFileReader(new File(bamFileName), new File(bamFileName + ".bai"))
    bamFile.setValidationStringency(SAMFileReader.ValidationStringency.SILENT)

    for ((refSeq, region) <- regions) {
      val snps = new SnpCaller(bamFile, ref.getSequence(refSeq), refSeq, region, weirdness).run()
    }
  }
}
