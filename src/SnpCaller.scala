package biggie

import scala.math.{max, min}

class SnpCaller(samFile: String, refFile: String, regionStart: Int = 1, regionEnd: Int = 2) {
  // [regionStart, regionEnd) with respect to reference, 1-indexed

  // Range of coverages for which to call bases, both in total and per direction.
  val TOTAL_COVERAGE_RANGE = 20 to 100
  val DIR_COVERAGE_RANGE = 8 to 50

  // What fraction of coverage the second-most-common base should have for us to report it;
  // you'd think this should be something like .4, but PCR imbalance messes it up.
  val SECOND_THRESHOLD = 0.2
  val SECOND_DIRECTIONAL_THRESHOLD = 0.01  // Ditto but per direction

  val ref: Array[Byte] = FASTA.read(refFile).pieces(0).data // 0-indexed
  val reader = new SamRegionReader(samFile, regionStart until regionEnd) // 1-indexed
  val regionSize = regionEnd - regionStart
  val baseCount = Array.ofDim[Int](2, 4, regionSize + 100)
  val coverage = Array.ofDim[Int](2, regionSize + 100)
  val snps = new Array[SNP](regionSize + 100)

  def run() {
    for (read <- reader.reads) {
      val dir = read.direction
      var posInRef = read.position
      var posInRead = 0
      // TODO: stop processing after end of region
      for ((count, op) <- read.parseCigar()) {
        op match {
          case 'M' =>
            for (i <- 0 until count) {
              // TODO: Filter based on Phred score
              val base = DNA.BASE_TO_CODE(read.sequence.charAt(posInRead))
              if ( (regionStart until regionEnd) contains posInRef) {
                val regionPos = posInRef - regionStart
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

          case other =>
            //println("Unhandled CIGAR element: " + other)
        }
      }
    }

    for (pos <- regionStart until regionEnd) {
      call(pos)
    }
    snps
  }

  // pos 0 is first base of region which is 1-indexed
  private def call(pos: Int) {
    val regionPos = pos - regionStart
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
    val regionPos = pos - regionStart
    if (base != ref(pos-1).toChar) {
      val snp = new SNP(ref(pos-1).toChar, base, base)
      snps(regionPos) = snp
      println("Calling SNP " + snp + " at " + (pos-1))
    } else {
      snps(regionPos) = SNP.MATCHING(base)
    }
  }

  private def callTwo(pos: Int, base1: Char, base2: Char) {
    val regionPos = pos - regionStart
    val snp = {
      if (base1 == ref(pos-1).toChar) {
        new SNP(ref(pos-1).toChar, base1, base2)
      } else if (base2 == ref(pos-1)) {
        new SNP(ref(pos-1).toChar, base2, base1)
      } else {
        new SNP(ref(pos-1).toChar, min(base1, base2).toChar, max(base1, base2).toChar)
      }
    }
    snps(regionPos) = snp
    //diffList += ((pos, snp))
    println("Calling SNP " + snp + " at " + (pos-1))
  }
}

object SnpCaller {
  def main(args: Array[String]) {
    val snps = new SnpCaller(args(0), args(1), args(2).toInt, args(3).toInt).run()
  }
}
