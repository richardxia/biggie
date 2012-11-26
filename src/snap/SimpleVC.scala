package snap

import scala.collection.mutable.ArrayBuffer
import scala.math.{max, min}

import DNA.{BASE_TO_CODE, CODE_TO_BASE}

/**
 * Pretty much the simplest possible variant caller, which calls only SNPs in regions
 * where the "weirdness" is not too high. Works on a single reference piece and a
 * sorted SAM file.
 */
class SimpleVC(reference: GenomePiece, samFile: String) extends Logging {
  // Weirdness parameters. We sum up the "weirdness" in a window around each base, defined as
  // a weighted count of substitutions and indels seen in that window divided by the coverage
  // of the window, and don't call bases in locations that exceed WEIRDNESS_THRESHOLD.
  val WEIRDNESS_WINDOW = 65       // Should be an odd number so it's balanced around our base.
  val SUB_WEIRDNESS = 3
  val INDEL_WEIRDNESS = 10
  val MULTI_WEIRDNESS = 3         // From a multiple-hits read mapping to a given place.
  val WEIRDNESS_THRESHOLD = 0.1
  val MIN_PHRED = 30              // Ignore bases mapped with score lower than this.

  // Range of coverages for which to call bases, both in total and per direction.
  val MIN_TOTAL_COVERAGE = 22
  val MAX_TOTAL_COVERAGE = 100
  val MIN_DIR_COVERAGE = 10
  val MAX_DIR_COVERAGE = 50

  // What fraction of coverage the second-most-common base should have for us to report it;
  // you'd think this should be something like .4, but PCR imbalance messes it up.
  val SECOND_THRESHOLD = 0.2
  val SECOND_DIRECTIONAL_THRESHOLD = 0.01  // Ditto but per direction

  // At what quality do we consider reads confident
  val READ_QUALITY_THRESHOLD = 20

  val bufSize = reference.data.length + 400

  val ref: Array[Byte] = reference.data
  val subCount = new Array[Int](bufSize)   // TODO: Could use Shorts here
  val insCount = new Array[Int](bufSize)
  val delCount = new Array[Int](bufSize)
  val coverage = Array.ofDim[Int](2, bufSize)     // (direction, location)
  val baseCount = Array.ofDim[Int](2, 4, bufSize) // (direction, base, location)
  val multiCount = new Array[Int](bufSize)
  val weirdness = new Array[Float](bufSize)
  val snps = Array.fill(bufSize)(SNP.UNCALLED)
  val diffList = new ArrayBuffer[(Int, SNP)]

  var nextToCall = 0  // Next position to make a call at

  def run(low: Int = 0, high: Int = reference.data.length) {
    // TODO: Maintain a set of reads at each position and eliminate duplicates
    val reads = SAM.read(samFile)
    for (read <- reads) {
      val readPos = read.position
      val dir = read.direction
      if (readPos > low) {
        if (readPos > high + WEIRDNESS_WINDOW) {
          callUpTo(high + WEIRDNESS_WINDOW + 1)
          return
        }
        callUpTo(readPos)
        // Decide whether to shouldIgnore the read (e.g. if it maps to another chromosome)
        if (shouldIgnoreRead(read)) {
          logDebug("Will not use " + read + " for calling")
          if (read.piece == reference.name && read.mapQuality < READ_QUALITY_THRESHOLD) {
            logDebug("Using it to increment multiRead count though")
            val len = getCoveredLength(read)
            var i = 0
            while (i < len) {
              multiCount(readPos + i) += 1
              i += 1
            }
          }
        } else {
          logDebug("Processing " + read + " @" + read.position + ": " + read.cigar + "\n" + read.sequence)
          // Read the CIGAR string and update the base counts with it
          var posInRead = 0
          var posInRef = readPos
          for ((count, op) <- parseCigar(read.cigar)) {
            op match {
              case 'S' =>
                posInRead += count

              case '=' =>
                var i = 0
                while (i < count) {
                  val base = BASE_TO_CODE(read.sequence.charAt(posInRead + i))
                  if (base != 'N') {
                    baseCount(dir)(base)(posInRef + i) += 1
                    coverage(dir)(posInRef + i) += 1
                  }
                  i += 1
                }
                posInRead += count
                posInRef += count

              case 'X' =>
                var i = 0
                while (i < count) {
                  if (Utils.parsePhred(read.quality.charAt(posInRead + i)) >= MIN_PHRED) {
                    val base = BASE_TO_CODE(read.sequence.charAt(posInRead + i))
                    if (base != 'N') {
                      baseCount(dir)(base)(posInRef + i) += 1
                      coverage(dir)(posInRef + i) += 1
                      subCount(posInRef + i) += 1
                    }
                  }
                  i += 1
                }
                posInRead += count
                posInRef += count

              case 'M' =>
                var i = 0
                while (i < count) {
                  if (Utils.parsePhred(read.quality.charAt(posInRead + i)) >= MIN_PHRED) {
                    val base = BASE_TO_CODE(read.sequence.charAt(posInRead + i))
                    if (base != 'N') {
                      baseCount(dir)(base)(posInRef + i) += 1
                      coverage(dir)(posInRef + i) += 1
                      if (read.sequence.charAt(posInRead + i) != ref(posInRef + i)) {
                        subCount(posInRef + i) += 1
                      }
                    }
                  }
                  i += 1
                }
                posInRead += count
                posInRef += count

              case 'I' =>
                posInRead += count
                insCount(posInRef) += count

              case 'D' =>
                var i = 0
                while (i < count) {
                  delCount(posInRef + i) += 1
                  i += 1
                }
                posInRef += count

              case other =>
                logError("Unhandled CIGAR element: " + other)
            }
          }
        }
      }
    }
    callUpTo(high + WEIRDNESS_WINDOW + 1)
  }

  // Call bases up to, but not including, the given position
  def callUpTo(position: Int) {
    while (nextToCall < min(position, bufSize)) {
      call(nextToCall)
      nextToCall += 1
    }
  }

  // Call the base at the given position
  def call(pos: Int) {
    weirdness(pos) = computeWeirdness(pos)
    val totalCoverage = coverage(0)(pos) + coverage(1)(pos)
    if (weirdness(pos) < WEIRDNESS_THRESHOLD &&
        totalCoverage >= MIN_TOTAL_COVERAGE && totalCoverage <= MAX_TOTAL_COVERAGE &&
        coverage(0)(pos) >= MIN_DIR_COVERAGE && coverage(0)(pos) <= MAX_DIR_COVERAGE &&
        coverage(1)(pos) >= MIN_DIR_COVERAGE && coverage(1)(pos) <= MAX_DIR_COVERAGE) {
      // Find the best and second-best possible bases at this position
      var best = 0
      var bestFrac = 0.0
      var second = 0
      var secondFrac = 0.0
      var i = 0
      while (i < 4) {
        val frac = if (baseCount(0)(i)(pos) == 0 || baseCount(1)(i)(pos) == 0) {
          0.0
        } else {
          (baseCount(0)(i)(pos) + baseCount(1)(i)(pos)).toDouble / totalCoverage
        }
        if (frac > bestFrac) {
          second = best
          secondFrac = bestFrac
          best = i
          bestFrac = frac
        } else if (frac > secondFrac) {
          second = i
          secondFrac = frac
        }
        i += 1
      }
      // Now we need to decide whether to call one or both; let's use some arbitrary thresholds
      val secondMinDirFrac = min(baseCount(0)(second)(pos).toDouble / coverage(0)(pos),
                                 baseCount(1)(second)(pos).toDouble / coverage(1)(pos))
      if (secondFrac >= SECOND_THRESHOLD && secondMinDirFrac >= SECOND_DIRECTIONAL_THRESHOLD) {
        callTwo(pos, CODE_TO_BASE(best), CODE_TO_BASE(second))
      } else {
        callSingle(pos, CODE_TO_BASE(best))
      }
    }
  }

  def callSingle(pos: Int, base: Char) {
    if (base != ref(pos).toChar) {
      val snp = new SNP(ref(pos).toChar, base, base)
      snps(pos) = snp
      diffList += ((pos, snp))
      logInfo("Calling SNP " + snp + " at " + pos)
    } else {
      snps(pos) = SNP.MATCHING(base)
    }
  }

  def callTwo(pos: Int, base1: Char, base2: Char) {
    val snp = {
      if (base1 == ref(pos).toChar) {
        new SNP(ref(pos).toChar, base1, base2)
      } else if (base2 == ref(pos)) {
        new SNP(ref(pos).toChar, base2, base1)
      } else {
        new SNP(ref(pos).toChar, min(base1, base2).toChar, max(base1, base2).toChar)
      }
    }
    snps(pos) = snp
    diffList += ((pos, snp))
    logInfo("Calling SNP " + snp + " at " + pos)
  }

  def computeWeirdness(pos: Int): Float = {
    var weirdness = 0.0f
    var totalCoverage = 0.0f
    var i = max(0, pos - WEIRDNESS_WINDOW / 2)
    while (i <= min(pos + WEIRDNESS_WINDOW / 2, bufSize - 1)) {
      if (i != pos) {  // Don't count the substitution weirdness of the purported SNP itself
        weirdness += SUB_WEIRDNESS * subCount(i)
      }
      weirdness += INDEL_WEIRDNESS * (insCount(i) + delCount(i))
      weirdness += MULTI_WEIRDNESS * multiCount(i)
      totalCoverage += coverage(0)(i) + coverage(1)(i)
      i += 1
    }
    return weirdness / totalCoverage
  }

  def shouldIgnoreRead(read: SAMEntry): Boolean = {
    if (read.mapQuality < READ_QUALITY_THRESHOLD || read.piece != reference.name) {
      return true
    }
    var amountClipped = 0
    for ((count, op) <- parseCigar(read.cigar) if op == 'S') {
      amountClipped += count
    }
    if (read.sequence.length - amountClipped < 70) {
      return true
    }
    return false
  }

  def parseCigar(cigar: String): Iterator[(Int, Char)] = new Iterator[(Int, Char)] {
    var pos = 0

    def hasNext: Boolean = pos < cigar.length

    override def next(): (Int, Char) = {
      if (!hasNext) {
        throw new java.util.NoSuchElementException("next on empty iterator")
      }
      var count = 0
      while (cigar.charAt(pos) >= '0' && cigar.charAt(pos) <= '9') {
        count = 10 * count + (cigar.charAt(pos) - '0')
        pos += 1
      }
      val op = cigar.charAt(pos)
      pos += 1
      return (count, op)
    }
  }

  // Get the length of the reference genome that the given read spans
  def getCoveredLength(read: SAMEntry): Int = {
    parseCigar(read.cigar).filter(p => "=XMD".contains(p._2)).map(_._1).sum
  }

  def totalBaseCount(base: Int, pos: Int): Int = baseCount(0)(base)(pos) + baseCount(1)(base)(pos)

  def print(start: Int, end: Int) {
    println(" pos:" + (start until end by 4).map(" %-14d|".format(_)).mkString)
    //println("-----" + "----" * (end - start))
    println("   A:" + (start until end).map(i => " %3d".format(totalBaseCount(DNA.A, i))).mkString)
    println("   C:" + (start until end).map(i => " %3d".format(totalBaseCount(DNA.C, i))).mkString)
    println("   G:" + (start until end).map(i => " %3d".format(totalBaseCount(DNA.G, i))).mkString)
    println("   T:" + (start until end).map(i => " %3d".format(totalBaseCount(DNA.T, i))).mkString)
    println(" ref:" + ref.slice(start, end).map(x => " %3c".format(x.toChar)).mkString)
    println(" cal:" + snps.slice(start, end).map(x => " %3s".format(x.diffString)).mkString)
    println(" ind:" + (start until end).map(i => " %3d".format(insCount(i) + delCount(i))).mkString)
    println(" wrd:" + weirdness.slice(start, end).map(x => " %3d".format((x * 100).toInt)).mkString)
  }

  def printAround(pos: Int, radius: Int = 8) {
    print(max(pos - radius, 0), min(pos + radius, bufSize))
  }
}

object SimpleVC {
  def main(args: Array[String]) {
    println("Loading genome...")
    val genome = FASTA.read(args(0))
    new SimpleVC(genome.pieces(0), args(1)).run()
  }
}
