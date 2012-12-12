package biggie

import java.io.File

import scala.collection.JavaConversions._
import scala.math.{max, min}

import java.io.FileOutputStream
import java.nio.ByteBuffer
import java.nio.FloatBuffer

import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMRecord

class SimpleClassifier(bamFile: SAMFileReader){
  
  val BASE_TO_CODE = new Array[Byte](128)
  BASE_TO_CODE('A') = 0
  BASE_TO_CODE('C') = 1
  BASE_TO_CODE('G') = 2
  BASE_TO_CODE('T') = 3

  val WEIRDNESS_WINDOW = 65       // Should be an odd number so it's balanced around our base.
  val SUB_WEIRDNESS = 3
  val INDEL_WEIRDNESS = 10
  val MULTI_WEIRDNESS = 3         // From a multiple-hits read mapping to a given place.
  val WEIRDNESS_THRESHOLD = 0.1
  val MIN_PHRED = 30              // Ignore bases mapped with score lower than this.
  val MIN_HIGH_COMPLEXITY_REGION_LENGTH = 1000
  val MIN_HIGH_COMPLEXITY_REGION_DENSITY = 0.15

  // Range of coverages for which to call bases, both in total and per direction.
  val MIN_TOTAL_COVERAGE = 22
  val MAX_TOTAL_COVERAGE = 10000
  val MIN_DIR_COVERAGE = 10
  val MAX_DIR_COVERAGE = 5000

  // At what quality do we consider reads confident
  val READ_QUALITY_THRESHOLD = 20

  val reads = bamFile

  val bufSize = reads.getFileHeader().getSequence(0).getSequenceLength() + 400

  val subCount = new Array[Int](bufSize)   // TODO: Could use Shorts here
  val insCount = new Array[Int](bufSize)
  val delCount = new Array[Int](bufSize)
  val coverage = Array.ofDim[Int](2, bufSize)     // (direction, location)
  val baseCount = Array.ofDim[Int](2, 4, bufSize) // (direction, base, location)
  val multiCount = new Array[Int](bufSize)
  val weirdnessBuf = new Array[Float](bufSize)
  var offset = -1000
  var last_pos = -1

  def run(): Array[Float] = {
    parseSam()
    callUpTo(last_pos)
    //writeSerializedWeirdness()
    weirdnessBuf
  }

  def writeSerializedWeirdness() {
    val floatBuf = FloatBuffer.wrap(weirdnessBuf)
    val byteBuf = ByteBuffer.allocate(4 * bufSize)
    byteBuf.asFloatBuffer.put(floatBuf)
    //val weirdnesses = ByteBuffer.wrap(weirdness)
    val fos = new FileOutputStream("weirdness.bin")
    val channel = fos.getChannel()
    channel.write(byteBuf)
    fos.close()
  }
  
  def parseSam()
  {
    // TODO: Maintain a set of reads at each position and eliminate duplicates
    for (read <- reads) {
      if(offset == -1000) {
        offset = read.getAlignmentStart() - 200
      }
      val readPos = read.getAlignmentStart() - offset
      val dir = if (read.getReadNegativeStrandFlag()) 1 else 0

      last_pos = max(last_pos,readPos)


      // Decide whether to shouldIgnore the read (e.g. if it maps to another chromosome)
      if (shouldIgnoreRead(read)) {
        //logDebug("Will not use " + read + " for calling")
        if (read.getMappingQuality() < READ_QUALITY_THRESHOLD) {
          //logDebug("Using it to increment multiRead count though")
          val len = getCoveredLength(read)
          var i = 0
          while (i < len) {
            multiCount(readPos + i) += 1
            i += 1
          }
        }
      } 
      else 
      {
        //logDebug("Processing " + read + " @" + read.position + ": " + read.cigar + "\n" + read.sequence)
        // Read the CIGAR string and update the base counts with it
        var posInRead = 0
        var posInRef = readPos
        for ((count, op) <- parseCigar(read.getCigarString())) {
          op match 
          {
            case 'S' =>
              posInRead += count

            case '=' =>
              var i = 0
              while (i < count) {
                val base = BASE_TO_CODE(read.getReadBases()(posInRead+i).asInstanceOf[Char])
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
                if (read.getBaseQualities()(posInRead + i) >= MIN_PHRED) {
                  val base = BASE_TO_CODE(read.getReadBases()(posInRead + i).asInstanceOf[Char])
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
                if (read.getBaseQualities()(posInRead + i) >= MIN_PHRED) {
                  val base = BASE_TO_CODE(read.getReadBases()(posInRead + i).asInstanceOf[Char])
                  if (base != 'N') {
                    baseCount(dir)(base)(posInRef + i) += 1
                    coverage(dir)(posInRef + i) += 1
                    //***if (read.sequence.charAt(posInRead + i) != ref(posInRef + i)) {
                      //***subCount(posInRef + i) += 1
                      //***}
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
                println("Unhandled CIGAR element: " + other)
                //logError("Unhandled CIGAR element: " + other)
            }
          }
        }
      }
    }

  // Call bases up to, but not including, the given position
  def callUpTo(position: Int) {
    var left = -1
    var last = -1
    var pos = 200
    var counter = 0
    var total_length = 0
    var covered_count = 0
    var total_count = 0
    while (pos < position) {
      val totalCoverage = coverage(0)(pos) + coverage(1)(pos)
      weirdnessBuf(pos+offset) = computeWeirdness(pos)
      if (computeWeirdness(pos) >= WEIRDNESS_THRESHOLD &&
        totalCoverage >= MIN_TOTAL_COVERAGE && totalCoverage <= MAX_TOTAL_COVERAGE &&
        coverage(0)(pos) >= MIN_DIR_COVERAGE && coverage(0)(pos) <= MAX_DIR_COVERAGE &&
        coverage(1)(pos) >= MIN_DIR_COVERAGE && coverage(1)(pos) <= MAX_DIR_COVERAGE) 
      {
        if(left == -1) {
          left = pos
        }
        else if( (counter+1)/(pos-left+1.0) <  MIN_HIGH_COMPLEXITY_REGION_DENSITY)
        {
          if( last-left+1 >= MIN_HIGH_COMPLEXITY_REGION_LENGTH ) {
            println("Region:\t"+(left+offset)+"\t--\t"+(last+offset)+"\tLength:\t"+(last-left+1)+"\tDensity:\t"+(counter/(last-left+1.0)))
            total_length += last-left+1
            covered_count += counter
          }
          counter = 0
          left = pos
        }
        last = pos
        counter += 1
        total_count += 1
      }
      pos += 1
    }

    if( last-left+1 >= MIN_HIGH_COMPLEXITY_REGION_LENGTH )
      {
      println("Region:\t"+(left+offset)+"\t--\t"+(last+offset)+"\tLength:\t"+(last-left+1)+"\tDensity:\t"+(counter/(last-left+1.0)))
      total_length += last-left+1
    }

    println("Classifying genome sequence of length "+(position-200)+" in ["+(200+offset)+", "+(position+offset)+"]")
    println("High complexity length = "+total_length)
    println("Covers "+covered_count+"/"+total_count+" weird bases")
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

  def shouldIgnoreRead(read: SAMRecord): Boolean = {
    if (read.getMappingQuality() < READ_QUALITY_THRESHOLD) {
      return true
    }
    var amountClipped = 0
    for ((count, op) <- parseCigar(read.getCigarString()) if op == 'S') {
      amountClipped += count
    }
    if (read.getReadLength() - amountClipped < 70) {
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
  def getCoveredLength(read: SAMRecord): Int = {
    parseCigar(read.getCigarString()).filter(p => "=XMD".contains(p._2)).map(_._1).sum
  }

  def totalBaseCount(base: Int, pos: Int): Int = baseCount(0)(base)(pos) + baseCount(1)(base)(pos)
}

object SimpleClassifier {
  def main(args: Array[String]) {
    println("Loading bam file...")
    val bamFileName = args(0)
    val bamFile = new SAMFileReader(new File(bamFileName), new File(bamFileName + ".bai"))
    new SimpleClassifier(bamFile).run()
  }
}
