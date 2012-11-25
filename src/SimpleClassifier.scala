import scala.collection.mutable.ArrayBuffer
import scala.math.{max, min}
import scala.io.Source
import java.util.concurrent.ConcurrentHashMap
import scala.collection.mutable.ArrayBuffer

object SimpleClassifier {
  def main(args: Array[String]) {
    println("Loading sam file...")
    new SimpleClassifier(args(0),args(1).toInt).run()
  }
}

class SimpleClassifier(samFile: String, length: Int){
  
  val BASE_TO_CODE = new Array[Byte](128)
  BASE_TO_CODE('A') = 0
  BASE_TO_CODE('C') = 1
  BASE_TO_CODE('G') = 2
  BASE_TO_CODE('T') = 3

  val WEIRDNESS_WINDOW = 65       // Should be an odd number so it's balanced around our base.
  val SUB_WEIRDNESS = 3
  val INDEL_WEIRDNESS = 10
  val MULTI_WEIRDNESS = 3         // From a multiple-hits read mapping to a given place.
  val WEIRDNESS_THRESHOLD = 0.05
  val MIN_PHRED = 30              // Ignore bases mapped with score lower than this.
  val HIGH_COMPLEXITY_GAP = 5
  val MIN_HIGH_COMPLEXITY_REGION_LENGTH = 10

  // Range of coverages for which to call bases, both in total and per direction.
  val MIN_TOTAL_COVERAGE = 22
  val MAX_TOTAL_COVERAGE = 10000
  val MIN_DIR_COVERAGE = 10
  val MAX_DIR_COVERAGE = 5000

  // At what quality do we consider reads confident
  val READ_QUALITY_THRESHOLD = 20

  val bufSize = length + 400

  val subCount = new Array[Int](bufSize)   // TODO: Could use Shorts here
  val insCount = new Array[Int](bufSize)
  val delCount = new Array[Int](bufSize)
  val coverage = Array.ofDim[Int](2, bufSize)     // (direction, location)
  val baseCount = Array.ofDim[Int](2, 4, bufSize) // (direction, base, location)
  val multiCount = new Array[Int](bufSize)
  var offset = -1000
  var last_pos = -1

  def run() 
  {
    parseSam()
    callUpTo(last_pos)
  }
  
  def parseSam()
  {
	// TODO: Maintain a set of reads at each position and eliminate duplicates
	val reads = SAM.read(samFile)
	for (read <- reads) 
	{
		if(offset == -1000)
		{
			offset = read.position - 200
		}
		val readPos = read.position - offset
		val dir = read.direction
		if(readPos - 200 > length)
		{
			println(readPos+" "+length);
			return
		}

		last_pos = max(last_pos,readPos)


		// Decide whether to shouldIgnore the read (e.g. if it maps to another chromosome)
		if (shouldIgnoreRead(read)) {
			  //logDebug("Will not use " + read + " for calling")
			  if (read.mapQuality < READ_QUALITY_THRESHOLD) {
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
			  for ((count, op) <- parseCigar(read.cigar)) 
			  {
				    op match 
				    {
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
    var total_length = 0
    while (pos < position) {
	val totalCoverage = coverage(0)(pos) + coverage(1)(pos)
	if (computeWeirdness(pos) >= WEIRDNESS_THRESHOLD &&
	totalCoverage >= MIN_TOTAL_COVERAGE && totalCoverage <= MAX_TOTAL_COVERAGE &&
	coverage(0)(pos) >= MIN_DIR_COVERAGE && coverage(0)(pos) <= MAX_DIR_COVERAGE &&
	coverage(1)(pos) >= MIN_DIR_COVERAGE && coverage(1)(pos) <= MAX_DIR_COVERAGE) 
	{
		if(left == -1)
		{
			left = pos;
		}
		else if(pos - last > HIGH_COMPLEXITY_GAP)
		{
			if( last-left+1 >= MIN_HIGH_COMPLEXITY_REGION_LENGTH )
			{
				println("High complexity region of length " + (last-left+1) + " in ["+(left+offset)+", "+(last+offset)+"]")				
				total_length += last-left+1
			}
			left = pos;
		}
		last = pos;
	}
	pos += 1
    }
    println("Classifying genome sequence of length "+(position-200)+" in ["+(200+offset)+", "+(position+offset)+"]")
    println("High complexity length = "+total_length)
  }

  // Call the base at the given position
  def call(pos: Int) {
    
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
    if (read.mapQuality < READ_QUALITY_THRESHOLD) {
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
}

object Utils {
  /** Split a string around instances of a given delimiter */
  def split(s: String, delimiter: Char): Seq[String] = {
    val buf = new ArrayBuffer[String]
    var i = 0
    while (i < s.length) {
      var j = i
      while (j < s.length && s.charAt(j) != delimiter) {
        j += 1
      }
      if (j > i) {
        buf += s.substring(i, j);
      }
      i = j
      while (i < s.length && s.charAt(i) == delimiter) {
        i += 1
      }
    }
    return buf
  }

  def parsePhred(score: Char): Int = score - 33
  def parsePhred(score: Byte): Int = score - 33
}

class SAMEntry(
    val readId: String,
    val flags: Int,
    val piece: String, 
    val position: Int,
    val mapQuality: Int,
    val cigar: String,
    val nextPiece: String,
    val nextPosition: Int,
    val templateLen: Int,
    val sequence: String,
    val quality: String) {
  override def toString(): String = readId

  def direction = if ((flags & SAM.REVERSE) != 0) 1 else 0
}

object SAM {
  def parseEntry(line: String): SAMEntry = {
    val fields = Utils.split(line, '\t')
    new SAMEntry(
      fields(0),                  // read ID
      fields(1).toInt,            // flags
      fields(2),                  // piece
      fields(3).toInt,            // position
      fields(4).toInt,            // map quality
      fields(5),                  // cigar
      fields(6),                  // next piece
      fields(7).toInt,            // next position
      if (fields(2) == fields(6)) fields(8).toInt else 0, // template len
      fields(9),                  // sequence
      fields(10)                  // quality
    )
  }

  def read(file: String): Iterator[SAMEntry] = {
    val lines = Source.fromFile(file).getLines
    lines.filter(!_.startsWith("@")).map(parseEntry)
  }

  // Bits in the flags field
  val REVERSE = 0x10
}
