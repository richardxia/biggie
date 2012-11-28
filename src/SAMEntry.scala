package biggie

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
  val endPosition = position + sequence.size // one position past the last base

  def direction = if ((flags & SAM.REVERSE) != 0) 1 else 0

  def parseCigar(): Iterator[(Int, Char)] = new Iterator[(Int, Char)] {
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
}
