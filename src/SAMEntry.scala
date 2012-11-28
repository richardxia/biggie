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
}
