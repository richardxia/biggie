package snap

abstract class Aligner {
  def align(read: Array[Byte]): AlignResult

  def align(read: String): AlignResult = align(read.getBytes)
}
