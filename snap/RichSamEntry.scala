package snap

import scala.math._

class RichSamEntry(samLine: String) {
  val fields = samLine.split("\t")
  def ID = fields(0)
  def flags = fields(1).toInt
  def readPiece = fields(2)
  def readPos = fields(3).toLong
  def readAbsPos = GenomeLoader.genome.getAbsPos(readPiece, readPos)
  def confident(threshold: Int = 10) = {
    // check for *
    val mappingQualityStr = fields(4)
    val confident = if (mappingQualityStr != "*") {
      val mappingQuality = mappingQualityStr.toInt
      mappingQuality >= threshold
    } else
      false
    confident
  }
  def pairPiece =
    if (fields(6) == "=") readPiece
    else fields(6)
  def pairPos = fields(7).toLong
  def pairAbsPos = GenomeLoader.genome.getAbsPos(pairPiece, pairPos)
  def pairSeparation = highEnd - lowEnd
  def lowEnd = min(readPos, pairPos)
  def highEnd = max(readPos, pairPos)
  def read = new Read(fields(0).getBytes, fields(9).getBytes, fields(10).getBytes)
  def cigarStr = fields(5)
  def readIsRC = {
    /*
    val reverse = flags & 0x10
    val forward = flags & 0x20
    
    if (forward != 0) false
    else if (reverse != 0) true
    else throw new IllegalArgumentException("Flags error")
    */
	  (flags & 0x10) != 0
  }
  def readIsUnmapped = {
    val unmapped = flags & 0x4
    
    unmapped != 0
  }
  def pairIsUnmapped = pairPiece == "*"
  def illuminaID = {
    val eID = ID.split("/")(0)
    val i = eID.indexOf(":")
    eID.substring(i + 1, eID.length)
  }
  def illuminaIDWithSuffix = {
    val i = ID.indexOf(":")
    ID.substring(i + 1, ID.length)
  }
  def readEditDist: Int = {
    var editDist = 0
    (0 until cigarStr.length).map(i => {
      // look for any I's, X's, or D's
      // count # before
      if (cigarStr(i) == 'X' || cigarStr(i) == 'I' || cigarStr(i) == 'D')
        editDist += new String(Array(cigarStr(i-1))).toInt
    })
    
    editDist
  }
  def numPrefixBasesClipped: Int = {
    val CIGAR_REGEX = """(\d+)S(.*)""".r
    
    val numClipped = cigarStr match {
      case CIGAR_REGEX(prefix, rest) => prefix.toInt
      case _ => 0
    }
    
    numClipped
  }
}