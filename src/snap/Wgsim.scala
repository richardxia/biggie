package snap

import scala.math._

object Wgsim {
  val ID_REGEX = """@*([^:]+)_(\d+)_(\d+)_\d+:.*""".r
  
  def getChr(id: String): String = {
    var chr = ""
    id match {
      case ID_REGEX(piece, start, end) => chr = piece
    }
    chr
  }
  
  def getLowEnd(id: String): Long = {
    var lowEnd = 0L
    id match {
      case ID_REGEX(piece, start, end) => lowEnd = min(start.toLong, end.toLong)
    }
    lowEnd
  }
  
  def getHighEnd(id: String): Long = {
    var highEnd = 0L
    id match {
      case ID_REGEX(piece, start, end) => highEnd = max(start.toLong, end.toLong)
    }
    highEnd
  }
  
  def isCorrect(id: String, pos: Long, maxDist: Int, zeroIndexing: Boolean = true): Boolean = {
    val posToCheck = 
      if (zeroIndexing) pos + 1
      else pos
    
    val (myPiece, myOffset) = GenomeLoader.genome.getLocation(posToCheck)
    
    (myPiece == getChr(id) && myOffset >= getLowEnd(id) - maxDist && myOffset <= getHighEnd(id) + maxDist)
  }
  
  def isCorrect(e: SAMEntry): Boolean = {
    (e.piece == getChr(e.readId) && e.position >= getLowEnd(e.readId) && e.position <= getHighEnd(e.readId))
  }
  
  def getWgsimId(originalId: String, pos1: Long, pos2: Long, zeroIndexing: Boolean = true): String = {
    val (myPiece1, myOffset1) = GenomeLoader.genome.getLocation(pos1)
    val (myPiece2, myOffset2) = GenomeLoader.genome.getLocation(pos2)
    
    val lowEnd = 
      if (zeroIndexing)
        min(myOffset1, myOffset2) + 1
      else
        min(myOffset1, myOffset2)
      
    val highEnd = 
      if (zeroIndexing)
        max(myOffset1, myOffset2) + 1
      else
        max(myOffset1, myOffset2)
    
    if (myPiece1 == myPiece2)
      List(myPiece1, lowEnd, highEnd, 0 + ":" + originalId).mkString("_")
    else
      null
  }
}