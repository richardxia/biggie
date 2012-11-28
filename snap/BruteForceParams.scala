package snap

import spark.SparkContext
import SparkContext._

class PairedBruteForceParams (val maxDist: Int, val confDiff: Int, separationDist: Int) extends Serializable {
  def getUpperBound = maxDist + confDiff - 1

  def getSnapSummary(snapRes1: AlignResult, snapRes2: AlignResult): (Int, Boolean) = {
    var upperBound = maxDist
    var firstIsRC = false

    snapRes1 match {
      case RichSingleHit(pos1, isRC1, editDist1) => 
        snapRes2 match {
          case RichSingleHit(pos2, isRC2, editDist2) => 
            if (scala.math.abs(pos1 - pos2) < separationDist && isRC1 != isRC2) {
              if (editDist1 + editDist2 < upperBound)
                upperBound = editDist1 + editDist2
              firstIsRC = isRC1
            }
          case _ =>
        }
      case _ =>
    }
    
    (upperBound, firstIsRC)
  }

  // should be absolute pos
  def getRead2StartPos(piece: GenomePiece, pos1: Long) = {
    if (pos1 - separationDist > piece.startIndex)
      pos1 - separationDist
    else
      piece.startIndex
  }
  
  // should be absolute pos
  def getRead2EndPos(piece: GenomePiece, pos1: Long) = {
    if (pos1 + separationDist < piece.endIndex)
      pos1 + separationDist
    else
      piece.endIndex
  }
}