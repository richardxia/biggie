package snap

sealed trait AlignResult {
  def found: Boolean
}

case class SingleHit(position: Long, isRC: Boolean) extends AlignResult {
  override def found = true
}

case object MultipleHits extends AlignResult {
  override def found = true
}

case object NotFound extends AlignResult {
  override def found = false
  override def toString = "NotFound"
}

case class RichSingleHit(position: Long, isRC: Boolean, editDistance: Int) extends AlignResult {
  override def found = true
  override def toString = {
    val (myPiece, myOffset) = GenomeLoader.genome.getLocation(position)
    
    List("RichSingleHit - bestPos: " + myPiece + "(" + myOffset + ")", "isRC: " + isRC, "bestScore: " + editDistance).mkString(", ")
  }
}

case class RichPairSingleHit(pos1: Long, pos2: Long, score: Int, isRC: Boolean) extends AlignResult {
  override def found = true
  override def toString = {
    val (myPiece1, myOffset1) = GenomeLoader.genome.getLocation(pos1)
    val (myPiece2, myOffset2) = GenomeLoader.genome.getLocation(pos2)
    
    List("RichPairSingleHit", "read1: " + myPiece1 + "(" + myOffset1 + ") - " + isRC, "read2: " + myPiece2 + "(" + myOffset2 + ") - " + !isRC, "score: " + score).mkString("\n")
  }
}

object RichPairSingleHitHelper {
  // note:  primarily for checking position; score is for pair, not just this read
  def getFirstHit(res: RichPairSingleHit) = RichSingleHit(res.pos1, res.isRC, res.score)
  
  def getSecondHit(res: RichPairSingleHit) = RichSingleHit(res.pos2, !res.isRC, res.score)
  
  def getPairHit(res1: RichSingleHit, res2: RichSingleHit) = {
    assert(res1.isRC != res2.isRC)
    RichPairSingleHit(res1.position, res2.position, res1.editDistance + res2.editDistance, res1.isRC)
  }
}

case class RichMultipleHits(bestPos: Long, bestScore: Int, isRC: Boolean, secondBestPos: Long, secondBestScore: Int) extends AlignResult {
  override def found = true
  override def toString = {
    val (myPiece1, myOffset1) = GenomeLoader.genome.getLocation(bestPos)
    val (myPiece2, myOffset2) = GenomeLoader.genome.getLocation(secondBestPos)

    List("RichMultipleHit - bestPos: " + myPiece1 + "(" + myOffset1 + ")", "isRC: " + isRC, "bestScore: " + bestScore, 
      "secondBestPos: " + myPiece2 + "(" + myOffset2 + ")", "secondBestScore: " + secondBestScore).mkString(", ")
  }
}

case class RichPairMultiHit(bestPos1: Long, bestPos2: Long, bestScore: Int, secondBestPos1: Long, secondBestPos2: Long, secondBestScore: Int, 
  isRC: Boolean) extends AlignResult {
  override def found = true
  override def toString = {
    val (myPiece1b, myOffset1b) = GenomeLoader.genome.getLocation(bestPos1)
    val (myPiece2b, myOffset2b) = GenomeLoader.genome.getLocation(bestPos2)
    
    val (myPiece1s, myOffset1s) = GenomeLoader.genome.getLocation(secondBestPos1)
    val (myPiece2s, myOffset2s) = GenomeLoader.genome.getLocation(secondBestPos2)
    
    List("RichPairMultipleHit", "Best: " + myPiece1b + "(" + myOffset1b + ", " + myOffset2b + ") at " + bestScore, 
      "Second Best: " + myPiece1s + "(" + myOffset1s + ", " + myOffset2s + ") at " 
      + secondBestScore).mkString("\n")
  }
} 


case class RichManyHits(numHits: Int, maxDist: Int) extends AlignResult {
  import scala.collection.mutable.{Map, Set}

  // Keep track of a bunch of hits per read
  // To be used for paired-end alignment
  override def found = true
  var worstEditDist = maxDist
  var matches = Map[Int, Set[Long]]()
  var numMatches = 0

  private def addMatch(pos: Long, editDist: Int) = {
    if (matches.keySet.contains(editDist)) {
      val s = matches.get(editDist)
      s match {
        case Some(posSet) => {
          matches += ((editDist, posSet + pos))
        }
        case _ =>
      }
    } else {
      matches += ((editDist, Set(pos)))
    }
  }
  
  def addNewMatch(pos: Long, editDist: Int) = {
    if (editDist < maxDist) {
      if (numMatches == 0) {
        addMatch(pos, editDist)
        worstEditDist = editDist
        numMatches += 1
      } else if (numMatches < numHits)  {
        addMatch(pos, editDist)
        if (editDist > worstEditDist)
          worstEditDist = editDist
        numMatches += 1
      } else {
        if (editDist < worstEditDist) {
          // add this tuple
          addMatch(pos, editDist)

          // remove previous worst & update worst dist
          matches -= worstEditDist
          worstEditDist = matches.keySet.max
        }
      }
    }
  }  
  override def toString = {
    "Matches (worst at edit distance " + worstEditDist + "): " + matchingPos.mkString(", ")
  }
  
  def matchingPos = matches.keySet.flatMap(k => 
    matches.get(k) match {
      case Some(s) => s.toList
      case None => List()
    }).toList.sortWith(_ < _)
}