package snap

class Score (maxDist: Int, confDiff: Int, initScore: Int = Int.MaxValue, initPos: Long = 0L) {
  var bestScore = initScore
  var bestPos1 = initPos
  var bestPos2 = initPos
  
  var secondBestScore = initScore
  var secondBestPos1 = initPos
  var secondBestPos2 = initPos
  
  var firstReadIsRC = false
  
  def update(newScore: Int, pos1: Long, pos2: Long, isRC1: Boolean, debug: Boolean = false) = {
    if (newScore >= 0) {
      // Check that this isn't the same pair with slightly different offsets
      if (math.abs(pos1 - bestPos1) <= maxDist && math.abs(pos2 - bestPos2) <= maxDist) {
        // Found an alignment at nearly the same position, likely due to indels; if it's better than best,
        // then have it replace best, but otherwise don't count it as second-best because it's probably
        // just happening because two of our seeds had an indel in-between them
        if (newScore < bestScore) {
          if (debug) println("Old best: " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")

          bestScore = newScore
          bestPos1 = pos1
          bestPos2 = pos2
          firstReadIsRC = isRC1
        
          if (debug) println("Best got updated via first case to " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
        } 
      } else if (newScore < bestScore) {
        if (debug) println("Old best: " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
        if (debug) println("Old 2nd best: " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")

        // Update the best score and count the old best score as second-best
        secondBestScore = bestScore
        secondBestPos1 = bestPos1
        secondBestPos2 = bestPos2

        bestScore = newScore
        bestPos1 = pos1
        bestPos2 = pos2
        firstReadIsRC = isRC1

        if (debug) println("Best got updated via 2nd case to " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
        if (debug) println("2nd best got updated via 2nd case to " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")
      } else if (newScore < secondBestScore) {
        if (debug) println("Old 2nd best: " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")

        secondBestScore = newScore
        secondBestPos1 = pos1
        secondBestPos2 = pos2

        if (debug) println("2nd best got updated via 3rd case to " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")
      }
    }
  }

  def getResult: AlignResult = {
    if (bestScore <= maxDist) {
      if (secondBestScore >= bestScore + confDiff)
        return RichPairSingleHit(bestPos1, bestPos2, bestScore, firstReadIsRC)
      else
        return RichPairMultiHit(bestPos1, bestPos2, bestScore, secondBestPos1, secondBestPos2, secondBestScore, firstReadIsRC)
    } else {
      return NotFound
    }
  }

  def getDistToCheck(snapUpperBound: Int) = {
    if (bestScore > snapUpperBound) snapUpperBound + confDiff - 1
    else if (secondBestScore < bestScore + confDiff) bestScore - 1 // We're ambiguous, so need to find a better unambiguous hit
    else bestScore + confDiff - 1 // No need to search for worse scores than this
  }
}