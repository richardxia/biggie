package snap

import scala.math._
import spark.SparkContext
import SparkContext._

object AlignmentVerifier {
  def getBestScore(read1: Array[Byte], read2: Array[Byte], pos1: Long, pos2: Long, maxDist: Int): Int = {
    val lv = new LandauVishkin(maxDist)
    val scores = scala.collection.mutable.Set[Int]()
    var score1 = 0
    var score2 = 0

    // try all combinations
    // read1 forward at first pos, read2 reverse at second pos
    score1 = alignRead(read1, false, pos1, maxDist)
    score2 = alignRead(read2, true, pos2, maxDist)
    println("score1: " + score1)
    println("score2: " + score2)
    if (score1 >= 0 && score2 >= 0) scores += score1 + score2

    // read1 reverse at first pos, read2 forward at second pos
    score1 = alignRead(read1, true, pos1, maxDist)
    score2 = alignRead(read2, false, pos2, maxDist)
    println("score1: " + score1)
    println("score2: " + score2)
    if (score1 >= 0 && score2 >= 0) scores += score1 + score2
    
    // read1 forward at second pos, read2 reverse at first pos
    score1 = alignRead(read1, false, pos2, maxDist)
    score2 = alignRead(read2, true, pos1, maxDist)
    println("score1: " + score1)
    println("score2: " + score2)
    if (score1 >= 0 && score2 >= 0) scores += score1 + score2
    
    // read1 reverse at second pos, read2 forward at first pos
    score1 = alignRead(read1, true, pos2, maxDist)
    score2 = alignRead(read2, false, pos1, maxDist)
    println("score1: " + score1)
    println("score2: " + score2)
    if (score1 >= 0 && score2 >= 0) scores += score1 + score2

    println(scores.toList.mkString(", "))

    if (scores.size > 0)
      scores.min
    else 
      -1
  }
  
  def alignRead(read: Array[Byte], isRC: Boolean, pos: Long, maxDist: Int): Int = {
    val lv = new LandauVishkin(maxDist)
    val readLen = read.length
    val ref = new Array[Byte](readLen)
    GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)

    val readToCheck = 
      if (isRC) {
        val rc = new Array[Byte](readLen)
        DNA.getReverseComplement(read, rc)
        rc
      } else read

      println("ref:         " + new String(ref))
      println("readToCheck: " + new String(readToCheck))
    
    lv.distance(ref, readLen, readToCheck, readLen, maxDist)
  }

  def checkForMatch(bfaRes: RichPairSingleHit, alignerRes: RichSingleHit, checkFirstRead: Boolean, maxDist: Int): Boolean = {
    if (checkFirstRead)
      abs(bfaRes.pos1 - alignerRes.position) <= maxDist
    else
      abs(bfaRes.pos2 - alignerRes.position) <= maxDist
  }

  def checkForMatch(bfaRes: RichPairSingleHit, alignerRes: RichSingleHit, checkFirstRead: Boolean, maxDist: Int, numClippedBases: Int): Boolean = {
    val alignerPos = alignerRes.position - numClippedBases
    
    if (checkFirstRead)
      abs(bfaRes.pos1 - alignerPos) <= maxDist
    else
      abs(bfaRes.pos2 - alignerPos) <= maxDist
  }
  
  def needToVerify(res1: AlignResult, res2: AlignResult, maxSeparation: Int): Boolean = {
    val needToVerify = if (res1.isInstanceOf[RichSingleHit] && res2.isInstanceOf[RichSingleHit]) {
      // def getLocation(pos: Long): (String, Long) = {
      val (piece1, pos1) = GenomeLoader.genome.getLocation(res1.asInstanceOf[RichSingleHit].position)
      val (piece2, pos2) = GenomeLoader.genome.getLocation(res2.asInstanceOf[RichSingleHit].position)
      (piece1 == piece2 && abs(pos1 - pos2) <= maxSeparation)
    } else false
    needToVerify
  }

  def verifyAlignerError(read1: Array[Byte], read2: Array[Byte], bfaRes: RichPairSingleHit, alignerRes: RichPairSingleHit, maxDist: Int): Boolean = {
    println("bfa:")
    val bfaScore = getBestScore(read1, read2, bfaRes.pos1, bfaRes.pos2, maxDist)
    println("bfa score: " + bfaScore)
    println
    
    println("aligner:")
    val alignerScore = getBestScore(read1, read2, alignerRes.pos1 - 1, alignerRes.pos2 - 1, maxDist)  // recall that aligner is 1-indexed => needs to be corrected
    println("aligner score: " + alignerScore)
    println

    assert(alignerScore != bfaScore)
    //assert(bfaScore != -1)
    
    // aligner & bfa both got a hit, but aligner hit is worse || aligner didn't get a hit
    (bfaScore != -1 && alignerScore > bfaScore) || alignerScore == -1
  }
}