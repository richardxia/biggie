package snap

import spark.SparkContext
import SparkContext._
import scala.math.{max, min}
import it.unimi.dsi.fastutil.longs.LongArrayList
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import scala.collection.mutable.Set

object BruteForceAligner {
  def snapAlign(index: Index, reads: IndexedSeq[Read], numSeeds: Int, seedsToTry: Int, maxDist: Int, confDiff: Int, maxHits: Int): Array[AlignResult] = {
    printf("Trying seedLen=%d numSeeds=%d seedsToTry=%d maxDist=%d confDiff=%d maxHits=%d%n",
      index.seedLen, numSeeds, seedsToTry, maxDist, confDiff, maxHits)
    val alnStart = System.currentTimeMillis
    val aligner = numSeeds match {
      case 1 => new SingleSeedAligner(GenomeLoader.genome, index, seedsToTry, maxDist, confDiff, maxHits)
      case 2 => new TwoSeedAligner(GenomeLoader.genome, index, seedsToTry, maxDist, confDiff, maxHits)
      case _ => throw new IllegalArgumentException("Invalid numSeeds: " + numSeeds)
    }
    val results = new Array[AlignResult](reads.size)
    var i = 0
    while (i < reads.size) {
      if (i % 1000000 == 0 && i > 0)
        printf("Progress: %d/%d%n", i, reads.size)
      results(i) = aligner.align(reads(i).data)
      i += 1
    }
    val alnTime = (System.currentTimeMillis - alnStart) / 1000.0
    printf("Aligning took %.3fs (%.0f reads/s)%n", alnTime, reads.size / alnTime)
    
    results
  }
  
  def align(read: Array[Byte], initDist: Int, maxDist: Int, confDiff: Int, snapRC: Boolean): AlignResult = {
    val genomeSize = GenomeLoader.genome.totalSize
    val readLen = read.length
    val lv = new LandauVishkin(maxDist + confDiff - 1)
    
    var pos = 0L
    var bestScore = Int.MaxValue
    var secondBestScore = Int.MaxValue
    var bestPos = 0L
    var secondBestPos = 0L
    var bestIsRC = snapRC
    
    // if snap says that the best alignment is using the reverse complement, try that first
    val readToCheck = if (snapRC) {
      val rc = new Array[Byte](readLen)
      DNA.getReverseComplement(read, rc)
      rc
    } else read
    
    for (p <- GenomeLoader.genome.pieces) {
      println("Aligning against " + p.name + "...")
      pos = p.startIndex
      while (pos < p.endIndex) {
        if (pos % 10000000 == 0)
          println("Pos: " + pos)
      
        val ref = new Array[Byte](readLen)
        GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)
      
        val distToCheck = {
          if (bestScore > initDist) initDist + confDiff - 1
          else if (secondBestScore < bestScore + confDiff) bestScore - 1 // We're ambiguous, so need to find a better unambiguous hit
          else bestScore + confDiff - 1 // No need to search for worse scores than this
        }
          
        val score = lv.distance(ref, readLen, readToCheck, readLen, distToCheck)
      
        if (score != -1 && score <= initDist + confDiff - 1) {
          if (math.abs(pos - bestPos) <= maxDist) {
            // Found an alignment at nearly the same position, likely due to indels; if it's better than best,
            // then have it replace best, but otherwise don't count it as second-best because it's probably
            // just happening because two of our seeds had an indel in-between them
            if (score < bestScore) {
              bestScore = score
              bestPos = pos
            } 
          } else if (score < bestScore) {
            // Update the best score and count the old best score as second-best
            secondBestScore = bestScore
            secondBestPos = bestPos
            bestPos = pos
            bestScore = score
          } else if (score < secondBestScore) {
            secondBestScore = score
            secondBestPos = pos
          }
        }
      
        pos += 1
      }
    
      println("Aligning rc...")
      val rc = new Array[Byte](readLen)
      DNA.getReverseComplement(readToCheck, rc)
    
      pos = p.startIndex
      while (pos < p.endIndex) {
        if (pos % 10000000 == 0)
          println("Pos: " + pos)
      
        val ref = new Array[Byte](readLen)
        GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)
      
        val distToCheck = {
          if (bestScore > initDist) initDist + confDiff - 1
          else if (secondBestScore < bestScore + confDiff) bestScore - 1 // We're ambiguous, so need to find a better unambiguous hit
          else bestScore + confDiff - 1 // No need to search for worse scores than this
        }
      
        val score = lv.distance(ref, readLen, rc, readLen, distToCheck)
      
        if (score != -1 && score <= initDist + confDiff - 1) {
          if (math.abs(pos - bestPos) <= maxDist) {
            // Found an alignment at nearly the same position, likely due to indels; if it's better than best,
            // then have it replace best, but otherwise don't count it as second-best because it's probably
            // just happening because two of our seeds had an indel in-between them
            if (score < bestScore) {
              bestScore = score
              bestPos = pos
              bestIsRC = !bestIsRC
            }
          } else if (score < bestScore) {
            // Update the best score and count the old best score as second-best
            secondBestScore = bestScore
            secondBestPos = bestPos
            bestPos = pos
            bestScore = score
            bestIsRC = !bestIsRC
          } else if (score < secondBestScore) {
            secondBestScore = score
            secondBestPos = pos
          }
        }
      
        pos += 1
      }
    }
        
    if (bestScore <= initDist) {
      if (secondBestScore >= bestScore + confDiff)
        return RichSingleHit(bestPos, bestIsRC, bestScore)
      else
        return RichMultipleHits(bestPos, bestScore, bestIsRC, secondBestPos, secondBestScore)
    } else {
      return NotFound
    }
  }

  def alignInRegion(read: Array[Byte], pieceName: String, pieceStartPos: Long, pieceEndPos: Long, initDist: Int, maxDist: Int, confDiff: Int, snapRC: Boolean = false): AlignResult = {
    val genomeSize = GenomeLoader.genome.totalSize
    val readLen = read.length
    val lv = new LandauVishkin(maxDist + confDiff - 1)
    
    var pos = 0L
    var bestScore = Int.MaxValue
    var secondBestScore = Int.MaxValue
    var bestPos = 0L
    var secondBestPos = 0L
    var bestIsRC = snapRC
    
    // if snap says that the best alignment is using the reverse complement, try that first
    val readToCheck = if (snapRC) {
      val rc = new Array[Byte](readLen)
      DNA.getReverseComplement(read, rc)
      rc
    } else read
    
    // figure out where to start
    val piece = GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId(pieceName))
        
    val startPos = piece.startIndex + pieceStartPos    
    val endPos = piece.startIndex + pieceEndPos
    
    pos = startPos
    while (pos < endPos) {
      if (pos % 10000000 == 0)
        println("Pos: " + pos)
    
      val ref = new Array[Byte](readLen)
      GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)
    
      val distToCheck = {
        if (bestScore > initDist) initDist + confDiff - 1
        else if (secondBestScore < bestScore + confDiff) bestScore - 1 // We're ambiguous, so need to find a better unambiguous hit
        else bestScore + confDiff - 1 // No need to search for worse scores than this
      }
        
      val score = lv.distance(ref, readLen, readToCheck, readLen, distToCheck)
    
      if (score != -1 && score <= initDist + confDiff - 1) {
        if (math.abs(pos - bestPos) <= maxDist) {
          // Found an alignment at nearly the same position, likely due to indels; if it's better than best,
          // then have it replace best, but otherwise don't count it as second-best because it's probably
          // just happening because two of our seeds had an indel in-between them
          if (score < bestScore) {
            bestScore = score
            bestPos = pos
          } 
        } else if (score < bestScore) {
          // Update the best score and count the old best score as second-best
          secondBestScore = bestScore
          secondBestPos = bestPos
          bestPos = pos
          bestScore = score
        } else if (score < secondBestScore) {
          secondBestScore = score
          secondBestPos = pos
        }
      }
    
      pos += 1
    }
  
    println("Aligning rc...")
    val rc = new Array[Byte](readLen)
    DNA.getReverseComplement(readToCheck, rc)
  
    pos = startPos
    while (pos < endPos) {
      if (pos % 10000000 == 0)
        println("Pos: " + pos)
    
      val ref = new Array[Byte](readLen)
      GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)
    
      val distToCheck = {
        if (bestScore > initDist) initDist + confDiff - 1
        else if (secondBestScore < bestScore + confDiff) bestScore - 1 // We're ambiguous, so need to find a better unambiguous hit
        else bestScore + confDiff - 1 // No need to search for worse scores than this
      }
    
      val score = lv.distance(ref, readLen, rc, readLen, distToCheck)
    
      if (score != -1 && score <= initDist + confDiff - 1) {
        if (math.abs(pos - bestPos) <= maxDist) {
          // Found an alignment at nearly the same position, likely due to indels; if it's better than best,
          // then have it replace best, but otherwise don't count it as second-best because it's probably
          // just happening because two of our seeds had an indel in-between them
          if (score < bestScore) {
            bestScore = score
            bestPos = pos
            bestIsRC = !bestIsRC
          }
        } else if (score < bestScore) {
          // Update the best score and count the old best score as second-best
          secondBestScore = bestScore
          secondBestPos = bestPos
          bestPos = pos
          bestScore = score
          bestIsRC = !bestIsRC
        } else if (score < secondBestScore) {
          secondBestScore = score
          secondBestPos = pos
        }
      }
    
      pos += 1
    }
        
    //if (bestScore <= initDist) {
      if (secondBestScore >= bestScore + confDiff)
        return RichSingleHit(bestPos, bestIsRC, bestScore)
      else
        return RichMultipleHits(bestPos, bestScore, bestIsRC, secondBestPos, secondBestScore)
    /*} else {
      else
      return NotFound
    }*/
  }

  def lookForLocalHit(lv: LandauVishkin, read: Array[Byte], centerPos: Long, separationDist: Int, maxDist: Int, upperBoundForRead2: Int, confDiff: Int, firstReadIsRC: Boolean, debug: Boolean = false): AlignResult = {
    val readLen = read.length

    // since one read in a pair will be reverse complemented and one won't, make this read the opposite of the first
    // if first read is reverse complemented, leave this one alone
    // else, reverse complement this one
    var readToCheck = 
    if (firstReadIsRC) {
      read
    } else {
      var rc = new Array[Byte](readLen)

      DNA.getReverseComplement(read, rc)
      rc
    }

    // check this read in the neighborhood of (pos - separationDist, pos + separationDist), using distToCheck as upper bound on edit distance
    // NOTE:  not exactly sure on the [] vs. (), whether I need a +/- 1 in there

    // figure out what piece pos is on
    val myPiece = GenomeLoader.genome.getPiece(centerPos)
    
    // assuming that reads won't span two chromosomes (truncate range if necessary)
    val startIndex = 
    if (centerPos - separationDist > myPiece.startIndex)
      centerPos - separationDist
    else
      myPiece.startIndex
      
    val endIndex = 
    if (centerPos + separationDist < myPiece.endIndex)
      centerPos + separationDist
    else
      myPiece.endIndex
    
    var pos = startIndex
    var bestScore = Int.MaxValue
    var secondBestScore = Int.MaxValue
    var bestPos = 0L
    var secondBestPos = 0L

    /*
    if (debug) {
      println("In lookForLocalHit...")
      println("centerPos: " + centerPos)
      println("startIndex: " + startIndex)
      println("endIndex: " + endIndex)
    }
    */
    
    while (pos < endIndex) {
      val distToCheck = {
        if (bestScore > upperBoundForRead2) upperBoundForRead2
        else if (secondBestScore < bestScore + confDiff) bestScore - 1 // We're ambiguous, so need to find a better unambiguous hit
        else bestScore + confDiff - 1 // No need to search for worse scores than this
      }
      
      val ref = new Array[Byte](readLen)
      GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)
            
      //val score = editDist(pos, read, lv, distToCheck, false)
      val score = lv.distance(ref, readLen, readToCheck, readLen, distToCheck)
      
      if (score != -1 && score <= upperBoundForRead2 + confDiff - 1) {
        if (math.abs(pos - bestPos) <= maxDist) {
          // Found an alignment at nearly the same position, likely due to indels; if it's better than best,
          // then have it replace best, but otherwise don't count it as second-best because it's probably
          // just happening because two of our seeds had an indel in-between them
          if (score < bestScore) {
            if (debug) println("  Local(" + centerPos + "): Single hit - Old best: " + bestScore + "(" + bestPos + ")")

            bestScore = score
            bestPos = pos
            
            if (debug) println("  Local(" + centerPos + "): Single hit - Best got updated via first case to " + bestScore + "(" + bestPos + ")")
          } 
        } else if (score < bestScore) {
          // Update the best score and count the old best score as second-best
          if (debug) println("  Local(" + centerPos + "): Single hit - Old best: " + bestScore + "(" + bestPos + ")")
          if (debug) println("  Local(" + centerPos + "): Single hit - Old 2nd best: " + secondBestScore + "(" + secondBestPos + ")")

          // Update the best score and count the old best score as second-best
          secondBestScore = bestScore
          secondBestPos = bestPos
          bestPos = pos
          bestScore = score

          if (debug) println("  Local(" + centerPos + "): Single hit - Best got updated via 2nd case to " + bestScore + "(" + bestPos + ")")
          if (debug) println("  Local(" + centerPos + "): Single hit - 2nd best got updated via 2nd case to " + secondBestScore + "(" + secondBestPos + ")")
        } else if (score < secondBestScore) {
          if (debug) println("  Local(" + centerPos + "): Single hit - Old 2nd best: " + secondBestScore + "(" + secondBestPos + ")")

          secondBestScore = score
          secondBestPos = pos

          if (debug) println("  Local(" + centerPos + "): Single hit - 2nd best got updated via 3rd case to " + secondBestScore + "(" + secondBestPos + ")")

        }
      }
            
      pos += 1
    }
    
    if (bestScore <= upperBoundForRead2) {
      if (secondBestScore >= bestScore + confDiff)
        return RichSingleHit(bestPos, !firstReadIsRC, bestScore)
      else
        return RichMultipleHits(bestPos, bestScore, !firstReadIsRC, secondBestPos, secondBestScore)
    } else {
      return NotFound
    }
  }

  def pairedAlignOneWay(p: GenomePiece, read1: Array[Byte], read2: Array[Byte], firstReadIsRC: Boolean, upperBound: Int, maxDist: Int, confDiff: Int, separationDist: Int, lv: LandauVishkin, initBestScore: Int, initBestPos1: Long, initBestPos2: Long, initSecondBestScore: Int, initSecondBestPos1: Long, initSecondBestPos2: Long, debug: Boolean = false) = {
    var bestScore = initBestScore
    var bestPos1 = initBestPos1
    var bestPos2 = initBestPos2
    var secondBestScore = initSecondBestScore
    var secondBestPos1 = initSecondBestPos1
    var secondBestPos2 = initSecondBestPos2

    val readLen = read1.length

    val dir = if (firstReadIsRC) "R" else "F"

    println("Aligning against " + p.name + "...")
    var pos = p.startIndex
    while (pos < p.endIndex) {
    //var pos = p.startIndex + 49694000
    //while (pos < p.startIndex + 49694500) {
      if (pos % 10000000 == 0)
        if (debug) println("Pos: " + pos)

      val ref = new Array[Byte](readLen)
      GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)
    
      val distToCheck = {
        if (bestScore > upperBound) upperBound + confDiff - 1
        else if (secondBestScore < bestScore + confDiff) bestScore - 1 // We're ambiguous, so need to find a better unambiguous hit
        else bestScore + confDiff - 1 // No need to search for worse scores than this
      }
        
      //val score1 = editDist(pos, readToCheck1, lv, distToCheck, false)
      val score1 = lv.distance(ref, readLen, read1, readLen, distToCheck)

      if (score1 != -1 && score1 <= upperBound + confDiff - 1) {
        // I found a pretty good hit for the 1st read, so I'm going to check nearby to see if I can find a good hit for the 2nd read
        val read2UpperBound = math.max(math.min(bestScore, upperBound) - score1, 0)
             // just distToCheck - score1?
        val res2 = lookForLocalHit(lv, read2, pos, separationDist, maxDist, read2UpperBound, confDiff, firstReadIsRC, debug)

        res2 match {
          case RichSingleHit(pos2, isRC, score2) => {
            val score = score1 + score2

            // Update best, 2nd best pos
            if (math.abs(pos - bestPos1) <= maxDist) {
              // Found an alignment at nearly the same position, likely due to indels; if it's better than best,
              // then have it replace best, but otherwise don't count it as second-best because it's probably
              // just happening because two of our seeds had an indel in-between them
              if (score < bestScore) {
                if (debug) println("Single hit (" + dir + ") - Old best: " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")

                bestScore = score
                bestPos1 = pos
                bestPos2 = pos2
                
                if (debug) println("Single hit (" + dir + ") - Best got updated via first case to " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
              } 
            } else if (score < bestScore) {
              if (debug) println("Single hit (" + dir + ") - Old best: " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
              if (debug) println("Single hit (" + dir + ") - Old 2nd best: " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")

              // Update the best score and count the old best score as second-best
              secondBestScore = bestScore
              secondBestPos1 = bestPos1
              secondBestPos2 = bestPos2

              bestScore = score
              bestPos1 = pos
              bestPos2 = pos2

              if (debug) println("Single hit (" + dir + ") - Best got updated via 2nd case to " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
              if (debug) println("Single hit (" + dir + ") - 2nd best got updated via 2nd case to " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")
            } else if (score < secondBestScore) {
              if (debug) println("Single hit (" + dir + ") - Old 2nd best: " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")

              secondBestScore = score
              secondBestPos1 = pos
              secondBestPos2 = pos2

              if (debug) println("Single hit (" + dir + ") - 2nd best got updated via 3rd case to " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")
            }
          }
          case RichMultipleHits(hitBestPos, hitBestScore, isRC, hitSecondBestPos, hitSecondBestScore) => {
            val scoreA = score1 + hitBestScore
            val scoreB = score1 + hitSecondBestScore
            
            // Update best, 2nd best pos
            if (math.abs(pos - bestPos1) <= maxDist) {
              // Found an alignment at nearly the same position, likely due to indels; if it's better than best,
              // then have it replace best, but otherwise don't count it as second-best because it's probably
              // just happening because two of our seeds had an indel in-between them
              if (scoreA < bestScore) {
                if (debug) println("Multi hit A (" + dir + ") - Old best: " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")

                bestScore = scoreA
                bestPos1 = pos
                bestPos2 = hitBestPos
                
                if (debug) println("Multi hit A (" + dir + ") - Best got updated via first case to " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
              } 
            } else if (scoreA < bestScore) {
              if (debug) println("Multi hit A (" + dir + ") - Old best: " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
              if (debug) println("Multi hit A (" + dir + ") - Old 2nd best: " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")

              // Update the best score and count the old best score as second-best
              secondBestScore = bestScore
              secondBestPos1 = bestPos1
              secondBestPos2 = bestPos2

              bestScore = scoreA
              bestPos1 = pos
              bestPos2 = hitBestPos

              if (debug) println("Multi hit A (" + dir + ") - Best got updated via 2nd case to " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
              if (debug) println("Multi hit A (" + dir + ") - 2nd best got updated via 2nd case to " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")
            } else if (scoreA < secondBestScore) {
              if (debug) println("Multi hit A (" + dir + ") - Old 2nd best: " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")

              secondBestScore = scoreA
              secondBestPos1 = pos
              secondBestPos2 = hitBestPos

              if (debug) println("Multi hit A (" + dir + ") - 2nd best got updated via 3rd case to " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")
            }
            
            // Update best, 2nd best pos
            if (math.abs(pos - bestPos1) <= maxDist) {
              // Found an alignment at nearly the same position, likely due to indels; if it's better than best,
              // then have it replace best, but otherwise don't count it as second-best because it's probably
              // just happening because two of our seeds had an indel in-between them
              if (scoreB < bestScore) {
                if (debug) println("Multi hit B (" + dir + ") - Old best: " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")

                bestScore = scoreB
                bestPos1 = pos
                bestPos2 = hitSecondBestPos
                
                if (debug) println("Multi hit B (" + dir + ") - Best got updated via first case to " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
              } 
            } else if (scoreB < bestScore) {
              if (debug) println("Multi hit B (" + dir + ") - Old best: " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
              if (debug) println("Multi hit B (" + dir + ") - Old 2nd best: " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")

              // Update the best score and count the old best score as second-best
              secondBestScore = bestScore
              secondBestPos1 = bestPos1
              secondBestPos2 = bestPos2

              bestScore = scoreB
              bestPos1 = pos
              bestPos2 = hitSecondBestPos

              if (debug) println("Multi hit B (" + dir + ") - Best got updated via 2nd case to " + bestScore + "(" + bestPos1 + "," + bestPos2 + ")")
              if (debug) println("Multi hit B (" + dir + ") - 2nd best got updated via 2nd case to " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")
            } else if (scoreB < secondBestScore) {
              if (debug) println("Multi hit B (" + dir + ") - Old 2nd best: " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")

              secondBestScore = scoreB
              secondBestPos1 = pos
              secondBestPos2 = hitSecondBestPos

              if (debug) println("Multi hit B (" + dir + ") - 2nd best got updated via 3rd case to " + secondBestScore + "(" + secondBestPos1 + "," + secondBestPos2 + ")")
            }
          }
          case _ =>
        }
      }
    
      pos += 1
    }
    
    (bestScore, bestPos1, bestPos2, secondBestScore, secondBestPos1, secondBestPos2)
  }

  def pairedAlign(read1: Array[Byte], read2: Array[Byte], snapUpperBound: Int, maxDist: Int, confDiff: Int, separationDist: Int, snapFirstReadIsRC: Boolean): AlignResult = {
    println("In pairedAlign...")
    println("snapFirstReadIsRC = " + snapFirstReadIsRC)
    
    val genomeSize = GenomeLoader.genome.totalSize
    val readLen = read1.length
    val lv = new LandauVishkin(2 * maxDist + confDiff - 1)
    
    var pos = 0L
    var bestScore = Int.MaxValue
    //var bestScore = 20
    var secondBestScore = Int.MaxValue
    var bestPos1 = 0L
    var bestPos2 = 0L
    var secondBestPos1 = 0L
    var secondBestPos2 = 0L
    var firstReadIsRC = snapFirstReadIsRC
    
    // if snap says that the best alignment is using the reverse complement, try that first
    var readToCheck1 = read1
    var readToCheck2 = read2
    
    if (firstReadIsRC) {
      var rc = new Array[Byte](readLen)

      DNA.getReverseComplement(read1, rc)
      readToCheck1 = rc
    }
    
    for (p <- GenomeLoader.genome.pieces) {
    //val pieces = List(GenomeLoader.genome.pieces(GenomeLoader.genome.getPieceId("chr5")))
    //for (p <- pieces) {
      // Align forward
      val (bestScoreF, bestPos1F, bestPos2F, secondBestScoreF, secondBestPos1F, secondBestPos2F) = pairedAlignOneWay(p, readToCheck1, readToCheck2, firstReadIsRC, snapUpperBound, maxDist, confDiff, separationDist, lv, bestScore, bestPos1, bestPos2, secondBestScore, secondBestPos1, secondBestPos2, true)

      // Align rc
      println("Aligning rc...")
      val rc1 = new Array[Byte](readLen)
      DNA.getReverseComplement(readToCheck1, rc1)

      val (bestScoreR, bestPos1R, bestPos2R, secondBestScoreR, secondBestPos1R, secondBestPos2R) = pairedAlignOneWay(p, rc1, readToCheck2, !firstReadIsRC, snapUpperBound, maxDist, confDiff, separationDist, lv, bestScoreF, bestPos1F, bestPos2F, secondBestScoreF, secondBestPos1F, secondBestPos2F, true)
    
      bestScore = bestScoreR
      bestPos1 = bestPos1R
      bestPos2 = bestPos2R
      secondBestScore = secondBestScoreR
      secondBestPos1 = secondBestPos1R
      secondBestPos2 = secondBestPos2R
    }

    if (bestScore <= snapUpperBound) {
      if (secondBestScore >= bestScore + confDiff)
        return RichPairSingleHit(bestPos1, bestPos2, bestScore, firstReadIsRC)
      else
        return RichPairMultiHit(bestPos1, bestPos2, bestScore, secondBestPos1, secondBestPos2, secondBestScore, firstReadIsRC)
    } else {
      return NotFound
    }
  }

  def parallelAlign(sc: SparkContext, outputDest: String, readsWithSnapResults: IndexedSeq[(Read, AlignResult)], maxDist: Int, confDiff: Int, outputPrefix: String = "bfaResults") = {
    sc.parallelize(readsWithSnapResults).map(r => {
      val (read, snapResult) = r
      
      val (initDist, snapRC) = 
        snapResult match {
          case RichSingleHit(pos, isRC, editDistance) => (min(editDistance, maxDist.toInt), isRC)
          case RichMultipleHits(bestPos, bestScore, isRC, secondBestPos, secondBestScore) => (min(bestScore, maxDist.toInt), isRC) // maybe secondBestScore?
          case _ => (maxDist.toInt, false)
        }

      (read, align(read.data, initDist, maxDist, confDiff, snapRC))
    })
    .saveAsObjectFile(outputDest + outputPrefix)
  }
  
  def parallelPairedAlign(sc: SparkContext, outputDest: String, outputPrefix: String, readsWithSnapResults: IndexedSeq[((Read, AlignResult), (Read, AlignResult))], maxDist: Int, confDiff: Int, separationDist: Int) = {
    // note that confDiff is per pair, not per read
    
    sc.parallelize(readsWithSnapResults).map(pair => {
      val ((read1, snapRes1), (read2, snapRes2)) = pair
      var snapFirstReadIsRC = false
      
      // set the upper bound for edit distance
      var snapUpperBound = 2 * maxDist
      var bestPairScore = verifyPairedHit(snapRes1, snapRes2, separationDist)
      if (bestPairScore < snapUpperBound)
        snapUpperBound = bestPairScore
      
      // align the reads
      ((read1, read2), pairedAlign(read1.data, read2.data, snapUpperBound, maxDist, confDiff, separationDist, snapFirstReadIsRC)) // weird naming; not clear what diff b/t upperBound & maxDist is
    })
    .saveAsObjectFile(outputDest + outputPrefix)
  }
  
  def verifyPairedSingleHits(hit1: RichSingleHit, hit2: RichSingleHit, separationDist: Int): Int = {
    val RichSingleHit(pos1, isRC1, editDistance1) = hit1
    val RichSingleHit(pos2, isRC2, editDistance2) = hit2

    var pairScore = 
      if (scala.math.abs(pos1 - pos2) < separationDist && isRC1 != isRC2)
        editDistance1 + editDistance2
      else
        Int.MaxValue

    pairScore
  }

  def verifyPairedMixedHits(hit1: RichMultipleHits, hit2: RichSingleHit, separationDist: Int): Int = {
    val RichMultipleHits(bestPos1, bestScore1, isRC1, secondBestPos1, secondBestScore1) = hit1
    val RichSingleHit(pos2, isRC2, score2) = hit2
    
    var bestPairScore = Int.MaxValue
    
    // check 2 hits for read1 with hit for read2
    val hit1Best = RichSingleHit(bestPos1, isRC1, bestScore1)
    val pairScoreA = verifyPairedSingleHits(hit1Best, hit2, separationDist)
    
    val hit1SecondBest = RichSingleHit(secondBestPos1, isRC1, secondBestScore1)
    val pairScoreB = verifyPairedSingleHits(hit1SecondBest, hit2, separationDist)
    
    if (pairScoreA < bestPairScore)
      bestPairScore = pairScoreA
    
    if (pairScoreB < bestPairScore)
      bestPairScore = pairScoreB
    
    bestPairScore
  }
  
  def verifyPairedHit(hit1: AlignResult, hit2: AlignResult, separationDist: Int): Int = {
    var bestPairScore = Int.MaxValue
    
    hit1 match {
      case RichSingleHit(pos1, isRC1, editDistance1) => {
        hit2 match {
          case RichSingleHit(pos2, isRC2, editDistance2) => {
            val pairScore = verifyPairedSingleHits(hit1.asInstanceOf[RichSingleHit], hit2.asInstanceOf[RichSingleHit], separationDist)
            
            if (pairScore < bestPairScore)
              bestPairScore = pairScore
          }
          case RichMultipleHits(bestPos2, bestScore2, isRC2, secondBestPos2, secondBestScore2) => {
            val pairScore = verifyPairedMixedHits(hit2.asInstanceOf[RichMultipleHits], hit1.asInstanceOf[RichSingleHit], separationDist)
            
            if (pairScore < bestPairScore)
              bestPairScore = pairScore
          }
          case _ =>
        }
      }
      case RichMultipleHits(bestPos1, bestScore1, isRC1, secondBestPos1, secondBestScore1) => {
        hit2 match {
          case RichSingleHit(pos2, isRC2, editDistance2) => {
            val pairScore = verifyPairedMixedHits(hit1.asInstanceOf[RichMultipleHits], hit2.asInstanceOf[RichSingleHit], separationDist)
            
            if (pairScore < bestPairScore)
              bestPairScore = pairScore
          }
          case RichMultipleHits(bestPos2, bestScore2, isRC2, secondBestPos2, secondBestScore2) => {
            val hit1Best = RichSingleHit(bestPos1, isRC1, bestScore1)
            val pairScoreA = verifyPairedMixedHits(hit2.asInstanceOf[RichMultipleHits], hit1Best.asInstanceOf[RichSingleHit], separationDist)
            
            val hit1SecondBest = RichSingleHit(secondBestPos1, isRC1, secondBestScore1)
            val pairScoreB = verifyPairedMixedHits(hit2.asInstanceOf[RichMultipleHits], hit1SecondBest.asInstanceOf[RichSingleHit], separationDist)
            
            if (pairScoreA < bestPairScore)
              bestPairScore = pairScoreA

            if (pairScoreB < bestPairScore)
              bestPairScore = pairScoreB
          }
          case _ =>
        }
      }
      case _ =>
    }
    
    bestPairScore
  }
  
  def parallelScoreBFA(readsWithBfaResults: IndexedSeq[(Read, AlignResult)], maxDist: Int) = {
    val (reads, results) = readsWithBfaResults.unzip
    scoreBFA(reads, results, maxDist, None)
  }

  def scoreBFA(reads: IndexedSeq[Read], results: IndexedSeq[AlignResult], maxDist: Int, snapResults: Option[Array[AlignResult]] /* used just for debugging */) = {
    var good = 0       // Total reads that are not nonsense (too many N's)
    var confident = 0  // Good reads for which we report a single hit
    var correct = 0    // Confident hits that we placed in the right position
    var multiHits = 0  // Good reads for which we report multiple hits
    var unmatched = 0  // Good reads which we didn't match
    
    // multiple hit cases
    var case1 = 0     // first choice is correct
    var case2 = 0     // 2nd choice is correct
    var case3 = 0     // neither is correct
    
    var case3Reads = List[(Read, AlignResult)]()
    
    val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r
    for (i <- 0 until reads.size) {
      if (i % 1000000 == 0 && i > 0)
        printf("Progress: %d/%d%n", i, reads.size) 
      val read = reads(i)
      if (read.data.count(_ == 'N') < read.data.length / 2) {
        good += 1
        results(i) match {
          case RichSingleHit(pos, rc, editDist) =>
            confident += 1
            // Get the true position of the read and compare it
            read.idStr match {
	          case ID_REGEX(piece, start, end) =>
	            val lowEnd = min(start.toLong, end.toLong)
	            val highEnd = max(start.toLong, end.toLong)
	            val (myPiece, myOffset) = GenomeLoader.genome.getLocation(pos)
	            if (myPiece == piece && myOffset >= lowEnd - maxDist && myOffset <= highEnd + maxDist)
	              correct += 1
	            else {
  	            println("For read " + i + ", wgsim says " + lowEnd + " - " + highEnd)
  	            println("BFA says:  " + results(i))	            
  	            snapResults match {
  	              case Some(results) => 
    	              println("SNAP says: " + results(i))
    	            case None => 
  	            }
  	            println
	            }	            
	          case _ =>
	            println("WARNING: Read ID " + read.idStr + " is not of the expected form! Are you using wgsim?")
	        }
          case RichMultipleHits(bestPos, bestScore, isRC, secondBestPos, secondBestScore) => {
            multiHits += 1

            read.idStr match {
	            case ID_REGEX(piece, start, end) => {
	              val lowEnd = min(start.toLong, end.toLong)
	              val highEnd = max(start.toLong, end.toLong)

  	            println("For read " + i + ", wgsim says " + lowEnd + " - " + highEnd)
  	            println("BFA says:  " + results(i))
  	            snapResults match {
  	              case Some(results) => 
    	              println("SNAP says: " + results(i))
    	            case None => 
  	            }
  	            println
  	            
  	            // characterize multiple hit
  	            val (myPiece, myOffset) = GenomeLoader.genome.getLocation(bestPos)
  	            if (myPiece == piece && myOffset >= lowEnd - maxDist && myOffset <= highEnd + maxDist)
  	              case1 += 1
                else {
    	            val (myPiece, myOffset) = GenomeLoader.genome.getLocation(secondBestPos)
    	            if (myPiece == piece && myOffset >= lowEnd - maxDist && myOffset <= highEnd + maxDist)
    	              case2 += 1
    	            else {
    	              case3 += 1
    	              case3Reads = (reads(i), results(i)) :: case3Reads
    	            }
                }
	            }
	          }
          }
          case NotFound =>
            unmatched += 1
          case _ =>
        }
      }
    }
    val errors = confident - correct
    val matched = good - unmatched
    printf("Good reads: %d/%d (%.2f%%)%n", good, reads.size, good * 100.0 / reads.size)
    printf("Total matched: %d/%d (%.2f%%)%n", matched, good, matched  * 100.0/ good)
    printf("Confident matches: %d/%d (%.2f%%)%n", confident, good, confident * 100.0 / good)
    printf("Alignment errors: %d/%d (%.2f%%)%n", errors, confident, errors * 100.0 / confident) 
    printf("Multiple hits: %d/%d (%.2f%%)%n", multiHits, reads.size, multiHits * 100.0 / reads.size) 
    printf("Multiple hits (case 1): %d/%d (%.2f%%)%n", case1, multiHits, case1 * 100.0 / multiHits)
    printf("Multiple hits (case 2): %d/%d (%.2f%%)%n", case2, multiHits, case2 * 100.0 / multiHits)
    printf("Multiple hits (case 3): %d/%d (%.2f%%)%n", case3, multiHits, case3 * 100.0 / multiHits)
  }

  def scorePairedBFA(reads1: IndexedSeq[Read], reads2: IndexedSeq[Read], results: IndexedSeq[AlignResult], maxDist: Int, snapResults1: Option[Array[AlignResult]], snapResults2: Option[Array[AlignResult]] /* used just for debugging */) = {
    if (reads1.length != reads2.length)
      throw new IllegalArgumentException("Read arrays 1 & 2 must be same length.")

    var good = 0       // Total reads that are not nonsense (too many N's)
    var confident = 0  // Good reads for which we report a single hit
    var correct = 0    // Confident hits that we placed in the right position
    var multiHits = 0  // Good reads for which we report multiple hits
    var unmatched = 0  // Good reads which we didn't match
    
    // multiple hit cases
    /*
    var case1 = 0     // first choice is correct
    var case2 = 0     // 2nd choice is correct
    var case3 = 0     // neither is correct
    
    var case3Reads = List[(Read, AlignResult)]()
    */
    
    val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r
    for (i <- 0 until reads1.size) {
      if (i % 1000000 == 0 && i > 0)
        printf("Progress: %d/%d%n", i, reads1.size) 
      val read1 = reads1(i)
      if (read1.data.count(_ == 'N') < read1.data.length / 2) {
        good += 1
        results(i) match {
          case RichPairSingleHit(pos1, pos2, score, isRC) => 
            confident += 1
            // Get the true position of the read and compare it
            // only need to look at id string of read1 b/c read2 will have same id string
            read1.idStr match {
	          case ID_REGEX(piece, start, end) =>
	            val lowEnd = min(start.toLong, end.toLong)
	            val highEnd = max(start.toLong, end.toLong)

	            val (myPiece1, myOffset1) = GenomeLoader.genome.getLocation(pos1)
	            val (myPiece2, myOffset2) = GenomeLoader.genome.getLocation(pos2)

	            if (
	              myPiece1 == piece && 
	              myPiece2 == piece &&
	              myOffset1 >= lowEnd - maxDist && myOffset1 <= highEnd + maxDist &&
	              myOffset2 >= lowEnd - maxDist && myOffset2 <= highEnd + maxDist
	            )
	              correct += 1
	            else {
  	            println("For read " + i + ", wgsim says " + lowEnd + " - " + highEnd)
  	            println("BFA says:  " + results(i))	            
  	            snapResults1 match {
  	              case Some(results) => 
    	              println("SNAP says: ")
    	              println("read1: " + results(i))
    	            case _ =>
  	            }
  	            snapResults2 match {
  	              case Some(results) =>
  	                println("read2: " + results(i))
    	            case _ =>
  	            }
  	            println
	            }	            
	          case _ =>
	            println("WARNING: Read ID " + read1.idStr + " is not of the expected form! Are you using wgsim?")
	        }
          case RichPairMultiHit(bestPos1, bestPos2, bestScore, secondBestPos1, secondBestPos2, secondBestScore, isRC) => {
            multiHits += 1

            read1.idStr match {
	            case ID_REGEX(piece, start, end) => {
	              val lowEnd = min(start.toLong, end.toLong)
	              val highEnd = max(start.toLong, end.toLong)

  	            println("For read " + i + ", wgsim says " + lowEnd + " - " + highEnd)
  	            println("BFA says:  " + results(i))
  	            snapResults1 match {
  	              case Some(results) => 
    	              println("SNAP says: ")
    	              println("read1: " + results(i))
    	            case _ => 
  	            }
  	            snapResults2 match {
  	              case Some(results) =>
  	                println("read2: " + results(i))
    	            case _ => 
  	            }
  	            println
  	            
  	            // characterize multiple hit
  	            /*
  	            val (myPiece, myOffset) = GenomeLoader.genome.getLocation(bestPos)
  	            if (myPiece == piece && myOffset >= lowEnd - maxDist && myOffset <= highEnd + maxDist)
  	              case1 += 1
                else {
    	            val (myPiece, myOffset) = GenomeLoader.genome.getLocation(secondBestPos)
    	            if (myPiece == piece && myOffset >= lowEnd - maxDist && myOffset <= highEnd + maxDist)
    	              case2 += 1
    	            else {
    	              case3 += 1
    	              case3Reads = (reads(i), results(i)) :: case3Reads
    	            }
                }
                */
	            }
	          }
          }
          case NotFound =>
            unmatched += 1
          case _ =>
        }
      }
    }
    val errors = confident - correct
    val matched = good - unmatched
    printf("Good reads: %d/%d (%.2f%%)%n", good, reads1.size, good * 100.0 / reads1.size)
    printf("Total matched: %d/%d (%.2f%%)%n", matched, good, matched  * 100.0/ good)
    printf("Confident matches: %d/%d (%.2f%%)%n", confident, good, confident * 100.0 / good)
    printf("Alignment errors: %d/%d (%.2f%%)%n", errors, confident, errors * 100.0 / confident) 
    printf("Multiple hits: %d/%d (%.2f%%)%n", multiHits, reads1.size, multiHits * 100.0 / reads1.size) 
    /*
    printf("Multiple hits (case 1): %d/%d (%.2f%%)%n", case1, multiHits, case1 * 100.0 / multiHits)
    printf("Multiple hits (case 2): %d/%d (%.2f%%)%n", case2, multiHits, case2 * 100.0 / multiHits)
    printf("Multiple hits (case 3): %d/%d (%.2f%%)%n", case3, multiHits, case3 * 100.0 / multiHits)
    */
  }

  // for finding similar regions
  def multiAlign(index: Index, currentPos: Long, readLen: Int, seedLen: Int, seedsToTry: Int, maxDist: Int, maxHits: Int): Set[Set[Long]] = {
    val read = new Array[Byte](readLen)
    GenomeLoader.genome.getSubstring(currentPos, currentPos + readLen, read)

    val matches = Set[Set[Long]]()
    
    val fwdHits = new LongArrayList
    val fwdSeen = new LongOpenHashSet
    val lv = new LandauVishkin(maxDist)
    val MAX_READ_LEN = 512 // Not a hard constraint; just for sizing the arrays
    val ref = new Array[Byte](MAX_READ_LEN)
        
    var seedNum = 0
    var offset = 0

    var fallBack = false
    // can count # of seeds with too many hits if want to watch for <some min #> of seeds violating maxHits

    // if any (or <some min #>) of the seeds violates maxHits, fall back to bfa (check each location)
    // -- gets around problem that index only includes some of the hits for popular seeds
    // else try candidates only
    while (seedNum < seedsToTry && !fallBack) {
      fwdHits.clear()
      if (!index.get(DNA.substringToLong(read, offset, offset + seedLen), fwdHits, maxHits)) {
        fallBack = true
      } else {
        var i = 0
        while (i < fwdHits.size) {
          val pos = fwdHits.getLong(i) - offset
          if (!fwdSeen.contains(pos)) {
            // Score the hit using Landau-Vishkin and remember it
            GenomeLoader.genome.getSubstring(pos, pos + read.length, ref)
            val score = lv.distance(ref, read.length, read, read.length, maxDist)
            fwdSeen.add(pos)

            // if within maxDist of read, record as a match
            if (score != -1)
              matches += Set(currentPos, pos)
          }
          i += 1
        }
      }      
      
      // Try next seed
      seedNum += 1
      offset = ((seedNum / (seedsToTry - 1.0)) * (read.length - seedLen)).round.toInt
    }
    
    if (fallBack)
      null
      // should actually be doing bfa (checking every position in the genome)
    else
      matches
  }
  
  // should go to fastq class
  def writeToFastq(readsWithBfaResults: Array[(Read, AlignResult)], fname: String) = {
    val fstream = new java.io.FileWriter(fname)
    val out = new java.io.BufferedWriter(fstream)
    
    readsWithBfaResults.foreach(r => 
      r match {
        case (read, res) => {
          res match {
            case RichSingleHit(pos, isRC, editDistance) => {
              out.write(asWgsimId(read.idStr, pos, read.data.length) + "\n")
              out.write(read.dataStr + "\n")
              out.write("+\n")
              if (read.quality != null)
                out.write(read.qualityStr + "\n")
              else
                out.write((1 to read.data.length).map(i => "2").mkString("") + "\n")
            }
            case _ => 
              // for now, only pass on reads that were unambiguously alignable
          }
        }
        case _ =>
      }
    )
    
    out.close
  }
  
  def writePairToFastq(readsWithBfaResults: Array[((Read, Read), AlignResult)], outputPrefix: String, confDiff: Int = 1) = {
    val fstream1 = new java.io.FileWriter(outputPrefix + "1.fq")
    val out1 = new java.io.BufferedWriter(fstream1)
    
    val fstream2 = new java.io.FileWriter(outputPrefix + "2.fq")
    val out2 = new java.io.BufferedWriter(fstream2)
    
    readsWithBfaResults.foreach(r => {
      val ((read1, read2), res) = r
      res match {
        case RichPairSingleHit(pos1, pos2, score, isRC) => {
          // make sure this is a valid hit (requires that pos1 & pos2 are on the same chromosome)
          if (asWgsimId(read1.idStr, pos1, pos2, read1.data.length) != null) {
            out1.write(asWgsimId(read1.idStr, pos1, pos2, read1.data.length) + "\n")
            out1.write(read1.dataStr + "\n")
            out1.write("+\n")
            if (read1.quality != null)
              out1.write(read1.qualityStr + "\n")
            else
              out1.write((1 to read1.data.length).map(i => "2").mkString("") + "\n")

            out2.write(asWgsimId(read2.idStr, pos1, pos2, read2.data.length) + "\n")
            out2.write(read2.dataStr + "\n")
            out2.write("+\n")
            if (read2.quality != null)
              out2.write(read2.qualityStr + "\n")
            else
              out2.write((1 to read2.data.length).map(i => "2").mkString("") + "\n")
          }
        }
        case RichPairMultiHit(bestPos1, bestPos2, bestScore, secondBestPos1, secondBestPos2, secondBestScore, isRC) => {
          //if (bestScore < secondBestScore) {
          if (bestScore + confDiff <= secondBestScore) {
            val pos1 = bestPos1
            val pos2 = bestPos2
            
            assert(asWgsimId(read1.idStr, pos1, pos2, read1.data.length) != null)
            
            out1.write(asWgsimId(read1.idStr, pos1, pos2, read1.data.length) + "\n")
            out1.write(read1.dataStr + "\n")
            out1.write("+\n")
            if (read1.quality != null)
              out1.write(read1.qualityStr + "\n")
            else
              out1.write((1 to read1.data.length).map(i => "2").mkString("") + "\n")

            out2.write(asWgsimId(read2.idStr, pos1, pos2, read2.data.length) + "\n")
            out2.write(read2.dataStr + "\n")
            out2.write("+\n")
            if (read2.quality != null)
              out2.write(read2.qualityStr + "\n")
            else
              out2.write((1 to read2.data.length).map(i => "2").mkString("") + "\n")
          }
        }
        case _ =>
          // for now, only pass on reads that were unambiguously alignable
      }
    })
    
    out1.close
    out2.close
  }
  
  def asWgsimId(id: String, pos: Long, readLen: Int): String = {
    val (myPiece, myOffset) = GenomeLoader.genome.getLocation(pos)
    List("@" + myPiece, myOffset, myOffset + readLen, 0 + ":" + id).mkString("_")
    // for now, since single-end, repeat pos rather than putting pos of other end
    // put 0 as stand in for rest of string (assuming that's for recording the changes from the reference made by wgsim)
    // we won't have that b/c this will be used for real data
  }

  def asWgsimId(id: String, pos1: Long, pos2: Long, readLen: Int): String = {
    val (myPiece1, myOffset1) = GenomeLoader.genome.getLocation(pos1)
    val (myPiece2, myOffset2) = GenomeLoader.genome.getLocation(pos2)
    
    val lowEnd = min(myOffset1, myOffset2)
    val highEnd = max(myOffset1, myOffset2)
    
    if (myPiece1 == myPiece2)
      List("@" + myPiece1, lowEnd + 1, highEnd + 1, 0 + ":" + id).mkString("_")
    else
      null
  }
}