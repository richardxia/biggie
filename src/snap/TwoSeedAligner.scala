package snap
import java.util.Random
import it.unimi.dsi.fastutil.longs.LongArrayList
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import scala.math.abs

/**
 * Align reads based on nearby matches of a two seeds.
 *
 * The algorithm parameters are as follows:
 * - seedsToTry: how many different seeds to test (should be at least 2).
 * - maxDist: maximum edit distance from reference genome to consider an alignment valid.
 * - confDiff: minimum difference in scores (edit distances from the reference genome)
 *   between the closest and second-closest alignments to return a confident result
 *   (as opposed to a "multiple hits" result).
 * - maxHits: caps the number of hits a seed can have where we will still try using all
 *   of them for alignment; seeds with too many hits are costly to test and are likely
 *   hitting repetitive sequences where we can't get a confident alignment anyway.
 * - lookForBest: if this is false, stop at the first alignment found with less than
 *   maxDist edit distance from the reference genome rather than finding the best
 *   alignment; this will return MultipleHits in this case because we can't know
 *   whether there are other hits that we haven't scored.
 */
class TwoSeedAligner(genome: Genome, index: Index, seedsToTry: Int, maxDist: Int, 
  confDiff: Int = 3, maxHits: Int = 40, lookForBest: Boolean = true)
extends Aligner {
  private val seedLen = index.seedLen
  private val rand = new Random(42)
  private val lv = new LandauVishkin(maxDist + confDiff - 1)
  private val fwdHits = new LongArrayList
  private val rcHits = new LongArrayList 
  private val fwdPos = new LongArrayList
  private val rcPos = new LongArrayList
  private val fwdSeen = new LongOpenHashSet
  private val rcSeen = new LongOpenHashSet
  private val toInsert = new LongArrayList

  // Byte arrays for the read's reverse complement and any genome substrings we compare with
  private val MAX_READ_LEN = 512 // Not a hard constraint; just for sizing the arrays
  private val rc = new Array[Byte](MAX_READ_LEN)
  private val ref = new Array[Byte](MAX_READ_LEN)
  
  def align(read: Array[Byte]): AlignResult = {
    if (read.length < seedLen)
      throw new IllegalArgumentException("Read too short (< seedLen)")
    if (read.length > MAX_READ_LEN)
      throw new IllegalArgumentException("Read too long (> " + MAX_READ_LEN + ")")
    if (Utils.count(read, 'N') > maxDist)
      return NotFound
    
    DNA.getReverseComplement(read, rc)
    
    fwdPos.clear()
    rcPos.clear()
    fwdSeen.clear()
    rcSeen.clear()
    if (rand.nextInt() % 128 == 0) {
      // A bit of a hack in case the arrays get too long on some hit
      fwdPos.clear()
      rcPos.clear()
      fwdSeen.trim()
      rcSeen.trim()
    }
    var bestScore = Int.MaxValue
    var bestPos = 0L
    var bestIsRC = false
    var secondBestScore = Int.MaxValue
    
    var offset = 0
    var foundUntestedRepeat = false
    
    var seedNum = 0
    while (seedNum < seedsToTry) {
      fwdHits.clear()
      if  (!index.get(DNA.substringToLong(read, offset, offset + seedLen), fwdHits, maxHits)) {
        foundUntestedRepeat = true
      } else {
        toInsert.clear()
        var i = 0
        while (i < fwdHits.size) {
          val pos = fwdHits.getLong(i) - offset
          if (!fwdSeen.contains(pos)) {
            // Search for another seed near it
            // TODO: take advantage of sorted lists to make this not be quadratic
            var j = 0
            var found = false
            while (j < fwdPos.size && !found) {
              if (abs(fwdPos.getLong(j) - pos) <= maxDist + confDiff - 1)
                found = true
              j += 1
            }
            if (found) {
              // Score the hit using Landau-Vishkin and remember it
              genome.getSubstring(pos, pos + read.length, ref)
              val distToCheck = {
                if (bestScore > maxDist) maxDist + confDiff - 1
                else if (secondBestScore < bestScore + confDiff) bestScore - 1 // We're ambiguous, so need to find a better unambiguous hit
                else bestScore + confDiff - 1 // No need to search for worse scores than this
              }
              val score = lv.distance(ref, read.length, read, read.length, distToCheck)
              fwdSeen.add(pos)
              //printf("Forward score for pos %d (from seed %s) was %d%n",
              //       pos, read.substring(offset, offset + seedLen), score)
              if (score != -1 && score <= maxDist + confDiff - 1) {
                if (!lookForBest)
                  return MultipleHits // We found that it aligns somewhere, so we're done
                // TODO: deal with nearby hits specially due to indels -- see SingleSeedAligner
                if (score < bestScore) {
                  secondBestScore = bestScore
                  bestPos = pos
                  bestScore = score
                  bestIsRC = false
                } else if (score < secondBestScore) {
                  secondBestScore = score
                }
                // Return early if our best scores so far are too close for an unambigous hit
                if (bestScore <= maxDist && (secondBestScore < confDiff || bestScore < confDiff && secondBestScore < bestScore + confDiff))
                  return MultipleHits
              }
            }
            toInsert.add(pos)
          }
          // Done searching for matches; also insert our new positions into fwdPos
          fwdPos.addAll(toInsert)
          i += 1
        }
      }
      
      rcHits.clear()
      if (!index.get(DNA.rcSubstringToLong(read, offset, offset + seedLen), rcHits, maxHits)) {
        foundUntestedRepeat = true
      } else {
        toInsert.clear()
        var i = 0
        while (i < rcHits.size) {
          val pos = rcHits.getLong(i) - offset
          if (!rcSeen.contains(pos)) {
            // Search for another seed near it
            // TODO: take advantage of sorted lists to make this not be quadratic
            var j = 0
            var found = false
            while (j < rcPos.size && !found) {
              if (abs(rcPos.getLong(j) - pos) <= maxDist + confDiff - 1)
                found = true
              j += 1
            }
            if (found) {
              // Score the hit using Landau-Vishkin and remember it
              genome.getSubstring(pos, pos + read.length, ref)
              val distToCheck = {
                if (bestScore > maxDist) maxDist + confDiff - 1
                else if (secondBestScore < bestScore + confDiff) bestScore - 1 // We're ambiguous, so need to find a better unambiguous hit
                else bestScore + confDiff - 1 // No need to search for worse scores than this
              }
              val score = lv.distance(ref, read.length, rc, read.length, distToCheck)
              rcSeen.add(pos)
              //printf("RC score for pos %d (from seed %s) was %d%n",
              //       pos, read.substring(offset, offset + seedLen), score)
              if (score != -1 && score <= maxDist + confDiff - 1) {
                if (!lookForBest)
                  return MultipleHits // We found that it aligns somewhere, so we're done
                // TODO: deal with nearby hits specially due to indels -- see SingleSeedAligner
                if (score < bestScore) {
                  secondBestScore = bestScore
                  bestPos = pos
                  bestScore = score
                  bestIsRC = false
                } else if (score < secondBestScore) {
                  secondBestScore = score
                }
                // Return early if our best scores so far are too close for an unambigous hit
                if (bestScore <= maxDist && (secondBestScore < confDiff || bestScore < confDiff && secondBestScore < bestScore + confDiff))
                  return MultipleHits
              }
            }
            toInsert.add(pos)
          }
          // Done searching for matches; also insert our new positions into rcPos
          rcPos.addAll(toInsert)
          i += 1
        }
      }
      
      // Try next seed
      seedNum += 1
      offset = ((seedNum / (seedsToTry - 1.0)) * (read.length - seedLen)).round.toInt
    }
    
    if (bestScore <= maxDist) {
      if (secondBestScore >= bestScore + confDiff)
        return SingleHit(bestPos, bestIsRC)
      else
        return MultipleHits
    } else if (!lookForBest && foundUntestedRepeat && fwdSeen.isEmpty && rcSeen.isEmpty) {
      // Kind of a conundrum: we found only untested repeats, but we didn't want to try
      // all of them for local alignment, so what do we return? Let's go with Ambiguous
      // for now, even though there's a small chance the read doesn't align anywhere.
      return MultipleHits
    } else {
      return NotFound
    }
  }
}
