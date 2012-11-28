package snap

import scala.math.min

/**
 * Performs the variant of Landau-Vishkin k-difference used in CloudBurst.
 * This takes O(k^2 * min(m, n)) time to find the best edit distance not exceeding
 * k between a pattern of length m and a prefix of a text of length n.
 */
class LandauVishkin(maxK: Int) {
  val L: Array[Array[Int]] = Array.fill(2, 2*maxK +1)(0)  // (numErrs % 2)(offset)

  /**
   * Returns the lowest edit distance between pattern and a prefix of text, assuming that
   * there exists a prefix for which the edit distance is at most k, or -1 otherwise.
   */
  def distance(text: String, pattern: String, k: Int): Int = {
    if (k > maxK)
      throw new IllegalArgumentException("k=" + k + " is bigger than maxK=" + maxK)
    
    val m = pattern.length
    val n = text.length
    
    if (m == 0)
      return 0

    var e = 0
    while (e <= k) {
      var d = -e
      while (d <= e) {
        // Pick the "best" (farthest) index in pattern with e-1 errors to continue from
        var best = -1
        if (e > 0) {
          if (d > -e && d < e) { // Try adding mismatch error
            val up = L((e-1)&1)(k+d) + 1
            if (up > best)
              best = up
          }
          if (d > -(e-1)) { // Try adding a shift-left error
            val left = L((e-1)&1)(k+d-1)
            if (left > best)
              best = left
          }
          if (d < e-1) { // Try adding a shift-right error
            val right = L((e-1)&1)(k+d+1) + 1
            if (right > best)
              best = right
          }
        } else {
          best = 0
        }

        var end = min(m, n-d)
        while (best < end && pattern.charAt(best) == text.charAt(best+d))
          best += 1
        L(e&1)(k+d) = best

        if (best == m) // Reached end of pattern, so this is the best value of e
          return e

        d += 1
      }
      e += 1
    }
    return -1 // No alignment found with <= k errors
  }

  /**
   * Returns the lowest edit distance between pattern and a prefix of text, assuming that
   * there exists a prefix for which the edit distance is at most k, or -1 otherwise.
   */
  def distance(text: Array[Byte], textLen: Int, pattern: Array[Byte], patternLen: Int, k: Int): Int = {
    if (k > maxK)
      throw new IllegalArgumentException("k=" + k + " is bigger than maxK=" + maxK)
    
    if (patternLen == 0)
      return 0

    var e = 0
    while (e <= k) {
      var d = -e
      while (d <= e) {
        // Pick the "best" (farthest) index in pattern with e-1 errors to continue from
        var best = -1
        if (e > 0) {
          if (d > -e && d < e) { // Try adding mismatch error
            val up = L((e-1)&1)(k+d) + 1
            if (up > best)
              best = up
          }
          if (d > -(e-1)) { // Try adding a shift-left error
            val left = L((e-1)&1)(k+d-1)
            if (left > best)
              best = left
          }
          if (d < e-1) { // Try adding a shift-right error
            val right = L((e-1)&1)(k+d+1) + 1
            if (right > best)
              best = right
          }
        } else {
          best = 0
        }

        var end = min(patternLen, textLen-d)
        while (best < end && pattern(best) == text(best+d))
          best += 1
        L(e&1)(k+d) = best

        if (best == patternLen) // Reached end of pattern, so this is the best value of e
          return e

        d += 1
      }
      e += 1
    }
    return -1 // No alignment found with <= k errors
  }
}
