package snap

import scala.math.min

/**
 * Computes Levenshtein distance, given two strings that are aligned.
 * Implementation from http://www.merriampark.com/ld.htm
 */
class Levenshtein {
  def minOfThree(a: Int, b: Int, c: Int): Int = {
    min(min(a, b), c)
  }
  
  def distance(s: String, t: String): Int = {
    val n = s.length
    val m = t.length
    
    if (n == 0) return m
    if (m == 0) return n
    
    val d = Array.fill(n + 1, m + 1)(0)
    
    (1 to n).foreach(i => d(i)(0) = i)
    (1 to m).foreach(j => d(0)(j) = j)
    
    var cost = 0
    
    (1 to n).foreach(i => {
      val s_i = s(i-1)
      
      (1 to m).foreach(j => {
        val t_j = t(j-1)
        
        if (s_i == t_j) cost = 0
        else cost = 1
        
        d(i)(j) = minOfThree(d(i-1)(j) + 1, d(i)(j-1) + 1, d(i-1)(j-1) + cost)
      })
    })
    
    d(n)(m)
  }
}