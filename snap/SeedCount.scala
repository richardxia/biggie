package snap

import scala.collection.JavaConversions._
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashMap
import scala.math.{abs, max, min}
import scala.io.Source
import scala.util.Random

import it.unimi.dsi.fastutil.ints._
import it.unimi.dsi.fastutil.longs._

class TestRead(val lowEnd: Int, val highEnd: Int, val data: String) {}

object SeedCount {
  def main(args: Array[String]) {
    args match {
      case Array(refFile, testFile, seedLen) =>
        run(refFile, testFile, seedLen.toInt)
      case _ =>
        System.err.println("Usage: SeedCount <refFile> <testFile> <seedLen>")
        System.exit(1)
    }
  }

  def run(refFile: String, testFile: String, seedLen: Int) {
    val index = buildIndex(refFile, seedLen)
  }

  // TODO: could also index the reverse complement
  def buildIndex(refFile: String, seedLen: Int): Long2ObjectMap[IntList] = {
    val lines = Source.fromFile(refFile).getLines
    val data = lines.filter(!_.startsWith(">")).map(_.toUpperCase).mkString
    val index = new Long2ObjectOpenHashMap[IntList](data.length)

    var lastN = -1   // last index where we saw an 'N'
    var key = 0L
    val mask = (1L << (seedLen * 2)) - 1
    for (i <- 0 until data.size - seedLen) {
      val base = data.charAt(i)
      key = ((key << 2) | DNA.BASE_TO_CODE(data.charAt(i))) & mask
      if (base == 'N')
        lastN = i
      if (lastN <= i - seedLen) {
        var list = index.get(key)
        if (list == null) {
          list = new IntArrayList(1)
          index.put(key, list)
        }
        list.add(i - seedLen + 1)
      }
      if (i % 1000000 == 0)
        println("Position: %d/%d".format(i, data.size))
    }
    return index
  }

  def align1(read: String, index: Long2ObjectMap[IntList], seedLen: Int, maxAttempts: Int): Int = {
    if (read.startsWith("NNN") && read.endsWith("NNN"))
      return -1 // There are lots of NNNN reads in simulated datasets
    val rc = DNA.reverseComplement(read)
    var offset = 0
    var attempt = 0
    var rand = new Random(42)
    while (attempt < maxAttempts) {
      val pos = alignWithSeed(read, index, offset, seedLen, 1)
      val pos2 = alignWithSeed(rc, index, offset, seedLen, 1)
      if (pos >= 0 && pos2 < 0)
        return pos - offset
      else if (pos2 >= 0 && pos < 0)
        return pos2 - offset
      // Either no hit, or too many hits; try another seed
      attempt += 1
      offset += seedLen
      while (offset > read.length - seedLen)
        offset = rand.nextInt(seedLen)
    }
    return -1
  }

  val fwdPos = new IntArrayList(128)
  val rcPos = new IntArrayList(128)
  val candidates = new IntArrayList(128)

  def align2(read: String, index: Long2ObjectMap[IntList], seedLen: Int, maxAttempts: Int, maxHits: Int = 10): Int = {
    if (read.startsWith("NNN") && read.endsWith("NNN"))
      return -1 // There are lots of NNNN reads in simulated datasets
    fwdPos.clear
    rcPos.clear
    var attempt = 0
    var offset = 0
    var rand = new Random(42)
    while (attempt < maxAttempts) {
      val fwdHits = index.get(DNA.substringToLong(read, offset, offset + seedLen))
      val rcHits = index.get(DNA.rcSubstringToLong(read, offset, offset + seedLen))
      // TODO: make this not be quadratic in maxHits by using sorted lists
      if (fwdHits != null && fwdHits.size() < maxHits) {
        // Check if any of our new hits are close to an old one
        var i = 0
        while (i < fwdPos.size()) {
          var j = 0
          while (j < fwdHits.size()) {
            if (abs(fwdPos.get(i) - (fwdHits.get(j) - offset)) < read.length / 4)
              return fwdPos.get(i)
            j += 1
          }
          i += 1
        }
        // Add in our new hits
        var j = 0
        while (j < fwdHits.size()) {
          fwdPos.add(fwdHits.get(j) - offset)
          j += 1
        }
      }
      if (rcHits != null && rcHits.size() < maxHits) {
        // Check if any of our new hits are close to an old one
        var i = 0
        while (i < rcPos.size()) {
          var j = 0
          while (j < rcHits.size()) {
            if (abs(rcPos.get(i) - (rcHits.get(j) - offset)) < read.length / 4)
              return rcPos.get(i)
            j += 1
          }
          i += 1
        }
        // Add in our new hits
        var j = 0
        while (j < rcHits.size()) {
          rcPos.add(rcHits.get(j) - offset)
          j += 1
        }
      }
      attempt += 1
      offset += seedLen
      while (offset > read.length - seedLen)
        offset = rand.nextInt(seedLen)
    }
    return -1
  }

  def align2b(read: String, index: Long2ObjectMap[IntList], rawData: String, maxDiff: Int, seedLen: Int, maxAttempts: Int, maxHits: Int = 10, secondDiff: Int = 2, maxGap: Int = 18): Int = {
    if (read.startsWith("NNN") && read.endsWith("NNN"))
      return -1 // There are lots of NNNN reads in simulated datasets
    fwdPos.clear
    rcPos.clear
    candidates.clear
    var attempt = 0
    var offset = 0
    var rand = new Random(42)
    while (attempt < maxAttempts) {
      val fwdHits = index.get(DNA.substringToLong(read, offset, offset + seedLen))
      val rcHits = index.get(DNA.rcSubstringToLong(read, offset, offset + seedLen))
      // TODO: make this not be quadratic in maxHits by using sorted lists
      if (fwdHits != null && fwdHits.size() < maxHits) {
        // Check if any of our new hits are close to an old one
        var i = 0
        while (i < fwdPos.size()) {
          var j = 0
          while (j < fwdHits.size()) {
            if (abs(fwdPos.get(i) - (fwdHits.get(j) - offset)) < maxGap)
              candidates.add(fwdPos.get(i))
            j += 1
          }
          i += 1
        }
        // Add in our new hits
        var j = 0
        while (j < fwdHits.size()) {
          fwdPos.add(fwdHits.get(j) - offset)
          j += 1
        }
      }
      if (rcHits != null && rcHits.size() < maxHits) {
        // Check if any of our new hits are close to an old one
        var i = 0
        while (i < rcPos.size()) {
          var j = 0
          while (j < rcHits.size()) {
            if (abs(rcPos.get(i) - (rcHits.get(j) - offset)) < read.length / 4)
              candidates.add(rcPos.get(i))
            j += 1
          }
          i += 1
        }
        // Add in our new hits
        var j = 0
        while (j < rcHits.size()) {
          rcPos.add(rcHits.get(j) - offset)
          j += 1
        }
      }
      if (candidates.size > 0) {
        var best = candidates.get(0)
        var bestScore = maxDiff + 10
        var secondBestScore = maxDiff + 10
        for (i <- 0 until candidates.size) {
          val score = checkAlignment(read, rawData, candidates.get(i), maxDiff)
          if (score != -1 && score < bestScore) {
            secondBestScore = bestScore
            bestScore = score
            best = candidates.get(i)
          } else if (score == bestScore) {
            secondBestScore = bestScore
          }
        }
        if (bestScore <= maxDiff && secondBestScore >= bestScore + secondDiff)
          return best
        else if (bestScore <= maxDiff)
          return -1 // got two or more good hits
        else
          candidates.clear
      }
      attempt += 1
      offset += seedLen
      if (offset > read.length - seedLen) {
        offset = rand.nextInt(seedLen)
        //fwdPos.clear
        //rcPos.clear
      }
    }
    return -1
  }

  def alignWithSeed(read: String, index: Long2ObjectMap[IntList], offset: Int, seedLen: Int, repeatCutoff: Int): Int = {
    // Ignore the string if there is an N
    val idx = read.indexOf('N', offset)
    if (idx != -1 && idx < offset + seedLen)
      return -1
    val list = index.get(DNA.substringToLong(read, offset, offset + seedLen))
    if (list == null) {
      //println("Seed did not hit")
      return -1
    } else if (list.size() > repeatCutoff) {
      //println("Seed hit a repeat with too many occurrences")
      return -1
      //return list.getInt(0)
    } else {
      //println("Seed hit a string with " + list.size() + " occurrences")
      return list.getInt(0)
    }
  }

  val lv = new LandauVishkin(16)

  def checkAlignment(read: String, text: String, pos: Int, k: Int): Int = {
    val substr = text.substring(max(0, pos), min(text.length, pos + read.length))
    val fwdScore = lv.distance(substr, read, k)
    if (fwdScore != -1)
      return fwdScore
    else
      return lv.distance(substr, DNA.reverseComplement(read), k)
  }

  def readReads(file: String, prefix: String): Array[TestRead] = {
    def toRead(t: Seq[String]): TestRead = {
      val end1 = t(0).split('_')(1).toInt
      val end2 = t(0).split('_')(2).toInt
      val lowEnd = min(end1, end2)
      val highEnd = max(end1, end2)
      val data = t(1)
      new TestRead(lowEnd, highEnd, data)
    }
    Source.fromFile(file).getLines.sliding(4, 4)
      .filter(x => x(0).startsWith(prefix) && Character.isDigit(x(0).charAt(prefix.length)))
      .map(toRead).toArray
  }

  def validate(rawData: String, reads: Array[TestRead], maxDiff: Int, func: String => Int) {
    def score(r: TestRead): Int = {
      val p = func(r.data);
      if (p == -1) 0
      else if (maxDiff != -1 && SeedCount.checkAlignment(r.data, rawData, p, maxDiff) == -1) 0
      else if(p >= r.lowEnd-15 && p <= r.highEnd+15) 2
      else 1
    }
    val scores = reads.map(score)
    val total = scores.size.toDouble
    println("Matched: " + scores.count(_ != 0) / total)
    println("Correct: " + scores.count(_ == 2) / total)
    println("Errors:  " + scores.count(_ == 1) / scores.count(_ != 0).toDouble)
  }
}
