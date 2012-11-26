package snap

import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import it.unimi.dsi.fastutil.longs.LongList
import scala.math.min

// Index for mapping from seed to # of occurrences in the genome
// Necessary b/c current indexes only allow keeping track of up to 100 hits per seed
class LongIntIndex(seedLen: Int, estimatedSize: Long) extends Index(seedLen) {
  private val NO_HIT = -1
  private val MULTI_HIT = -2
  
  // Contains counts for how many times each seed occurred (treat single- & multi-hit cases the same way)
  private val map = new Long2IntOpenHashMap
  map.defaultReturnValue(NO_HIT)
  
  def add(seed: Long): Unit = {
    val old = map.get(seed)
    if (old != NO_HIT) {
      map.put(seed, old + 1)
    } else {
      map.put(seed, 1)
    }
  }
  
  // must implement since I'm inheriting from Index, but don't plan to use
  override def get(seed: Long, dest: LongList, max: Int): Boolean = {
    true
  }
  
  def getNumOccurrences(seed: Long): Int = {
    map.get(seed)
  }
  
  def printStats() = {
    import scala.collection.JavaConversions._
    println("Total entries: " + map.keySet.size())
    val values = map.values
    println("2 Occurrences: " + values.count(_ == 2))
    println("3 Occurrences: " + values.count(_ == 3))
    println("4-10 Occurrences: " + values.count(x => x >= 4 && x <= 10))
    println("11-100 Occurrences: " + values.count(x => x >= 11 && x <= 100))
    println("101-1000 Occurrences: " + values.count(x => x >= 101 && x <= 1000))
    println("> 1000 Occurrences: " + values.count(_ > 1000))
  }
  
  def getPopularSeeds(numHits: Int) {
    import scala.collection.JavaConversions._

    val popularSeeds = map.keySet.toList.filter(seed => map.get(seed) >= numHits)
    
    popularSeeds.foreach(seed => {
      println(DNA.longToString(seed, seedLen) + " occurs " + map.get(seed) + " times.")
    })
  }
}

class LongIntIndexBuilder(seedLen: Int, estimatedSize: Long = 1024) extends IndexBuilder(seedLen) {
  val index = new LongIntIndex(seedLen, estimatedSize)
  
  override def add(seed: Long, pos: Long): Unit = index.add(seed)
  
  override def build(): Index = {index.printStats(); index}
}
