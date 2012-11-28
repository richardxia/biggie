package snap

import it.unimi.dsi.fastutil.longs.LongList
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongArrayList
import it.unimi.dsi.fastutil.longs.LongLists
import scala.math.min

class LongSeedIndex(seedLen: Int, estimatedSize: Long) extends Index(seedLen) {
  private val NO_HIT = -1
  private val MULTI_HIT = -2
  
  // Contains a position for single-hit seeds, or MULTI_HIT for multiple
  private val map = new HugeLongLongMap(emptyValue = NO_HIT, estimatedSize = this.estimatedSize)
  
  // Contains only lists of positions for multi-hit seeds
  private val overflow = new Long2ObjectOpenHashMap[LongList]
  
  def add(seed: Long, pos: Long): Unit = {
    val old = map.getOrUpdate(seed, pos)
    if (old != NO_HIT) {
      if (old == MULTI_HIT) {
        overflow.get(seed).add(pos)
      } else {
        val list = new LongArrayList(2)
        list.add(old)
        list.add(pos)
        overflow.put(seed, list)
        map(seed) = MULTI_HIT
      }
    }
  }

  override def get(seed: Long, dest: LongList, max: Int): Boolean = {
    val p = map(seed)
    if (p == NO_HIT || max == 0) {
      return true
    } else if (p == MULTI_HIT) {
      val list = overflow.get(seed)
      if (list.size() > max)
        return false
      var i = 0
      while (i < list.size()) {
        dest.add(list.get(i))
        i += 1
      }
      return true
    } else {
      dest.add(p)
    }
  }
  
  def printStats() {
    import scala.collection.JavaConversions._
    println("Total entries: " + map.maps.map(_.size).sum)
    println("Overflow entries: " + overflow.size())
    val values = overflow.values
    println("Overflow of 2: " + values.count(_.size() == 2))
    println("Overflow of 3: " + values.count(_.size() == 3))
    println("Overflow of 4-10: " + values.count(x => x.size() >= 4 && x.size() <= 10))
    println("Overflow of 11-100: " + values.count(x => x.size() >= 11 && x.size() <= 100))
    println("Overflow of 101-1000: " + values.count(x => x.size() >= 101 && x.size() <= 1000))
    println("Overflow of > 1000: " + values.count(x => x.size() >= 1001))
  }
}

class LongSeedIndexBuilder(seedLen: Int, estimatedSize: Long = 1024) extends IndexBuilder(seedLen) {
  val index = new LongSeedIndex(seedLen, estimatedSize)
  
  override def add(seed: Long, pos: Long): Unit = index.add(seed, pos)
  
  override def build(): Index = {index.printStats(); index}
}
