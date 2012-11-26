package snap

import it.unimi.dsi.fastutil.longs.LongList
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongArrayList
import it.unimi.dsi.fastutil.longs.LongLists
import scala.math.min

/**
 * A simple implementation of Index that just has a hashtable from seeds to
 * lists of positions. This works well when the number of seeds is small (for
 * example, when the seed length is small) but is memory-inefficient for
 * larger numbers of seeds with fewer hits per seed.  
 */
class SimpleIndex(seedLen: Int) extends Index(seedLen) {
  val map = new Long2ObjectOpenHashMap[LongList]

  def add(seed: Long, pos: Long): Unit = {
    val list = map.get(seed)
    if (list != null) {
      list.add(pos)
    } else {
      val newList = new LongArrayList
      newList.add(pos)
      map.put(seed, newList)
    }
  }

  override def get(seed: Long, dest: LongList, max: Int): Boolean = {
    val list = map.get(seed)
    if (list == null || list.size() > max) {
      return false
    } else {
      dest.addAll(list.subList(0, min(list.size(), max)))
      return true
    }
  }
}

class SimpleIndexBuilder(seedLen: Int) extends IndexBuilder(seedLen) {
  val index = new SimpleIndex(seedLen)
  
  override def add(seed: Long, pos: Long): Unit = index.add(seed, pos)
  
  override def build() = index
}
