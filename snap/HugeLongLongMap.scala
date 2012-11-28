package snap

import scala.collection.Traversable
import it.unimi.dsi.fastutil.HashCommon

/**
 * Like a LongLongMap, but bigger (can hold more than 2 billion items).
 * Because Java arrays are limited in size to 2 billion entries, sometimes
 * LongLongMap is not enough, so this class holds multiple LongLongMaps
 * and hashes each key to one of them.
 */
class HugeLongLongMap(emptyValue: Long = 0L, estimatedSize: Long = 16, loadFactor: Double = 0.7)
extends Traversable[(Long, Long)] {
  val SHIFT = 6
  val SHARDS = 1 << SHIFT
  val MASK = SHARDS - 1
  val estShardSize = (1.05 * estimatedSize / SHARDS).toInt
  val maps = Array.tabulate(SHARDS)(_ => new LongLongMap(emptyValue, estShardSize, loadFactor))
  
  private def mapId(key: Long): Int = {
    key.toInt & MASK
  }
  
  def apply(key: Long): Long = {
    maps(key.toInt & MASK)(key >> SHIFT)
  }
  
  def update(key: Long, value: Long): Unit = {
    maps(key.toInt & MASK)(key >> SHIFT) = value
  }
  
  def getOrUpdate(key: Long, value: Long): Long = {
    maps(key.toInt & MASK).getOrUpdate(key >> SHIFT, value)
  }
  
  override def foreach[T](f: ((Long, Long)) => T) {
    maps.foreach(_.foreach(f))
  }
}
