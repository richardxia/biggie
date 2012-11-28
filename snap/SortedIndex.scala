package snap

import it.unimi.dsi.fastutil.longs.LongList
import java.util.Arrays
import scala.math.max
import scala.math.min
import it.unimi.dsi.fastutil.longs.LongArrayList
import it.unimi.dsi.fastutil.longs.Long2ObjectMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.ints.IntArrayList

class SortedIndex(seedLen: Int, shardBits: Int, arrays: Array[Array[Long]], overflow: Long2ObjectMap[IntArrayList])
extends Index(seedLen) {
  val numShards = 1 << shardBits
  val shardMask = numShards - 1
  val keyShift = 64 - seedLen * 2 + shardBits
  val keyMask = (-1L) << keyShift
  val invKeyMask = ~keyMask
  val MULTI_HIT = invKeyMask
  
  override def get(seed: Long, dest: LongList, maxResults: Int): Boolean = {
    val shard = (seed & shardMask).toInt
    val key = seed >> shardBits
    val shiftedKey = key << keyShift
    val array = arrays(shard)
    var pos = Arrays.binarySearch(array, shiftedKey)
    if (pos < 0) {
      // binarySearch returns -insertionPoint-1 if the value is not in the index at which
      // it would be placed; this is exactly where we want to start searching from
      pos = -(pos + 1)
    }
    if (pos < array.length && (array(pos) & keyMask) == shiftedKey && (array(pos) & invKeyMask) == MULTI_HIT) {
      var list = overflow.get(seed)
      if (list.size() > maxResults)
        return false
      var i = 0
      while (i < min(list.size(), maxResults)) {
        dest.add(list.get(i) & 0xffffffffL)
        i += 1
      }
      return true
    }
    val maxPos = min(array.length, min(pos.toLong + maxResults, Int.MaxValue).toInt)
    while (pos < maxPos && (array(pos) & keyMask) == shiftedKey) {
      dest.add(array(pos) & invKeyMask)
      pos += 1
    }
    if (pos+1 < maxResults && (array(pos+1) & keyMask) == shiftedKey) {
      // Too many results -- clear the list and return false
      dest.clear()
      return false
    } else {
      return true
    }
  }
}

class SortedIndexBuilder(seedLen: Int, estimatedSize: Long = 1024, shardBits: Int = 10) extends IndexBuilder(seedLen) {
  val numShards = 1 << shardBits
  val shardMask = numShards - 1
  val keyShift = 64 - seedLen * 2 + shardBits
  val keyMask = (-1L) << keyShift
  val invKeyMask = ~keyMask
  val MULTI_HIT = invKeyMask
  
  val MAX_REPEATS = 300

  val overflow = new Long2ObjectOpenHashMap[IntArrayList]
  
  val arrays = Array.tabulate(numShards)(_ => new GrowingLongArray((estimatedSize / numShards).toInt))
  
  override def add(seed: Long, pos: Long): Unit = {
    val shard = (seed & shardMask).toInt
    val key = seed >> shardBits
    val entry = (key << keyShift) | pos
    arrays(shard).add(entry)
  }
  
  override def build(): Index = {
    val list = new LongArrayList
    val plainArrays = arrays.zipWithIndex.map { case (a, shard) =>
      Arrays.sort(a.data, 0, a.size)
      list.clear()
      var i = 0
      while (i < a.size) {
        var key = a.data(i) >>> keyShift
        var j = i + 1
        while (j < a.size && a.data(j) >>> keyShift == key)
          j += 1
        var repCount = j - i
        if (repCount < 3) {
          var k = i
          while (k < j) {
            list.add(a.data(k))
            k += 1
          }
        } else {
          list.add((key << keyShift) | MULTI_HIT)
          val overList = new IntArrayList
          var k = i
          while (k < j) {
            overList.add((a.data(k) & invKeyMask).toInt)
            k += 1
          }
          overflow.put((key << shardBits) | shard, overList)
        }
        i = j
      }
      a.data = null // Let it get GCed
      list.toLongArray()
    }
    println("Total entries = " + plainArrays.map(_.length.toLong).sum)
    println("Longest array = " + plainArrays.map(_.length).max)
    new SortedIndex(seedLen, shardBits, plainArrays, overflow)
  }
}

class GrowingLongArray(initialCapacity: Int = 0) {
  var data = new Array[Long](max(initialCapacity, 2))
  var size = 0
  
  def add(n: Long) {
    if (size >= data.length) {
      val newData = new Array[Long]((data.length * 1.6 + 1).toInt)
      System.arraycopy(data, 0, newData, 0, data.length)
      data = newData
    }
    data(size) = n
    size += 1
  }
}
