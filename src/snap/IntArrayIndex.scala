package snap

import it.unimi.dsi.fastutil.longs.LongList
import java.util.Arrays
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.ints.IntArrayList

class IntArrayIndexBuilder(seedLen: Int, estimatedSize: Long) extends IndexBuilder(seedLen) {
  val shardShift = math.max(seedLen*2 - 32, 1) // Always used in order to turn seed to 32-bit
  val shift = math.max(computeShift(estimatedSize), shardShift)
  val mask = (1 << shift) - 1
  val array = new Array[Array[Int]](mask + 1)
  val overflow = new Long2ObjectOpenHashMap[IntArrayList]
  
  val MAX_REPEATS = 300
  
  override def add(seed: Long, pos: Long): Unit = {
    val shard = seed & ((1 << shardShift) - 1)
    val shiftedKey = (seed >> shardShift).toInt
    val h = hash(shiftedKey)
    val bucket = (((h << shardShift) | shard) & mask).toInt
    //printf("For key %d: s=%d, sk=%d, h=%d, buck=%d%n", seed, shard, shiftedKey, h, bucket)
    if (array(bucket) == null) {
      array(bucket) = new Array[Int](2)
      array(bucket)(0) = shiftedKey
      array(bucket)(1) = pos.toInt
    } else {
      val a = array(bucket)
      var i = 0
      while (i < a.length / 2) {
        if (a(2*i) == shiftedKey) {
          if (a(2*i+1) == -1) {
            //println("seed: " + seed)
            //println("overflow: " + overflow)
            //println("array: " + array.deep.mkString(","))
            val list = overflow.get(seed)
            if (list.size() < MAX_REPEATS)
              list.add(pos.toInt)
          } else {
            val list = new IntArrayList(4)
            list.add(a(2*i+1))
            list.add(pos.toInt)
            overflow.put(seed, list)
            a(2*i+1) = -1
          }
          return
        }
        i += 1
      }
      // Not found, so expand the array
      array(bucket) = Arrays.copyOf(a, a.length + 2)
      array(bucket)(a.length) = shiftedKey
      array(bucket)(a.length+1) = pos.toInt
    }
  }
  
  override def build(): Index = new Index(seedLen) {
    override def get(seed: Long, dest: LongList, maxEntries: Int): Boolean = {
      val shard = seed & ((1 << shardShift) - 1)
      val shiftedKey = (seed >> shardShift).toInt
      val h = hash(shiftedKey)
      val bucket = (((h << shardShift) | shard) & mask).toInt
      val a = array(bucket)
      if (a != null) {
        var i = 0
        while (i < a.length / 2) {
          if (a(2*i) == shiftedKey) {
            if (a(2*i+1) == -1) {
              val list = overflow.get(seed)
              val toAdd = math.min(maxEntries + 1, list.size())
              var j = 0
              while (j < toAdd) {
                dest.add(list.getInt(j))
                j += 1
              }
            } else {
              dest.add(a(2*i+1))
            }
            if (dest.size() > maxEntries) {
              // Found too many entries
              dest.clear()
              return false
            } else {
              return true
            }
          }
          i += 1
        }
      }
      return true
    }
  }
  
  def computeShift(estimatedSize: Long): Int = {
    var shift = 1
    while (shift < 30 && (1 << shift) < estimatedSize / 8)
      shift += 1
    shift
  }

  // MurmurHash3 from http://sites.google.com/site/murmurhash
  private final def hash(key: Int): Int = {
    var x = key
    x ^= x >>> 16;
    x *= 0x85ebca6b;
    x ^= x >>> 13;
    x *= 0xc2b2ae35;
    x ^= x >>> 16;
    return x;
  }
}
