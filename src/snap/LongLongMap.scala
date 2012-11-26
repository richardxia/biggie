package snap

import java.util.Arrays

/** 
 * A simple insert-only hashtable with keys and values as longs. Requires an "empty value"
 * that won't be used by any real entries to be set in order to identify missing entries.
 */
class LongLongMap(emptyValue: Long = 0L, estimatedSize: Int = 16, loadFactor: Double = 0.7)
extends Traversable[(Long, Long)]
{
  private var elements = 0     // number of elements currently in table
  private var capacity = nextPowerOf2(estimatedSize / loadFactor)
  private var mask = capacity - 1
  private var data = new Array[Long](2*capacity)
  private var growThreshold = (capacity * loadFactor).toInt
  Arrays.fill(data, emptyValue)
  
  def update(key: Long, value: Long): Unit = {
    if (value == emptyValue)
      throw new IllegalArgumentException("update() called with empty value")
    var pos = hash(key) & mask
    var i = 1
    while (data(2*pos) != key && data(2*pos+1) != emptyValue) {
      pos = (pos + i) & mask
      i += 1
    }
    if (data(2*pos) != key || data(2*pos+1) == emptyValue) {
      elements += 1
      if (elements > growThreshold)
        growTable()
    }
    data(2*pos) = key
    data(2*pos+1) = value
  }
  
  def apply(key: Long): Long = {
    var pos = hash(key) & mask
    var i = 1
    while (data(2*pos+1) != emptyValue) {
      if (data(2*pos) == key)
        return data(2*pos+1)
      pos = (pos + i) & mask
      i += 1
    }
    return emptyValue
  }
  
  def getOrUpdate(key: Long, value: Long): Long = {
    if (value == emptyValue)
      throw new IllegalArgumentException("getOrUpdate() called with empty value")
    var pos = hash(key) & mask
    var i = 1
    while (data(2*pos) != key && data(2*pos+1) != emptyValue) {
      pos = (pos + i) & mask
      i += 1
    }
    if (data(2*pos) == key) {
      return data(2*pos+1)
    } else {
      data(2*pos) = key
      data(2*pos+1) = value
      elements += 1
      if (elements > growThreshold)
        growTable()
      return emptyValue
    }
  }
  
  // MurmurHash3 from http://sites.google.com/site/murmurhash
  private final def hash(key: Long): Int = {
    var x = key
    x ^= x >>> 33
    x *= 0xff51afd7ed558ccdL
    x ^= x >>> 33
    x *= 0xc4ceb9fe1a85ec53L
    x ^= x >>> 33
    x.toInt
  }
  
  // Double the capacity and re-hash everything into a new table
  private def growTable() {
    var oldData = data
    var oldCapacity = capacity
    capacity *= 2
    mask = capacity - 1
    data = new Array[Long](2 * capacity)
    Arrays.fill(data, emptyValue)
    elements = 0
    for (i <- 0 until oldCapacity)
      if (oldData(2*i + 1) != emptyValue)
        update(oldData(2*i), oldData(2*i + 1))
    growThreshold = (capacity * loadFactor).toInt
  }
  
  private def nextPowerOf2(x: Double): Int = {
    var i = 1
    while (i < x && i != 0)
      i *= 2
    if (i == 0) // We overflowed and there's no bigger power of 2 than x!
      throw new IllegalArgumentException("Argument too big: " + x)
    return i
  }
  
  override def size = elements
  
  override def foreach[T](func: ((Long, Long)) => T) {
    for (i <- 0 until capacity)
      if (data(2*i+1) != emptyValue)
        func((data(2*i), data(2*i+1)))
  }
}
