package snap

import scala.collection.mutable.Map
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import scala.collection.JavaConversions._
import java.math.BigInteger
import java.security.MessageDigest

class GenomeHash(genome: Genome, readLen: Int, numMaps: Int = 100) {
  val map = Map[Long, List[Long]]()
  //val counts = new Long2IntOpenHashMap()
  //val counts = Map[Long, Int]()
  val countMaps = Array.fill(numMaps)(new Long2IntOpenHashMap())
  var popularHashCodes = Set[java.lang.Long]()
  
  def hashRead(read: String): Long = {
    import java.math.BigInteger
    import java.security.MessageDigest
    
    val digest = MessageDigest.getInstance("MD5")
    digest.update(read.getBytes)
    val md5sum = digest.digest  // matei said that md5 returns a byte array & that i could just take the first 8 bytes
    val bigInt = new BigInteger(1, md5sum)
    bigInt.longValue
  }
  
  def hashRead(read: String, randomColumns: List[Int]): Long = {
    // pull those bases out of the read string
    var bases = List[Char]()
    randomColumns.foreach(c => {
      bases = read.charAt(c) :: bases
    })
    bases = bases.reverse
    
    hashRead(bases.mkString)
  }

  def hashRead(read: Array[Byte], randomColumns: List[Int]): Long = {
    // pull those bases out of the read string
    var bases = List[Char]()
    var numColumns = randomColumns.length
    var c = 0
    var hash = 0L
    
    while (c < numColumns) {
      //bases = read(randomColumns(c)).toChar :: bases
      hash = (hash << 2) | DNA.BASE_TO_CODE(read(randomColumns(c)))
      c += 1
    }
    //bases = bases.reverse
    
    //hashRead(bases.mkString)
    hash
  }

  def getRandomColumns(numColumns: Int): List[Int] = {
    scala.util.Random.shuffle((0 until readLen).toList).slice(0, numColumns)
  }
  
  /*
  def hashGenome(readLen: Int) = {
    var (read, pos) = genome.getNextReadWithPos(readLen)
    while (read != "notFound") {
      map += ((hashRead(read), (read, pos)))
      //(read, pos) = genome.getNextReadWithPos(readLen)
    }
  }
  */

  def hashGenome = {
    genome.scanReads(readLen, {
      case pos => {
        val read = genome.substring(pos, pos + readLen)
        val h = hashRead(read)
        val l = map.get(h)
        if (l != None) {
          map.put(h, pos :: l.get)
        } else {
          map.put(h, List(pos))
        }
      }
    })
  }
  
  def hashGenome(numColumns: Int) = {
    val randomColumns = getRandomColumns(numColumns)
    genome.scanReads(readLen, {
      case pos => {
        val read = genome.substring(pos, pos + readLen)
        val h = hashRead(read, randomColumns)
        val l = map.get(h)
        if (l != None) {
          map.put(h, pos :: l.get)
        } else {
          map.put(h, List(pos))
        }
      }
    })
  }
  
  def updateCountMap(hash: Long) = {
    val whichMap = getMapNum(hash)
    val l = countMaps(whichMap).get(hash)
    try {
      if (l != countMaps(whichMap).defaultReturnValue) {
        countMaps(whichMap).add(hash, 1)
      } else {
        countMaps(whichMap).put(hash, 1)
      }
    } catch {
      case e => {
        println(hash)
        println(l)
        e.printStackTrace
      }
    }
  }
  
  def updateMap(hash: Long, pos: Long) = {
    val l = map.get(hash)
    try {
      if (l != None) {
        map.put(hash, pos :: l.get)
      } else {
        map.put(hash, List(pos))
      }
    } catch {
      case e => {
        println(hash)
        println("pos: " + pos)
        println(l)
        e.printStackTrace
      }
    }
  }
  
  def getMapNum(hash: Long): Int = {
    scala.math.abs((hash % (numMaps * 10)).toInt) / 10
  }
  
  def hashGenome(numColumns: Int, minBucketSize: Int) = {
    val randomColumns = getRandomColumns(numColumns)
    
    // first pass:  count how many times each hash code occurs
    println("counting hash code occurrences...")
    var start = System.currentTimeMillis
    
    var i = 0
    genome.scanReads(readLen, {
      case pos => {
        val read = new Array[Byte](readLen)
        genome.getSubstring(pos, pos + readLen, read)
        val h = hashRead(read, randomColumns)
        updateCountMap(h)
      }
      
      i += 1
    })
    println("Hashing took " + (System.currentTimeMillis - start)/1000.0 + "s")
    
    println("# buckets (unfiltered): " + countMaps.map(m => BigInt(m.keySet.size)).sum)
    
    // get popular hash codes
    println("Getting popular hash codes")
    
    start = System.currentTimeMillis
    i = 0
    while (i < numMaps) {
      println("Processing map " + i + "...")
      popularHashCodes ++= countMaps(i).keySet.toSet.filter(hash => {
        countMaps(i).get(hash) >= minBucketSize
      })
      i += 1
    }
    println("Getting popular hash codes took " + (System.currentTimeMillis - start)/1000.0 + "s")
    
    
    // second pass:  store positions for substrings that hash to popular hash codes
    println("storing positions that hash to popular hash codes...")
    start = System.currentTimeMillis
    genome.scanReads(readLen, {
      case pos => {
        val read = new Array[Byte](readLen)
        genome.getSubstring(pos, pos + readLen, read)
        val h = hashRead(read, randomColumns)
        if (popularHashCodes.contains(h)) {
          updateMap(h, pos)
        }
      }
    })
    println("2nd pass took " + (System.currentTimeMillis - start)/1000.0 + "s")
  }
  
  def getHashStats = {
    println("# Buckets: " + map.size)
    println("# Buckets with size < 10: " + map.values.count(_.size < 10))
    println("# Buckets with size 11-20: " + map.values.count(x => x.size >= 11 && x.size <= 20))
    println("# Buckets with size 21-30: " + map.values.count(x => x.size >= 21 && x.size <= 30))
    println("# Buckets with size 31-40: " + map.values.count(x => x.size >= 31 && x.size <= 40))
    println("# Buckets with size 41-50: " + map.values.count(x => x.size >= 41 && x.size <= 50))
    println("# Buckets with size 51-75: " + map.values.count(x => x.size >= 51 && x.size <= 75))
    println("# Buckets with size 76-100: " + map.values.count(x => x.size >= 76 && x.size <= 100))
    println("# Buckets with size 101-1000: " + map.values.count(x => x.size >= 101 && x.size <= 1000))
    println("# Buckets with size > 1000: " + map.values.count(_.size > 1000))
  }

  def getHashHistogram(maxBuckets: Int): Array[Int] = {
    var h = Array.fill(maxBuckets)(0)

    (0 until maxBuckets - 1).foreach(size => {
      h(size) = map.values.count(_.size == size)
      
      if (size % 10 == 0)
        println("Filling in histogram for bucket size of " + size + "/" + maxBuckets)
    })
    
    h(maxBuckets - 1) = map.values.count(_.size >= maxBuckets)
    h
  }
  
  def printBuckets(minSize: Int, maxSize: Int) = {
    val hashCodes = map.keySet.filter(hash => {
      val size = map.get(hash).get.size
      size >= minSize && size <= maxSize
    })
    
    hashCodes.foreach(hash => {
      println(hash)
      val l = map.get(hash)
      if (l != None) {
        l.get.foreach(pos => println(genome.substring(pos, pos + readLen) + " at " + pos))
      }
    })
  }
  
  def printBuckets(fname: String) = {
    val fstream = new java.io.FileWriter(fname)
    val out = new java.io.BufferedWriter(fstream)
    
    var h = 0
    val numHashCodes = popularHashCodes.size
    
    popularHashCodes.foreach(hash => {
      val l = map.get(hash)
      if (l != None) {
        out.write(hash + ": " + l.get.mkString(", ") + "\n")
      }
      
      if (h % 1000000 == 0)
        println("Printed " + h + "/" + numHashCodes + " buckets to file.")
      
      h += 1
    })
    
    out.close
  }
  
  def printNBuckets(minSize: Int, maxSize: Int, numBuckets: Int, editDistanceThreshold: Int, numNodeThreshold: Int) = {
    val fname = "/root/similarRegions.txt"
    val fstream = new java.io.FileWriter(fname)
    val out = new java.io.BufferedWriter(fstream)

    val hashCodes = map.keySet.filter(hash => {
      val size = map.get(hash).get.size
      size >= minSize && size <= maxSize
    }).toList
    
    val maxBucket = 
      if (numBuckets <= hashCodes.size) numBuckets
      else {
        println("Printing all buckets with size in given range")
        hashCodes.size
      }
    
    (1 to maxBucket).foreach(i => {
      if (i % 10 == 0)
        println("Considering bucket " + i + "/" + numBuckets)
        
      val hash = hashCodes(i)
      val l = map.get(hash)
      if (l != None) {
        val vertices = l.get.map(pos => (genome.substring(pos, pos + readLen), pos))
        val components = AdjacencyGraphHandler.findConnectedComponents(vertices, editDistanceThreshold)
        AdjacencyGraphHandler.printConnectedComponentsSummary(vertices, components, numNodeThreshold, Some(out))
      }
    })

    out.close
  }
}