package snap

import it.unimi.dsi.fastutil.longs.LongList
import scala.collection.mutable.ListBuffer
import scala.collection.JavaConversions._
import scala.util.control.Breaks.{break, breakable}
import scala.collection.mutable.{Map, Set}

// TODO:  make (seed, pos) tuple into a case class

class BucketIndex(seedLen: Int, prefixLen: Int, estimatedSize: Long) extends Index(seedLen) {
  var tupleArray = new Array[(String, Long)](estimatedSize.toInt)
  var lastIdx = 0
  
  var bucketIndex = scala.collection.immutable.Map[String, List[(String, Long)]]()
  
  var prefixByNumOccurrences = Map[Int, List[String]]()
  
  // unimplemented for now; don't plan to use this for looking up
  override def get(seed: Long, dest: LongList, maxResults: Int): Boolean = {
    true
  }
  
  // unimplemented b/c we add seeds as Strings (too long to fit in a Long)
  def add(seed: Long, pos: Long): Unit = {
  }
  
  def add(seed: String, pos: Long): Unit = {
    tupleArray(lastIdx) = (seed, pos)
    lastIdx += 1
  }
  
  def groupByPrefix: Unit = {
    val listBuf = new ListBuffer[(String, Long)]
    tupleArray.foreach(tuple => 
      if (tuple != null)
        listBuf += tuple
    )

  	bucketIndex = listBuf.toList.groupBy(tuple => tuple._1.substring(0, prefixLen))

  	// make reverse mapping from # occurrences to prefix
  	bucketIndex.keySet.foreach(prefix => {
	  val hits = bucketIndex.get(prefix).get
	  val numOccurrences = hits.size
	  val prevContents = prefixByNumOccurrences.get(numOccurrences)

	  if (prevContents == None)
	    prefixByNumOccurrences.put(numOccurrences, List(prefix))
	  else
	    prefixByNumOccurrences.put(numOccurrences, prefix :: prevContents.get)
  	})
  }
  
  def getBucket(targetNumHits: Int): Unit = {
    val (prefix, numHits) = getPrefixWithTargetNumHits(targetNumHits)
    
    println("Prefix " + prefix + " occurs " + numHits + " times:")

    val hits = bucketIndex.get(prefix).get
    hits.foreach(hit => {
      println(hit._1 + " at " + hit._2)
    })
  }
  
  def getPrefixWithTargetNumHits(targetNumHits: Int, hitDiff: Int = 5): (String, Int) = {
    var prefix = ""
    var numHits = 0
    
    breakable {
      bucketIndex.keySet.toList.foreach(p => {
        val h = bucketIndex.get(p).get.size
        if (h >= (targetNumHits - hitDiff) && h <= (targetNumHits + hitDiff)) {
          prefix = p
          numHits = h
          break()
        }
      })
    }
    
    (prefix, numHits)
  }
  
  def printClusterSummary(prefix: String, editDistanceThreshold: Int): Unit = {
    val clusterMap = AdjacencyGraphHandler.getNeighborhoods(bucketIndex.get(prefix).get, editDistanceThreshold)
    
    // print a summary of the cluster map
    clusterMap.keySet.foreach(t => 
      println(t + " has " + clusterMap.get(t).get.size + " neighbors.")
    )
  }

  def getBucketSummaries(numOccurrences: Int, numNeighborsThreshold: Int, editDistanceThreshold: Int): Unit = {
    // For each prefix occurring numOccurrences times,
    // Get its bucket
    // Check each bucket member against the others in the bucket (does it have >= numNeighborsThreshold?)
    // Print bucket summary if at least one of its members had enough neighbors
    
    if (prefixByNumOccurrences.get(numOccurrences) == None) {
      println("No prefixes have the desired numOccurrences.")
      return
    }
    
    val prefixes = prefixByNumOccurrences.get(numOccurrences).get
    
    var prefixNum = 0
    val numPrefixes = prefixes.size
    
    prefixes.foreach(prefix => {
      val clusterMap = AdjacencyGraphHandler.getNeighborhoods(bucketIndex.get(prefix).get, editDistanceThreshold)
      
      // Do any of the seeds have more than numNeighborsThreshold?
      val seedsWithNeighborhood = clusterMap.keySet.filter(seed => clusterMap.get(seed).get.size > numNeighborsThreshold)
      
      if (seedsWithNeighborhood.size > 0) {
        println("prefix " + prefix + " has " + seedsWithNeighborhood.size + " seeds with at least " + numNeighborsThreshold 
          + " neighbors (within distance of " + editDistanceThreshold + ").")
      }
      
      if (prefixNum % 1000 == 0)
        println("processed " + prefixNum + "/" + numPrefixes + " prefixes.")

      prefixNum += 1
    })
  }
  
  def getSingleBucketSummary(prefix: String, editDistanceThreshold: Int): Unit = {
    // For each seed, print its neighbors (along with the edit distance between them)
    
    if (bucketIndex.get(prefix) == None) {
      println("Prefix not found.")
      return
    }
    
    println("prefix: " + prefix)
    val l = new Levenshtein
        
    val neighborhoods = AdjacencyGraphHandler.getNeighborhoods(bucketIndex.get(prefix).get, editDistanceThreshold)
    neighborhoods.keySet.foreach(seed => {
      val neighbors = neighborhoods.get(seed).get
      if (neighbors.size > 1) {
        println("neighbors of " + seed._1 + "(" + seed._2 + "):")
        neighbors.foreach(neighbor => {
          println("  " + neighbor._1 + "(" + neighbor._2 + ")" + " with distance of " + l.distance(seed._1, neighbor._1))
        })
      }
    })
  }
  
  def validatePrefix(prefix: String): Boolean = {
    val found = bucketIndex.get(prefix) != None
    if (!found)
      println("Prefix not found.")
      
    return found
  }
  
  def findConnectedComponentsAndPrintSummary(prefix: String, editDistanceThreshold: Int, numNodeThreshold: Int): Unit = {
    if (!validatePrefix(prefix)) return
    
    println("prefix: " + prefix)
    
    // Find connected components
    val vertices = bucketIndex.get(prefix).get
    val components = AdjacencyGraphHandler.findConnectedComponents(vertices, editDistanceThreshold)
    
    // Print summary of connected components (only those with at least numNodeThreshold nodes)
    AdjacencyGraphHandler.printConnectedComponentsSummary(vertices, components, numNodeThreshold, Some(new java.io.PrintWriter(System.out)))
  }
}

class BucketIndexBuilder(seedLen: Int, prefixLen: Int, estimatedSize: Long = 1024) extends IndexBuilder(seedLen) {
  val index = new BucketIndex(seedLen, prefixLen, estimatedSize)
  
  override def add(seed: Long, pos: Long): Unit = index.add(seed, pos)
  
  def add(seed: String, pos: Long): Unit = index.add(seed, pos)
  
  override def build(): Index = {/*index.printStats(); */index}
}
