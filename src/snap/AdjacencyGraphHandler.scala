package snap

import scala.collection.mutable.{Set, Map}

object AdjacencyGraphHandler {
  def getNeighborhoods(vertices: List[(String, Long)], editDistanceThreshold: Int): Map[(String, Long), List[(String, Long)]] = {
    val clusterMap = Map[(String, Long), List[(String, Long)]]()
    val l = new Levenshtein
    
    vertices.foreach(vertex => {
      // find each vertex's neighbors
      clusterMap.put(vertex, vertices.filter(h => {
        l.distance(vertex._1, h._1) <= editDistanceThreshold
      }))
    })
    
    clusterMap
  }
  def constructGraph(vertices: List[(String, Long)], editDistanceThreshold: Int): Set[Set[Int]] = {
    // get edges
    val edges = Set[Set[Int]]()
    
    val neighborhoods = getNeighborhoods(vertices, editDistanceThreshold)
    neighborhoods.keySet.foreach(seed => {
      val neighbors = neighborhoods.get(seed).get
      neighbors.foreach(neighbor => {
        if (seed != neighbor)
          edges += Set(vertices.indexOf(seed), vertices.indexOf(neighbor))
      })
    })
    
    edges
  }

  def findConnectedComponents(vertices: List[(String, Long)], editDistanceThreshold: Int): Set[Set[Int]] = {
    val components = Set[Set[Int]]()
    val edges = constructGraph(vertices, editDistanceThreshold)
    
    // initialize components
    vertices.foreach(vertex => components += Set(vertices.indexOf(vertex)))
    
    edges.foreach(edge => {
      val endpoints = edge.toList
      val v1 = endpoints(0)
      val v2 = endpoints(1)

      var set1:Set[Int] = null
      var set2:Set[Int] = null
      
      // get component in which each vertex occurs
      components.foreach(component => {
        if (component.contains(v1)) set1 = component
        if (component.contains(v2)) set2 = component
      })
      
      // if in same component, do nothing
      // otherwise, merge those components
      if (set1 != set2) {
        components -= set1
        components -= set2
        components += set1 ++ set2
      }
    })
    
    components
  }
  
  // Print summary of connected components (only those with at least numNodeThreshold nodes)
  def printConnectedComponentsSummary(vertices: List[(String, Long)], components: Set[Set[Int]], numNodeThreshold: Int, out: Option[java.io.Writer], editDistanceThreshold: Int = 10) = {
    /*
    val toFile = 
    out match {
      case Some(writer) => true
      case None => false
    }
    */
    
    // want to print consensus str & avg # diffs from consensus

    var componentNum = 0
    components.foreach(component => {
      if (component.size >= numNodeThreshold) {
        componentNum += 1
        /*
        if (toFile) {
          out.get.write("component " + componentNum + " has " + component.size + " nodes:\n")
          out.get.write("Avg # diffs from consensus: " + getAvgNumDiffs(vertices, component, editDistanceThreshold) + "\n")
          out.get.write("Consensus:\n")
          out.get.write("  " + getConsensusStr(vertices, component) + "\n")
          out.get.write("Members:\n")
        }
        else {
        */
          println("component " + componentNum + " has " + component.size + " nodes:")
          println("Avg # diffs from consensus: " + getAvgNumDiffs(vertices, component, editDistanceThreshold))
          println("Consensus:")
          println("  " + getConsensusStr(vertices, component))
          println("Members:")
        //}
        component.foreach(nodeNum => {
          /*
          if (toFile)
            out.get.write("  " + vertices(nodeNum)._1 + " at " + vertices(nodeNum)._2 + "\n")
          else
          */
            println("  " + vertices(nodeNum)._1 + " at " + vertices(nodeNum)._2)
        })
        /*
        if (toFile)
          out.get.write("\n")
        else
        */
          println
      }
    })
  }

  def getConsensusStr(vertices: List[(String, Long)], component: Set[Int]): String = {
    val (verticesStr, verticesLong) = vertices.unzip
    val readLen = verticesStr.head.length
    
    // get component as str list
    val componentAsStrList = getComponentAsStrList(vertices, component)
    
    // get consensus string
    val consensusStr = Array.fill(readLen)(' ')
    
    val ct = Array.fill(readLen)(Array.fill(4)(0))
    (0 until readLen).foreach(r => {
      // get most popular char
      (0 until componentAsStrList.length).foreach(h => {
        ct(r)(DNA.BASE_TO_CODE(componentAsStrList(h).charAt(r))) += 1
      })
    })
    
    // figure out which char is most popular at each base
    (0 until readLen).foreach(r => {
      consensusStr(r) = DNA.CODE_TO_BASE(argMax(ct(r)))
    })

    new String(consensusStr)
  }
  
  def getComponentAsStrList(vertices: List[(String, Long)], component: Set[Int]): List[String] = {
    val (verticesStr, verticesLong) = vertices.unzip
    
    var componentAsStrList = List[String]()
    
    (0 until verticesStr.length).foreach(s => {
      if (component.contains(s))
        componentAsStrList = verticesStr(s) :: componentAsStrList
    })
    
    componentAsStrList
  }
  
  def argMax(baseCts: Array[Int]): Int = {
    val max = baseCts.max
    val argMax = baseCts.indexOf(max)
    argMax
  }
  
  def getAvgNumDiffs(vertices: List[(String, Long)], component: Set[Int], editDistanceThreshold: Int): Double = {
    val consensusStr = getConsensusStr(vertices, component)
    
    val hitStrList = getComponentAsStrList(vertices, component)
    
    val numDiffs = hitStrList.map(hit => getNumDiffs(hit, consensusStr, editDistanceThreshold))
    
    (numDiffs.sum.toDouble / numDiffs.length).round
  }
  
  def getNumDiffs(str1: String, str2: String, editDistanceThreshold: Int): Int = {
    assert(str1.length == str2.length)
    
    var numDiffs = 0
    
    /*
    (0 until str1.length).foreach(i => {
      if (str1.charAt(i) != str2.charAt(i)) numDiffs += 1
    })
    */
    val lv = new LandauVishkin(str1.length)
    numDiffs = lv.distance(str1, str2, str1.length)
    
    numDiffs
  }

  
}