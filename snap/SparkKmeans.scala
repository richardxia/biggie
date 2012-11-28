package snap

import snap._

object SparkKmeans {
  def main(args: Array[String]) : Unit = {
    args match {
      case Array("kmeansRaw") => {
        val s1 = "GGAACAGAATGGAGTCATCGAATGAAATCGAATGGAATCATCATCAAATGAAATCAAATGGAGTCATCGTATGGACTCCAACGGAATCATCATCGACTGG"
        val s2 = "ATTCCATTCCTTTCTTTTGACAGGGTATCATTGTGTCACTGAGGCTGGAGTACAGTGGCACAATCTCAGCTCACATTGCATGTCAACTTTCCATTTCATT"
        
        println("dist: " + (new Levenshtein).distance(s1, s2))  // 59
        
        // create two clusters
        var points = Array.fill(20)(ReadPoint(""))
        (0 until 10).foreach(i => {
          points(i) = ReadPoint(perturbByN(s1, 10))
        })
        
        (10 until 20).foreach(i => {
          points(i) = ReadPoint(perturbByN(s2, 10))
        })
        
        kmeans(points, 2)
      }
    }
  }

  def perturbByN(str: String, n: Int): String = {
    val strArray = str.toCharArray

    val randomSubset = scala.util.Random.shuffle(0 until str.length).take(n)
    randomSubset.foreach(println)
    
    // flip n random bases in the string
    (0 until n).foreach(i => {
      val idx = randomSubset(i)
      println(i)
      println(idx)
      strArray(idx) = DNA.COMPLEMENT(strArray(idx)).toChar
    })
    
    new String(strArray)
  }
  
  def kmeans(points: Seq[Point], k: Int) = {
    // initialize the centroids to be random points
    var centroids = Array.fill(k)(points(scala.util.Random.nextInt(k)))
    var membership = Array.fill(points.length)(0)

    // iterate
    (0 until 10).foreach(i => {
      // E step:  for each point, find the closest centroid
      membership = updateMembership(points, centroids)

      // M step:  for each cluster, find the new centroid
      centroids = updateCentroids(points, membership, k)
    })
    
    // print out summary of clusters
    (0 until k).foreach(c => {
      println("Cluster " + c + ":")
      println("  Consensus string: " + centroids(c).toString)
      println("  Membership: " + membership.count(_ == c))
    })
  }
  
  def updateMembership(points: Seq[Point], centroids: Seq[Point]): Array[Int] = {
    points.toArray.map(p => {
      (0 until centroids.length).map(c => {
        (c, centroids(c).distance(p))
      })
      .sortWith(_._2 < _._2)
      .head
      ._1
    })
  }
  
  def updateCentroids(points: Seq[Point], membership: Seq[Int], k: Int): Array[Point] = {
    (0 until k).map(c => {
      val clusterPoints = (0 until points.length).flatMap(j => {
        // if point's membership is in c, include it
        if (membership(j) == c)
          List(points(j))
        else
          Nil
      })
      new ReadPoint(getConsensusStr(clusterPoints))
    }).toArray
  }

  def getConsensusStr(points: IndexedSeq[Point]) = {
    assert(points(0).isInstanceOf[ReadPoint])

    val readLen = points(0).asInstanceOf[ReadPoint].read.length
    val consensusStr = Array.fill(readLen)(' ')
    
    val ct = Array.fill(readLen)(Array.fill(4)(0))
    (0 until readLen).foreach(r => {
      // get most popular char
      (0 until points.length).foreach(h => {
        ct(r)(DNA.BASE_TO_CODE(points(h).asInstanceOf[ReadPoint].read.charAt(r))) += 1
      })
    })
    
    // figure out which char is most popular at each base
    (0 until readLen).foreach(r => {
      consensusStr(r) = DNA.CODE_TO_BASE(argMax(ct(r)))
    })

    new String(consensusStr)
  }
  
  def argMax(baseCts: Array[Int]): Int = {
    val max = baseCts.max
    val argMax = baseCts.indexOf(max)
    argMax
  }
}
