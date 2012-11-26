package snap

import spark.SparkContext
import SparkContext._

class Duplicates (master: String, accessKey: String, secretKey: String, folder: String) {
  val sc = new SparkContext("1@" + master + ":5050", "analyzeDuplicates")
  val duplicates = sc.objectFile[(String, Array[Long])]("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + folder)

  def findDuplicates(readLen: Int, minBucketSize: Int): Unit = {
    val genomeSize = GenomeLoader.genome.totalSize
    
    sc.parallelize(0 until (genomeSize/2).toInt).flatMap(pos => {
      var l = List[(String, Long)]()

      val pos1 = 2 * pos.toLong
      val read1 = GenomeLoader.genome.substring(pos1, pos1 + readLen)
      if (!read1.contains('N'))
        l = (read1, pos1) :: l

      val pos2 = 2 * pos.toLong + 1
      val read2 = GenomeLoader.genome.substring(pos2, pos2 + readLen)
      if (!read2.contains('N'))
        l = (read2, pos2) :: l
        
      l
    })
    .groupByKey
    .filter(group => group._2.size > minBucketSize)
    .map(group => (group._1, group._2.toArray))
    .saveAsObjectFile("s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + folder)
    
    System.exit(0)
  }
  
  def sizeHistogram(fname: String) = {
    val sizeHist = duplicates.map(i => (i._2.size, 1)).reduceByKey(_ + _).collect

    val fstream = new java.io.FileWriter(fname)
    val out = new java.io.BufferedWriter(fstream)

    sizeHist.map(i => out.write(i._1 + ", " + i._2 + "\n"))

    out.close
    System.exit(0)
  }
  
  def posList(fname: String) = {
    val fstream = new java.io.FileWriter(fname)
    val out = new java.io.BufferedWriter(fstream)

    // print out all positions
    duplicates.collect.map(group => {
      val posArray = group._2
      posArray.foreach(pos => out.write(pos + "\n"))
    })
    
    out.close
    System.exit(0)
  }

  def popularDuplicates(fname: String, minRepeats: Int, sampleFraction: Double) = {
    val fstream = new java.io.FileWriter(fname)
    val out = new java.io.BufferedWriter(fstream)
    
    // look at popular strings
    out.write("popular:\n")
    val popularStrings = duplicates.filter(group => group._2.size >= minRepeats).sample(false, sampleFraction, System.currentTimeMillis.toInt).collect
    popularStrings.map(i => out.write(i._1 + ", " + i._2.size + "\n"))
    out.write("\n")
    
    out.close
    System.exit(0)
  }
}