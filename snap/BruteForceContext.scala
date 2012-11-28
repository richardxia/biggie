package snap

import spark.SparkContext
import SparkContext._

class PairedBruteForceContext (master: String, sparkMaster: String, accessKey: String, secretKey: String, fastqFile1: String, fastqFile2: String, 
  snapFile1: String, snapFile2: String, outputDest: String, outputBucket: String, params: PairedBruteForceParams) {
  // Setup

  // load raw reads
  println("Loading raw reads...")
  val reads1 = FASTQ.read(fastqFile1, true /* keep quality */)
  val reads2 = FASTQ.read(fastqFile2, true /* keep quality */)

  // load snap results
  println("Loading SNAP results...")
  var fileIn = new java.io.FileInputStream(snapFile1)
  var in = new java.io.ObjectInputStream(fileIn)
  val snapResults1 = in.readObject.asInstanceOf[Array[AlignResult]]
  val readsWithSnapResults1 = reads1.zip(snapResults1)

  fileIn = new java.io.FileInputStream(snapFile2)
  in = new java.io.ObjectInputStream(fileIn)
  val snapResults2 = in.readObject.asInstanceOf[Array[AlignResult]]
  val readsWithSnapResults2 = reads2.zip(snapResults2)

  val readsWithSnapResults = readsWithSnapResults1.zip(readsWithSnapResults2)

  assert(reads1.length == reads2.length && reads1.length == snapResults1.length && reads1.length == snapResults2.length)

  // set output destination
  var outputPath = ""
  outputDest match {
    case "hdfs" => outputPath = "hdfs://" + master + ":9000/" + outputBucket
    case "s3" => outputPath = "s3n://" + accessKey + ":" + secretKey + "@hashingGenome/" + outputBucket
    case _ => {
      println("Incorrect output destination.  Should be hdfs or s3.")
      System.exit(0)
    }
  }
  
  // set up spark context
  val sparkDest = 
    if (sparkMaster == "local") "local[8]"
    else "1@" + master + ":5050"

  //def getParams: PairedBruteForceParams = params
}