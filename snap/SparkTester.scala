package snap

import spark.SparkContext
import SparkContext._

object SparkTester {
  def main(args: Array[String]): Unit = {
    val spark = new SparkContext("local", "wc")
    val file = spark.textFile("/root/test.txt")
    val counts = file.flatMap(line => line.split(" ")).
      map(word => (word, 1)).
      reduceByKey(_ + _)
    counts.saveAsTextFile("/root/test")
  }
}