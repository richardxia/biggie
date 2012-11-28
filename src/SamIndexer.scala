package biggie

import scala.io.Source

class SamIndexer(fileName: String) {
  val file = Source.fromFile(fileName)
  val positions: Array[(Int, Long)] = file.getLines().map(line => {
    val values = line.split('\t')
    val coord: Int = values(0).toInt
    val offset: Long = values(1).toLong
    (coord, offset)
  }).toArray
  //val posMap = new collection.mutable.HashMap[Int, Long]
  //for (line <- file.getLines()) {
  //  val values = line.split('\t')
  //  val coord = values(0).toInt
  //  val offset = values(1).toLong
  //  posMap(coord) = offset
  //}

  def apply(coord: Int): Long = {
    //posMap(coord)
    val index = binarySearch(coord, 0, positions.size)
    positions(index)._2
  }

  private def binarySearch(coord: Int, min: Int, max: Int): Int = {
    if (max <= min) {
      return max
    }
    val mid = (min + max)/2
    if (positions(mid)._1 > coord) {
      return binarySearch(coord, min, mid-1)
    } else if (positions(mid)._1 < coord) {
      return binarySearch(coord, mid, max)
    } else {
      return mid
    }
  }
}

object SamIndexer {
  def main(args: Array[String]) {
    if (args.size != 1) {
      println("Usage: SamIndexer file.sam")
      println("Creates a file named file.sam.sai")
      println("Each line is 'genome-coordinate<TAB>byte-offset'")
      System.exit(1)
    }

    val samFile = new java.io.RandomAccessFile(args(0), "r")
    val indexFile = new java.io.FileWriter(args(0) + ".sai")
    //val indexFile = new java.io.FileOutputStream(args(0) + ".sai")

    try {
      var lastCoord = 0
      //val buf = java.nio.ByteBuffer.allocate(16)
      while (true) {
        val pos = samFile.getFilePointer()
        val line = samFile.readLine()

        if (line == null) {
          indexFile.close()
          samFile.close()
          return
        }

        if (!line.startsWith("@")) {
          val coord = SAM.parseEntry(line).position
          if (coord != lastCoord) {
            indexFile.write(coord.toString + "\t" + pos.toString + "\n")
            //buf.putLong(0, coord)
            //buf.putLong(8, pos)
            //indexFile.write(buf.array())
            lastCoord = coord
          }
        }
      }
    } catch {
      case e: java.io.IOException => e.printStackTrace()
    }
  }
}
