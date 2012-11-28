package biggie

object SamIndexer {
  def main(args: Array[String]) {
    if (args.size != 1) {
      println("Usage: SamIndexer file.sam")
      println("Creates a file named file.sam.sai")
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
