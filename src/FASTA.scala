package biggie

import scala.io.Source

object FASTA {
  def read(file: String): Genome = {
    var curName: String = null
    var curData = new ByteArrayBuilder
    val genome = new Genome
    for (line <- Source.fromFile(file).getLines) {
      if (line.startsWith(">")) {
        if (curName != null) {
          genome.addPiece(curName, curData.toArray)
          curData.clear()
        }
        curName = line.substring(1)
        //curData.append('N'.toByte) // To make pieces 1-indexed  /* for now, want it to be 0-indexed so it's compatible with C++ version of SNAP */
      } else {
        curData.append(line.toUpperCase.getBytes)
      }
    }
    if (curName != null) {
      genome.addPiece(curName, curData.toArray)
      curData.clear()
    }
    genome
  }
}
