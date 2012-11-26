package snap

import scala.io.Source

import java.util.concurrent.ConcurrentHashMap

object SAM {
  def parseEntry(line: String): SAMEntry = {
    val fields = Utils.split(line, '\t')
    new SAMEntry(
      fields(0),                  // read ID
      fields(1).toInt,            // flags
      fields(2),                  // piece
      fields(3).toInt,            // position
      fields(4).toInt,            // map quality
      fields(5),                  // cigar
      fields(6),                  // next piece
      fields(7).toInt,            // next position
      if (fields(2) == fields(6)) fields(8).toInt else 0, // template len
      fields(9),                  // sequence
      fields(10)                  // quality
    )
  }

  def read(file: String): Iterator[SAMEntry] = {
    val lines = Source.fromFile(file).getLines
    lines.filter(!_.startsWith("@")).map(parseEntry)
  }

  // Bits in the flags field
  val REVERSE = 0x10
}
