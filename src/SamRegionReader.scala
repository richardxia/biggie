package biggie

class SamRegionReader(samFile: String, var region: Range) {
  /* regions should be non-empty array of non-overlapping Ranges in ascending
   * order */
  //val file = SAM.read(samFile)
  val file = new java.io.RandomAccessFile(samFile, "r")
  val index = new SamIndexer(samFile + ".sai")
  //val reads = file.filter(read => rangesIntersect(region, read))

  def reads(): Iterator[SAMEntry] = new Iterator[SAMEntry] {
    var pos = region.start
    file.seek(index(pos))
    //val lastIndex = index(region.end)

    def hasNext(): Boolean = pos < region.end

    override def next(): SAMEntry = {
      val entry = SAM.parseEntry(file.readLine())
      pos = entry.position
      entry
    }
  }

  private def rangesIntersect(range: Range, entry: SAMEntry): Boolean = {
    return ! ( range.start > entry.position + entry.sequence.size ||
               range.end < entry.position )
  }
}
