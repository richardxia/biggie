package biggie

class SamRegionReader(samFile: String, region: Range) {
  /* regions should be non-empty array of non-overlapping Ranges in ascending
   * order */
  val file = SAM.read(samFile)
  val reads = file.filter(read => rangesIntersect(region, read))

  private def rangesIntersect(range: Range, entry: SAMEntry): Boolean = {
    return ! ( range.start > entry.position + entry.sequence.size ||
               range.end < entry.position )
  }
}
