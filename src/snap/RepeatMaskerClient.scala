package snap

object RepeatMaskerClient {
  def browse(repeatLine: String) = {
    val parser = new RepeatMaskerParser(repeatLine)
    
    // get genome at the location indicated
    val absBeginPos = GenomeLoader.genome.getAbsPos(parser.querySequence, parser.beginInQuery)
    val absEndPos = GenomeLoader.genome.getAbsPos(parser.querySequence, parser.endInQuery)
    GenomeLoader.genome.substring(absBeginPos, absEndPos)
  }
}