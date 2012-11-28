package biggie

class SnpCaller(samFile: String, start: Int = 0, end: Int = 1) {
  val reader = new SamRegionReader(samFile, start until end)
  val size = end - start
  val baseCount = Array.ofDim[Int](2, 4, size)
  val coverage = Array.ofDim[Int](2, size)
  val baseCalls = new Array[Int](size)

  // [start, end) with respect to reference
  def run() {
    //  For each read in region
    //    build hash table of counts
    for (read <- reader.reads) {
      var regionPos = read.position - start
      if (regionPos < 0)
        regionPos = 0
      var max = read.position + read.sequence.size - start
      if (max >= size)
        max = size
      val dir = read.direction
      while (regionPos < max) {
        val base = DNA.BASE_TO_CODE(read.sequence.charAt(regionPos - read.position))
        baseCount(dir)(base)(regionPos) += 1
      }
    }

    for (pos <- 0 until size) {
      call(pos)
    }
    baseCalls
  }

  // pos 0 is first base of region
  private def call(pos: Int) {
    val totalCoverage = coverage(0)(pos) + coverage(1)(pos)
    var best = 0
    var bestFrac = 0.0
    var second = 0
    var secondFrac = 0.0
    for (base <- 0 until 4) {
        val frac = if (baseCount(0)(base)(pos) == 0 || baseCount(1)(base)(pos) == 0) {
          0.0
        } else {
          (baseCount(0)(base)(pos) + baseCount(1)(base)(pos)).toDouble / totalCoverage
        }
        if (frac > bestFrac) {
          second = best
          secondFrac = bestFrac
          best = base
          bestFrac = frac
        } else if (frac > secondFrac) {
          second = base
          secondFrac = frac
        }
    }

    // Just call best for now
    baseCalls(pos) = best
  }
}

object SnpCaller {
  def main(args: Array[String]) {
    val baseCalls = new SnpCaller(args(0), args(1).toInt, args(2).toInt).run()
  }
}
