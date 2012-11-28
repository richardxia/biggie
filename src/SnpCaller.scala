package biggie

class SnpCaller(samFile: String, refFile: String, regionStart: Int = 1, regionEnd: Int = 2) {
  // [regionStart, regionEnd) with respect to reference, 1-indexed
  val ref: Array[Byte] = FASTA.read(refFile).pieces(0).data // 0-indexed
  val reader = new SamRegionReader(samFile, regionStart until regionEnd) // 1-indexed
  val size = regionEnd - regionStart
  val baseCount = Array.ofDim[Int](2, 4, size)
  val coverage = Array.ofDim[Int](2, size)
  val baseCalls = new Array[Int](size)

  def run() {
    for (read <- reader.reads) {
      val dir = read.direction
      for (pos <- math.max(read.position, regionStart) until math.min(read.endPosition, regionEnd)) {
        val base = DNA.BASE_TO_CODE(read.sequence.charAt(pos - read.position))
        val regionPos = pos - regionStart
        baseCount(dir)(base)(regionPos) += 1
        coverage(dir)(regionPos) += 1
      }
    }

    for (pos <- regionStart until regionEnd) {
      call(pos)
    }
    baseCalls
  }

  // pos 0 is first base of region which is 1-indexed
  private def call(pos: Int) {
    val regionPos = pos - regionStart
    val totalCoverage = coverage(0)(regionPos) + coverage(1)(regionPos)
    var best = 0
    var bestFrac = 0.0
    var second = 0
    var secondFrac = 0.0
    for (base <- 0 until 4) {
        val frac = if (baseCount(0)(base)(regionPos) == 0 || baseCount(1)(base)(regionPos) == 0) {
          0.0
        } else {
          (baseCount(0)(base)(regionPos) + baseCount(1)(base)(regionPos)).toDouble / totalCoverage
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
    baseCalls(regionPos) = best
    val base = if (bestFrac == 0.0) "N" else DNA.CODE_TO_BASE(best)
    if (base != ref(pos-1).toChar)
      print("$(" + ref(pos-1).toChar + "->" + base + ")")
    else
      print(base)
    //print(" ")
    //print(bestFrac)
    //print(" ")
  }
}

object SnpCaller {
  def main(args: Array[String]) {
    val baseCalls = new SnpCaller(args(0), args(1), args(2).toInt, args(3).toInt).run()
  }
}
