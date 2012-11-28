package snap.apps

object HammingTest {
  def benchmark(name: String)(block: => Unit) {
    for (i <- 0 until 3) {
      val start = System.nanoTime()
      block
      val end = System.nanoTime()
      println("%s: %.2fs".format(name, (end - start) / 1.0e9))
    }
  }

  def hamming1(s1: Array[Byte], off1: Int, s2: Array[Byte], off2: Int, length: Int, limit: Int): Int = {
    if (off1 < 0 || off2 < 0 || length < 0 || off1 + length > s1.length || off2 + length > s2.length) {
      throw new IllegalArgumentException("Bad array bounds")
    }
    var distance = 0
    var i = 0
    while (i < length) {
      if (s1(off1 + i) != s2(off2 + i)) {
        distance += 1
        if (distance > limit) {
          return -1
        }
      }
      i += 1
    }
    return distance
  }

  def hamming2(s1: Array[Byte], off1: Int, s2: Array[Byte], off2: Int, length: Int, limit: Int): Int = {
    if (off1 < 0 || off2 < 0 || length < 0 || off1 + length > s1.length || off2 + length > s2.length) {
      throw new IllegalArgumentException("Bad array bounds")
    }
    var distance = 0
    var i1 = off1
    var i2 = off2
    var end1 = off1 + length
    while (i1 < end1) {
      if (s1(i1) != s2(i2)) {
        distance += 1
        if (distance > limit) {
          return -1
        }
      }
      i1 += 1
      i2 += 1
    }
    return distance
  }

  def hamming2b(s1: Array[Byte], s2: Array[Byte], limit: Int): Int = {
    if (s1.length != s2.length) {
      throw new IllegalArgumentException("Arrays not of equal size")
    }
    var distance = 0
    var i = 0
    while (i < s1.length) {
      if (s1(i) != s2(i)) {
        distance += 1
        if (distance > limit) {
          return -1
        }
      }
      i += 1
    }
    return distance
  }

  def main(args: Array[String]) {
    val s0 = "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG".getBytes
    val s0b = "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG".getBytes
    val s1 = "ACTGACTAACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG".getBytes
    val s2 = "ACTGACTAACTGACTGACTGACTGACTGACTGAATGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG".getBytes
    val s3 = "ACTGACTAACTGACTGACTGACTGACTGACTGAATGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTAACTGACTGACTGACTG".getBytes
    val s4 = "ACGATGCATGTGCATGTCTAGCTGTGCTAGCTAGCTAGTCAGTCGTAGTAGTCAGTGCTGATCGATCGTAGTAGCATGCTAGTCGATGCATGCTAGTAGC".getBytes

    val N = 10000000

    benchmark("hamming1, s0-s0") {
      for (i <- 1 to N) { hamming1(s0, 0, s0b, 0, 100, 3) }
    }
    benchmark("hamming1, s0-s3") {
      for (i <- 1 to N) { hamming1(s0, 0, s3, 0, 100, 3) }
    }
    benchmark("hamming1, s0-s4") {
      for (i <- 1 to N) { hamming1(s0, 0, s4, 0, 100, 3) }
    }

    benchmark("hamming2, s0-s0") {
      for (i <- 1 to N) { hamming2(s0, 0, s0b, 0, 100, 3) }
    }
    benchmark("hamming2, s0-s3") {
      for (i <- 1 to N) { hamming2(s0, 0, s3, 0, 100, 3) }
    }
    benchmark("hamming2, s0-s4") {
      for (i <- 1 to N) { hamming2(s0, 0, s4, 0, 100, 3) }
    }

    benchmark("hamming2b, s0-s0") {
      for (i <- 1 to N) { hamming2b(s0, s0b, 3) }
    }
    benchmark("hamming2b, s0-s3") {
      for (i <- 1 to N) { hamming2b(s0, s3, 3) }
    }
    benchmark("hamming2b, s0-s4") {
      for (i <- 1 to N) { hamming2b(s0, s4, 3) }
    }
  }
}
