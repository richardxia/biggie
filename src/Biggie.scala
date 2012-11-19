object Biggie {
  val usage = """Usage: biggie <command> [args, ...]

    Commands:
    all <SAM-file>
    sort <SAM-file>
    classify <sorted-BAM-file>
    snp <sorted-BAM-file>
  """

  def unimplemented() {
    println("Unimplemented")
    System.exit(2)
  }

  def runAll(args: Array[String]) {
    unimplemented()
  }

  def runSort(args: Array[String]) {
  }

  def runClassify(args: Array[String]) {
    unimplemented()
  }

  def runSnp(args: Array[String]) {
    unimplemented()
  }

  def main(args: Array[String]) {
    if (args.length == 0 || args.length == 1) {
      println(usage)
      System.exit(1)
    }

    val cmd = args(0)
    cmd match {
      case "all" => runAll(args.slice(1, args.length))
      case "sort" => runSort(args.slice(1, args.length))
      case "classify" => runClassify(args.slice(1, args.length))
      case "snp" => runSnp(args.slice(1, args.length))
    }
  }
}
