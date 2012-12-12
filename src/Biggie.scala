package biggie

import java.io.File
import net.sf.samtools.SAMFileReader

object Biggie {
  val usage = """Usage: biggie <command> [args, ...]

    Commands:
    all <SAM-file>
    //sort <SAM-file>
    classify <sorted-BAM-file>
    snp <sorted-BAM-file>
  """

  def unimplemented() {
    println("Unimplemented")
    System.exit(2)
  }

  def runAll(args: Array[String]) {
    val bamFileName = args(0)
    val bamFile = new SAMFileReader(new File(bamFileName), new File(bamFileName + ".bai"))
    val weirdness = new SimpleClassifier(bamFile).run()
    val regions = 1 until 2
    new SnpCaller(bamFile, new Array[Byte](0), args(1), regions, weirdness).run()
  }

  def runSort(args: Array[String]) {
  }

  def runClassify(args: Array[String]) {
    unimplemented()
  }

  def runSnp(args: Array[String]) {
    unimplemented()
  }

  def invalidCommand(arg: String) {
    println("Invalid command: " + arg)
    println(usage)
    System.exit(3)
  }

  def main(args: Array[String]) {
    if (args.length == 0) {
      println(usage)
      System.exit(1)
    }

    val cmd = args(0)
    val cmdArgs = if (args.length > 0) args.slice(1, args.length)
                  else Array[String]()
    cmd match {
      case "all" => runAll(cmdArgs)
      //case "sort" => runSort(cmdArgs)
      case "classify" => runClassify(cmdArgs)
      case "snp" => runSnp(cmdArgs)
      case _ => invalidCommand(cmd)
    }
  }
}
