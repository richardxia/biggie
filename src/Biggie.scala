package biggie

import java.io.File
import net.sf.samtools.SAMFileReader
import net.sf.picard.reference.IndexedFastaSequenceFile

object Biggie {
  val usage = """Usage: biggie <command> [args, ...]

    Commands:
    all <alignment.bam> <reference.fa>
    //sort <SAM-file>
    classify <sorted-BAM-file>
    snp <alignment.bam> <reference.fa> <regions.txt>
  """

  def unimplemented() {
    println("Unimplemented")
    System.exit(2)
  }

  def runAll(args: Array[String]) {
    val bamFileName = args(0)
    val refFileName = args(1)
    val bamFile = new SAMFileReader(new File(bamFileName), new File(bamFileName + ".bai"))
    bamFile.setValidationStringency(SAMFileReader.ValidationStringency.SILENT)
    val weirdness = new SimpleClassifier(bamFile).run()
    val refSeq = bamFile.getFileHeader().getSequence(0).getSequenceName()
    val ref = new IndexedFastaSequenceFile(new File(refFileName)).getSequence(refSeq) // 0-indexed
    val regions = 0 to ref.length()
    new SnpCaller(bamFile, ref, refSeq, regions, weirdness).run()
  }

  def runSort(args: Array[String]) {
  }

  def runClassify(args: Array[String]) {
    SimpleClassifier.main(args)
  }

  def runSnp(args: Array[String]) {
    SnpCaller.main(args)
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
