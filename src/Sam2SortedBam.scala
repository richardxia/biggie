/**
* input: unsorted sam file (first commandline argument)
* output: sorted bam file (second commandline argument)
* temporary files will created in the current directory (not deleted right now)
*/

import java.io._
import java.util.ArrayList
import java.util.concurrent._
import scala.collection.JavaConversions._
import scala.sys.process.Process
import scala.sys.process.ProcessIO

object Sam2SortedBam {

	// globals
	val NUM_LINES = 20000000 // number of lines in each file created by splitting the inSamFile
	val PREFIX = "unsorted" // prefix for the intermediate files
	val NUM_THREADS = 4;    // number of threads to use

	// write a set of strings to a new file
	def writeFile(fileName: String, data: Iterator[String]) = {
		val fstream = new FileWriter(fileName)
  	val out = new BufferedWriter(fstream)
  	data.foreach(s => out.write(s + "\n"))
  	out.close()
	}

	// get a list of files starting with PREFIX in the current directory
	def getFileList():ArrayList[String] = {
 
  	val folder = new File(".");
  	val listOfFiles = folder.listFiles()
  	val nameArray = new ArrayList[String]
 
  	for (file <- listOfFiles) {

  		if (file.isFile) {
  			var name = file.getName()
  			if (name.startsWith(PREFIX)) {
  				nameArray.add(file.getName())
  			}
  		}
  	}
  	return nameArray
	}

	// our file processing thread, returns type String
	class FileThread(fileName: String) extends Callable[String] {

		// TODO: currently not deleting temp files
		def call():String = {

			// add header to file if necessary (first file already has it)
			if (fileName != PREFIX + "aa") {
				val pb3 = Process("cat header.txt " + fileName) 
				val pio3 = new ProcessIO(_ => (), stdout => writeFile(fileName + "_header.sam", scala.io.Source.fromInputStream(stdout).getLines), _ => ())
				val p3 = pb3.run(pio3)
				println(p3.exitValue())
			} else {
				val pb3 = Process("mv " + fileName + " " + fileName + "_header.sam")
				val p3 = pb3.run()
				println(p3.exitValue())
			}

			// convert file to bam
			val pb4 = Process("samtools view -b -S -o " + fileName + ".bam " + fileName + "_header.sam")
			val p4 = pb4.run()
			println(p4.exitValue())

			// sort bam file
			val pb5 = Process("samtools sort " + fileName + ".bam " + fileName + "_sorted") // TODO: sort can take in a max memory flag (-m)
			val p5 = pb5.run()
			println(p5.exitValue())

			return "finished " + fileName
		}
	}

	def main(args: Array[String]) {

		val inSamFile = args(0)
		val outBamFile = args(1)

		// create the header file
		val pb0 = Process("samtools view -H -S " + inSamFile)
		val pio0 = new ProcessIO(_ => (), stdout => writeFile("header.txt", scala.io.Source.fromInputStream(stdout).getLines), _ => ())
		val p0 = pb0.run(pio0)
		println(p0.exitValue())

		// split up the sam file
		val pb1 = Process("split -l " + NUM_LINES + " " + inSamFile + " " + PREFIX)
		val p1 = pb1.run()
		println(p1.exitValue())

		// get list of files
		val fileList = getFileList()
		val numFiles = fileList.size()
		var nameString = ""

		// create list of threads (loop over haplotypes so we leave each out in turn)
		val allHapThreads = new ArrayList[FileThread]()
		for (i <- 0 to (numFiles-1)) {
			val fileName = fileList.get(i)
			allHapThreads.add(new FileThread(fileName))
			nameString += fileName + "_sorted.bam "
		}
		
		// run our threads
		val taskExecutor = Executors.newFixedThreadPool(NUM_THREADS)
		val futures = taskExecutor.invokeAll(allHapThreads)
		taskExecutor.shutdown()
		for (i <- 0 to (numFiles-1)) {
			println(futures.get(i).get())
		}

		// merge all the sorted bams
		val pb6 = Process("samtools merge " + outBamFile + " " + nameString)
		val p6 = pb6.run()
		println(p6.exitValue())
	}
}