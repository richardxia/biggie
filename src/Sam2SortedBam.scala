/**
* input: unsorted sam file (first commandline argument)
* output: sorted bam file (second commandline argument)
* temporary files will created and deleted in the current directory
*/

import java.io._
import java.util.ArrayList
import java.util.concurrent._
import scala.collection.JavaConversions._
import scala.sys.process.Process
import scala.sys.process.ProcessIO

object Sam2SortedBam {

	// globals
	val NUM_LINES = 8700000  // number of lines in each file created by splitting the inSamFile
	val PREFIX = "unsorted"  // prefix for the intermediate files, no other files in the directory should begin with PREFIX
	val NUM_THREADS = 7      // number of threads to use
	val MAX_MEM = 1000000000 // max memory for sorting, NUM_THREADS * MAX_MEM should not exceed total memory 

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
	class FileThread(fileName: String, refName: String) extends Callable[String] {

		def call():String = {

			// convert file to bam
			val pb2 = Process("samtools view -h -b -S -u -t " + refName + " -o " + fileName + ".bam " + fileName)
			val p2 = pb2.run()
			println("sam to bam: " + fileName + " " + p2.exitValue())

			// sort bam file
			val pb3 = Process("samtools sort -m " + MAX_MEM + " " + fileName + ".bam " + fileName + "_sorted")
			val p3 = pb3.run()
			println("sorted bam: " + fileName + " " + p3.exitValue())

			return "finished " + fileName
		}
	}

	def main(args: Array[String]) {

		val inSamFile = args(0)  // unsorted SAM file
		val refName = args(1)    // reference index file (fai format)
		val outBamFile = args(2) // output BAM file

		// split up the sam file
		val pb1 = Process("split -l " + NUM_LINES + " " + inSamFile + " " + PREFIX)
		val p1 = pb1.run()
		println("split sam file " + p1.exitValue())

		// get list of files
		val fileList = getFileList()
		val numFiles = fileList.size()
		var nameString = ""

		// create list of threads (loop over haplotypes so we leave each out in turn)
		val allHapThreads = new ArrayList[FileThread]()
		for (i <- 0 to (numFiles-1)) {
			val fileName = fileList.get(i)
			allHapThreads.add(new FileThread(fileName, refName))
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
		val pb4 = Process("samtools merge " + outBamFile + " " + nameString)
		val p4 = pb4.run()
		println("merged sorted bams " + p4.exitValue())

		// delete temporary files
		val pb5 = Process("rm " + PREFIX + "*")
		val p5 = pb5.run()
		println("deleted temp files " + p5.exitValue())
	}
}