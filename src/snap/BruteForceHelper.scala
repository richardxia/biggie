package snap

import snap.apps._
import scala.io.Source
import scala.collection.mutable.ArrayBuffer

object BruteForceHelper {
  def main(args: Array[String]) : Unit = {
    args match {
      case Array("prepSoapResults", accessKey, secretKey, bfaBucket, confDiff, soapSam) => {
        // load bfa results
        val readPairsWithBfaResults = AlignerJudge.loadReadPairsAndBfaResults(None, accessKey, secretKey, None, bfaBucket)
        val singleHitReadPairsWithBfaResults = readPairsWithBfaResults.flatMap(t => {
          t._2 match {
            case RichPairSingleHit(_, _, _, _) => List(t)
            case RichPairMultiHit(bestPos1, bestPos2, bestScore, secondBestPos1, secondBestPos2, secondBestScore, isRC) => {
              //if (bestScore < secondBestScore)
              if (bestScore + confDiff.toInt <= secondBestScore)
                List((t._1, RichPairSingleHit(bestPos1, bestPos2, bestScore, isRC)))
              else
                Nil
            }
            case _ => Nil
          }
        })
    
        // load soap results
        val soapEntries = AlignerJudge.loadSam(soapSam).map(samLine => new RichSamEntry(samLine))
        val soapIds = soapEntries.map(samEntry => samEntry.illuminaIDWithSuffix)
    
        // create output file by appending to input fname
        val outFile = soapSam + ".unmapped"
        val out = new java.io.File(outFile)
        if (!out.exists)
          out.createNewFile
      
        val fw = new java.io.FileWriter(out.getName)
        val bw = new java.io.BufferedWriter(fw)
    
        // for each read pair in bfa results
        // present in soap results?
        // if not, create a dummy entry (indicating unmapped)
        // write out dummy entries to file
        singleHitReadPairsWithBfaResults.foreach(t => {
          val (readPair, bfaRes) = t
    	    val (read1, read2) = readPair
	    
    	    // unmapped pair:
    	    //chrY_58985524_58986089_0:@HS2000-907_85:5:1108:19820:91616/12	101	*	0	0	*	*	0	0	GTGATGTAACTCTTGTCTTGGCTCTTCTTACAGGGTTATTGTGA GATATCTCTGCACTGATCACCCAGGTGATGTAACTCTTGTCTAGCCTCTGCCTACA	@?@DDFDDHDFFHGH:EGGIGCEHIIEIIIIIIIICFHCGIDGGIEEGCGHHGHEHE?FGHIIGHHGIA.==CED=C7CCECEEHFEEFBBCC>CAAADC
          //chrY_58985524_58986089_0:@HS2000-907_85:5:1108:19820:91616/22	133	*	0	0	*	*	0	0	GGATCTGCTTACAGGTGGCATTGTGACATATCTCTGCCCTGATC TCCCAGGTGATGTAACATTTGTCTAGGCTCTGCCAAAAGGAGGCATTGTGACATAA	:>><DCAC@@;DDDDCCCCCEDCDCB<DCCEHEJIGFFFC@3<IJIHIHFGFB9IIGHFFIIJIGHCGHIGIGHGJGIIJJJGJIJHBFFFFDDAFDB@?
          
          /*
	        val id1Entries = read1.idStr.split("/")
	        val id1 = id1Entries.toList.slice(0, id1Entries.size - 1).mkString("/")
    	    // is read1 in soap results?
    	    if (!soapIds.contains(id1)) {
    	      val id = Wgsim.getWgsimId(read1.idStr, bfaRes.asInstanceOf[RichPairSingleHit].pos1, bfaRes.asInstanceOf[RichPairSingleHit].pos2, false)
    	      val e = List(id, 101, "*", 0, 0, "*", "*", 0, 0, read1.dataStr, read1.qualityStr).mkString("\t")  // not sure about flags
    	      bw.write(e)
    	      bw.newLine
    	    }
	    
          val id2Entries = read2.idStr.split("/")
          val id2 = id2Entries.toList.slice(0, id2Entries.size - 1).mkString("/")
    	    // is read2 in soap results?
    	    if (!soapIds.contains(id2)) {
    	      val id = Wgsim.getWgsimId(read2.idStr, bfaRes.asInstanceOf[RichPairSingleHit].pos1, bfaRes.asInstanceOf[RichPairSingleHit].pos2, false)
    	      val e = List(id, 133, "*", 0, 0, "*", "*", 0, 0, read2.dataStr, read2.qualityStr).mkString("\t")
    	      bw.write(e)
    	      bw.newLine
    	    }
    	    */
    	    
    	    val entries1 = read1.idStr.split("/")

    	    val id1 = 
      	    if (entries1(entries1.size - 1).length == 1) entries1.toList.slice(0, entries1.size - 1).mkString("/")
      	    else read1.idStr
    	    if (!soapIds.contains(id1)) {
	          val id = Wgsim.getWgsimId(read1.idStr, bfaRes.asInstanceOf[RichPairSingleHit].pos1, bfaRes.asInstanceOf[RichPairSingleHit].pos2, false)
    	      val e = List(id, 101, "*", 0, 0, "*", "*", 0, 0, read1.dataStr, read1.qualityStr).mkString("\t")  // not sure about flags
    	      bw.write(e)
    	      bw.newLine
    	    } else {
	          val id = soapIds(soapIds.indexOf(id1))
	          val entries = id.split("/")
	          if (entries.length == 1) {
              // means that there is no suffix => can't tell if this is 1st or 2nd read from id
	            // so, have to check the read
		          val num = soapIds.count(_ == id1)
		          if (num == 1) {
		            val i = soapIds.indexOf(id1)
		            val e = soapEntries(i)
		            val r = e.read
		            if (read1.dataStr != r.dataStr && DNA.reverseComplement(read1.dataStr) != r.dataStr) {
		              // there's only one read with this id in the soap output, and this isn't it (it's read2), so print to file
		              val id = Wgsim.getWgsimId(read1.idStr, bfaRes.asInstanceOf[RichPairSingleHit].pos1, bfaRes.asInstanceOf[RichPairSingleHit].pos2, false)
                  val e = List(id, 101, "*", 0, 0, "*", "*", 0, 0, read1.dataStr, read1.qualityStr).mkString("\t")  // not sure about flags
                  bw.write(e)
                  bw.newLine
		            }
		          }
	          }
	        }
    	    
    	    val entries2 = read2.idStr.split("/")
    	    val id2 = 
    	      if (entries2(entries2.size - 1).length == 1) entries2.toList.slice(0, entries2.size - 1).mkString("/")
    	      else read2.idStr
    	    if (!soapIds.contains(id2)) {
    	      val id = Wgsim.getWgsimId(read2.idStr, bfaRes.asInstanceOf[RichPairSingleHit].pos1, bfaRes.asInstanceOf[RichPairSingleHit].pos2, true)
    	      val e = List(id, 133, "*", 0, 0, "*", "*", 0, 0, read2.dataStr, read2.qualityStr).mkString("\t")
    	      bw.write(e)
    	      bw.newLine
    	    } else {
	          val id = soapIds(soapIds.indexOf(id2))
	          val entries = id.split("/")
	          if (entries.length == 1) {
	            val num = soapIds.count(_ == id2)
		          println("num: " + num)
		          if (num == 1) {
		            val i = soapIds.indexOf(id2)
		            val e = soapEntries(i)
		            val r = e.read
		            if (read2.dataStr != r.dataStr && DNA.reverseComplement(read2.dataStr) != r.dataStr) {
		              val id = Wgsim.getWgsimId(read2.idStr, bfaRes.asInstanceOf[RichPairSingleHit].pos1, bfaRes.asInstanceOf[RichPairSingleHit].pos2, true)
                  val e = List(id, 133, "*", 0, 0, "*", "*", 0, 0, read2.dataStr, read2.qualityStr).mkString("\t")
                  bw.write(e)
                  bw.newLine
		            }
		          }
	          }
	        }
        })
    
        bw.close
      }
      case Array("unmappedSoapResultsToDummySam", /*accessKey, secretKey, bfaBucket, confDiff,*/ soapFile) => {
        // load bfa results
        /*
        val readPairsWithBfaResults = AlignerJudge.loadReadPairsAndBfaResults(None, accessKey, secretKey, None, bfaBucket)
        val singleHitReadPairsWithBfaResults = readPairsWithBfaResults.flatMap(t => {
          t._2 match {
            case RichPairSingleHit(_, _, _, _) => List(t)
            case RichPairMultiHit(bestPos1, bestPos2, bestScore, secondBestPos1, secondBestPos2, secondBestScore, isRC) => {
              //if (bestScore < secondBestScore)
              if (bestScore + confDiff.toInt <= secondBestScore)
                List((t._1, RichPairSingleHit(bestPos1, bestPos2, bestScore, isRC)))
              else
                Nil
            }
            case _ => Nil
          }
        })
        */
        
        // load unmapped 
        val unmappedReads = readUnmappedSoapOutput(soapFile)
        
        // for each, print to console
        // requires corresponding bfa result (so I can make the Wgsim-style ID)
        unmappedReads.foreach(u => {
          /*
          val readPairWithBfaRes = getCorrespondingBfaResult(u.idStr, singleHitReadPairsWithBfaResults)
          val (readPair, bfaRes) = readPairWithBfaRes
          val (read1, read2) = readPair
          
          // make Wgsim-style id
          val id = Wgsim.getWgsimId(u.idStr, bfaRes.asInstanceOf[RichPairSingleHit].pos1, bfaRes.asInstanceOf[RichPairSingleHit].pos2, false)
          
          val first = u.idStr(u.idStr.length - 1) == 1
          if (first)  println(List(id, 101, "*", 0, 0, "*", "*", 0, 0, read1.dataStr, read1.qualityStr).mkString("\t"))  // not sure about flags
          else println(List(id, 101, "*", 0, 0, "*", "*", 0, 0, read2.dataStr, read2.qualityStr).mkString("\t"))  // not sure about flags
          */
          println(List(u.idStr.substring(1, u.idStr.length), 101, "*", 0, 0, "*", "*", 0, 0, u.dataStr, u.qualityStr).mkString("\t"))
        })
      }
      case _ =>
    }
  }
  
  // adds a dummy quality string
  def readUnmappedSoapOutput(file: String): IndexedSeq[Read] = {
    val lines = Source.fromFile(file).getLines
    val reads = new ArrayBuffer[Read]
    for (Seq(id, data) <- lines.filter(_ != "").sliding(2, 2)) {
      reads += new Read(id.getBytes, data.getBytes, 
        Array.fill(id.length)("2").mkString.getBytes)
    }
    reads
  }
  
  def getCorrespondingBfaResult(readId: String, readPairsWithBfaResults: Array[((Read, Read), AlignResult)]): ((Read, Read), AlignResult) = {
    val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
    val (reads1, reads2) = readPairs.unzip
    
    val read1 = reads1.filter(_.idStr == readId).head

    // pull out corresponding bfa result
    val readNum = reads1.indexOf(read1)
    readPairsWithBfaResults(readNum)
  }
}