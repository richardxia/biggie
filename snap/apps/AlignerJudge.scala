package snap.apps

import snap._
import scala.math._
import spark.SparkContext
import SparkContext._

object AlignerJudge {
  def main(args: Array[String]) : Unit = {
    args match {
      case Array("scoreAligner", master, accessKey, secretKey, bfaResultsBucket, alignerResultsSam, maxDist, maxSeparation, minSeparation, clipMax) => {
        val samLines = loadSam(alignerResultsSam)
        
        val masterStr = 
          master match {
            case "None" => None
            case host => Option(host)
          }
        
        val readPairsWithBfaResults = loadReadPairsAndBfaResults(masterStr, accessKey, secretKey, None /* use default base bucket */, bfaResultsBucket)
        
        val verified = scala.collection.mutable.Set[String]()
        
        var nTotal = 0
        var nConfident = 0
        var nError = 0
        
        samLines.foreach(l => {
          //println(l)
          val entry = new RichSamEntry(l)
          
          if (includeRead(entry, maxSeparation.toInt, minSeparation.toInt, clipMax.toInt)) {
            nTotal += 1
            if (!entry.readIsUnmapped) {
              nConfident += 1
              
              val errorType = categorizeError(entry, maxDist.toInt, maxSeparation.toInt, minSeparation.toInt)

              if (errorType == 1) {
                //println(l)
                nError += 1
                //println("Skipping verification b/c error type is " + errorType + ".")
              }
              else if (errorType == 2) {
                println(l)
                nError += 1
                if (!verified.contains(entry.illuminaID)) {
                  verifyError(entry, getReadPairAndBfaResult(entry, readPairsWithBfaResults), maxDist.toInt)
                  verified += entry.illuminaID
                } else
                  println("Already verified.")
              } else if (errorType != 0) {
                println(l)
                println("Unrecognized error type.")
              }
            }
          }
        })
        
        println("% Aligned: " + nConfident.toDouble / nTotal * 100)
        println("% Error: " + nError.toDouble / nTotal * 100)
      }
      case Array("checkBfaErrors", accessKey, secretKey, bfaResultsBucket, snapFile1, snapFile2) => {
        /*
        val samEntryList = List(
          "chr17_42258055_42258349_0:S2000-907_85:5:2206:1250:40803/122 97  chr3  129081156 0 41=1X13=1X4=1X1=1X22=15S  chr3  129081911 0 CTGGCCGGGTGGGGGGCTGACCCCCCACCTCCCTCCCGGACAGGGCGGCTGGCCGGGCAGTGTGGCTCCTCACTTCCCAGTAGGGGCGGCTGGGGGGAGG  ?@@FFFDDH?@FFGEDDDBD?CADDD6B@<AABCCB?651B@@BBDB@888@?A050005+3(3+2<A?33:A:>44?>@AC@:@###############",
          "chr16_46428350_46428519_0:S2000-907_85:8:2305:1802:45543/1/2     99     chr16     46429059     60     4=1X8=1X10=1X20=1X16=1X5=1X7=1X1=1X2=1X14=1X3=     chr16     46428519     -640     GGAACAGAATGGAGTCATCGAATGAAATCGAATGGAATCATCATCAAATGAAATCAAATGGAGTCATCGTATGGACTCCAACGGAATCATCATCGACTGG     CCCFFDDFGHFFFIFHIIJJJJJJJIIIIIIIHIJJJJJJJJJJJIJJIJJJIGHHJIIIC>DBAGHGIFICHA;ACEEE>BE>BCD<CCCCCDCC@DBB",
          "chr2_234085021_234085330_0:S2000-907_85:5:2106:19432:25244/22     163     chr2     234085041     29     20S60M2D3M2D12M5S     =     234085330     390     GGGTGGGTAGGTGGGTGGATGGATGGATGGATGGATGGATGTGTGGGTGGGTGGAGGGATGGATAAATGGATGGATGGATGTGTGGGTGGGTGGAGGGAT     ?@BADDF8ADF:CEHBCG?DGIDGHJIJJJ9DHGBGCH8CF)==CHI5@EEB=B=BDD;?BD:@AACCDDCCC@9>AB3:@(:4<BD55><09@(8<@##     XT:A:M     NM:i:12     SM:i:29     AM:i:29     XM:i:8     XO:i:2     XG:i:4     MD:Z:8A3A20A1T8G0G8G5^AG3^GA3A8",
          "chr16_33865993_33866174_0:S2000-907_85:7:1102:18208:45583/12     83     chrY     13447775     17     100M     =     13447610     -265     ATTCCATTCCTTTCTTTTGACAGGGTATCATTGTGTCACTGAGGCTGGAGTACAGTGGCACAATCTCAGCTCACATTGCATGTCAACTTTCCATTTCATT     CDDDCCCEDECEFFFFFHHFHHHHECHIGIIIJJIIHFIJIJIIGHEJJJIIJIIJIHGDIJIGJJIJJJIIFIJJJIJJJJIIIGJGHHGHFFFDFC@C     XT:A:U     NM:i:5     SM:i:17     AM:i:17     ",
          "chr16_33865993_33866174_0:S2000-907_85:7:1102:18208:45583/22     163     chrY     13447610     17     36S30M4D1M1D33M     =     13447775     265     AATCTCAGCTCACATTTCTTTTCACCATTCCACTGTATTCCATTCTATCCCACTCCATTCCACTCCACTCCACTCCACTCCATTCCACTCCATCCCATTC     C@@FFFFFHGHHHGIIJJIJJJJJJJJIJJIGIJJIJJIIHJJIJJJIHGFHHIGIGIJJHIJIJIHIJJIIFIJJJJJJJJGGAEHFHFFFFFEDCECF     XT:A:M     NM:i:11     SM:i:17     AM:i:17     XM:i:6     XO:i:2     XG:i:5     MD:Z:9C6T3G9^TGCG1^G15C6G3T6",
          "chr17_22254554_22254752_0:S2000-907_85:5:1306:7533:124184/12     99     chr17     22252175     12     100M     =     22252373     289     CTCAGAAACTTCTCTGTGATGATTGCATTCAACTCACAGAGTTGAACCCTCCTATGGATAGAGCAGTGTTGAAACTCTCTTTTTGTGGAATCTGCAAGTG     CCCFFFFFHGHHHJJJGHIJJJJJJJJIIJJJJJJJJJIJJHIJIJJHHIIJJGIJJJJJJIJFJJHIJIIGGHIHGIIGHHHHFBFDCDECCDDDDC;@     XT:A:U     NM:i:4     SM:i:12     AM:i:12     X0:i:1     X1:i:2     XM:i:4     XO:i:0     XG:i:0     MD:Z:36A2C27T5T26     XA:Z:chr17,+22247420,100M,5;chr17,+22245041,100M,5;",
          "chr17_22254554_22254752_0:S2000-907_85:5:1306:7533:124184/22     147     chr17     22252373     12     91M9S     =     22252175     -289     TTCAACTCCCAGAGTTTCACGTTGCTTTTCATAGAGTAGTTCTGAAACATGCTTTTCGTAGTGTCTGCAAGTGGACATTTGGAGCGCTTTCAGGCCTGTG     CADC@BADDDDDDDDDADDDEEEEFFFFFFHHHHHHJJJJJJIJIHFGIIJJJIGJJJIJJJJIJJJJJIJJJIFJJJIIJJIJJJJHHHHHFFFFFCCC     XT:A:M     NM:i:9     SM:i:12     AM:i:12     XM:i:9     XO:i:0     XG:i:0     MD:Z:8A17A9C2A1G0A18G22A1A4",
          "chr19_27731991_27732158_0:S2000-907_85:6:2104:14956:91377/12     99     chr19     27731818     29     14M1I85M     =     27731991     273     TTGTGGAATTGGCAGGGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAGTGATTCTCAGAAACT     BCCFFFDFHHHHGFGHJ2CGFGFGHJHIJIIGGIJIIIIIDDGHIII=CFHGJJIJIJGHHHHHHFFCFFEEEEEDDDDDDDDDD:ACDDEDDDCCDCDD     XT:A:M     NM:i:11     SM:i:29     AM:i:29     XM:i:10     XO:i:1     XG:i:1     MD:Z:10T0C2A16A5A6C11G18A7A5T9",
          "chrY_13842655_13842675_0:S2000-907_85:6:1207:17609:18687/12     83     chrY     13842760     12     81M19S     =     13842655     -186     ATGGAATGGAATGGAATGCAATGGAATGGAATCAACCCGAGTGCAATGGAATGGAGAGGAATGGAATGGAATGGAATGGAAACTACCCGAATGGAATGGA     B:DFFEEEHHHEEHFHAIIIIJIHGJHHJIHFEGIIHJJJIJJJJIJIJJIIJJIIGJJJIIIGEEJJJJJJJJJJIJIEHHHFCJJHGHGHFFFFFCCC     XT:A:M     NM:i:6     SM:i:12     AM:i:12     XM:i:6     XO:i:0     XG:i:0     MD:Z:18G18T5G11A0T1C22"//,
          //"chr18_18515511_18516193_0:S2000-907_85:5:1108:12529:133864     147     chr18     18515680     17     99M1S     =     18515511     -268     ACCTAGACAGAAGCACTCTCAGAAACTACTTTGTGATAGCTGCATTGATATCAGAGAGTTGAATATTCCCTTTCTAAGGGCAGGCTTGAAAGCGTCTTTT     ##@@>;3A;@>A@@66(.3FDFDE?7)77.7;@:===.)CFC=8)8=0**9D??9*AGGHIGDB3F?@@F?B4G>GIIGHBGC3E?BFA2)8?BBBD@??     XT:A:M     NM:i:18     SM:i:17     AM:i:17     XM:i:18     XO:i:0     XG:i:0     MD:Z:12A2T11G4A4G0T7C1C0G3C8T0C10A0T2A5T6C1C5",
          //"chr10_42529309_42529465_0:S2000-907_85:7:2104:12601:137806     163     chr10     42529124     15     100M     =     42529309     285     AAACTGCTCTATGAAAAGAAAGGTTAAACTCTGTGAGTTGAACGCACACACCACAAAGTATTTGTTGAGAATGATTCTGTGTAGTTTTTATACGAAGATA     CCCFFFFFHHHHHJJJJJJJJIJJJJJJJJJJJHIJJJJIJJJJJIHIJJJJJJIJJJCFHIJJJJIJJHHHHHFFFFFFFFFECDEEDDDEDDDDDDDD     XT:A:M     NM:i:14     SM:i:15     AM:i:15     XM:i:14     XO:i:0     XG:i:0     MD:Z:5C11T24T0A6T0A0G7G2T0C15C0C8C1A7"
        )
        */
        val samEntryList = List("chr2_234085021_234085330_0:S2000-907_85:5:2106:19432:25244/22     163     chr2     234085041     29     20S60M2D3M2D12M5S     =     234085330     390     GGGTGGGTAGGTGGGTGGATGGATGGATGGATGGATGGATGTGTGGGTGGGTGGAGGGATGGATAAATGGATGGATGGATGTGTGGGTGGGTGGAGGGAT     ?@BADDF8ADF:CEHBCG?DGIDGHJIJJJ9DHGBGCH8CF)==CHI5@EEB=B=BDD;?BD:@AACCDDCCC@9>AB3:@(:4<BD55><09@(8<@##     XT:A:M     NM:i:12     SM:i:29     AM:i:29     XM:i:8     XO:i:2     XG:i:4     MD:Z:8A3A20A1T8G0G8G5^AG3^GA3A8")
        
        val samEntries = samEntryList.map(e => new RichSamEntry(e))
        
        // Load brute force results
        val readPairsWithBfaResults = loadReadPairsAndBfaResults(None, accessKey, secretKey, None, "NA12878_2000Reads_confDiff1_v2")
        val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
        val (reads1, reads2) = readPairs.unzip
        
        // Load snap results
        println("Loading SNAP results...")
        var fileIn = new java.io.FileInputStream(snapFile1)
        var in = new java.io.ObjectInputStream(fileIn)
        val snapResults1 = in.readObject.asInstanceOf[Array[AlignResult]]
        val readsWithSnapResults1 = reads1.zip(snapResults1)

        fileIn = new java.io.FileInputStream(snapFile2)
        in = new java.io.ObjectInputStream(fileIn)
        val snapResults2 = in.readObject.asInstanceOf[Array[AlignResult]]
        val readsWithSnapResults2 = reads2.zip(snapResults2)

        val readsWithSnapResults = readsWithSnapResults1.zip(readsWithSnapResults2)

        val maxDist = 10
        val separationDist = 1000
        val confDiff = 1
        
        // Try alignment
        samEntries.foreach(e => {
          println(e.ID)
          val readPairWithBfaResult = getReadPairAndBfaResult(e, readPairsWithBfaResults)
          val (readPair, bfaResult) = readPairWithBfaResult
          val (read1, read2) = readPair
          println("bfaResult: " + bfaResult)

          val readPairWithSnapResult = getReadPairAndSnapResult(e, readsWithSnapResults)
          val ((read1a, snapRes1), (read2a, snapRes2)) = readPairWithSnapResult

          var snapUpperBound = 2 * maxDist
          var bestPairScore = BruteForceAligner.verifyPairedHit(snapRes1, snapRes2, separationDist)
          println("bestPairScore: " + bestPairScore)
          if (bestPairScore < snapUpperBound)
            snapUpperBound = bestPairScore

          // NOTE:  brute force doesn't seem to be doing anything smart for "snapFirstReadIsRC" => potential spot for optimization but shouldn't affect correctness
          BruteForceAligner.pairedAlign(read1.data, read2.data, snapUpperBound, maxDist, confDiff, separationDist, false)
          //println("snapUpperBound: " + snapUpperBound)
          println
        })
      }
      case Array("genomePrint") => {
        // chr2_234085021_234085330_0
        val ref = new Array[Byte](122)
        val chr2Piece = GenomeLoader.genome.getPiece("chr2")
        GenomeLoader.genome.getSubstring(chr2Piece.startIndex + 234085020, chr2Piece.startIndex + 234085141, ref)
        
        val rc = new Array[Byte](122)
        DNA.getReverseComplement(ref, rc)
        
        println("Forward: " + new String(ref))
        println("Reverse: " + new String(rc))
      }
      case Array("convertPos", absPos) => {
        val (myPiece, myOffset) = GenomeLoader.genome.getLocation(absPos.toLong)
        println(myPiece + "(" + myOffset + ")")
      }
      case _ =>
    }
  }

  def loadSam(fname: String): List[String] = {
    val source = scala.io.Source.fromFile(fname)
    val lines = source.getLines

    lines.filter(l => !l.startsWith("@")).toList
  }

  // Note that in the future, I might have to pass in the sc instead of creating it inside this function if I'm using a sc elsewhere
  def loadReadPairsAndBfaResults(master: Option[String], accessKey: String, secretKey: String, baseBucketName: Option[String], subBucketName: String): Array[((Read, Read), AlignResult)] = {
    val baseBucket = 
      baseBucketName match {
        case Some(name) => name
        case None => "hashingGenome"
      }
    
    val s3Path = "s3n://" + accessKey + ":" + secretKey + "@" + baseBucket + "/" + subBucketName

    val scDest = 
      master match {
        case Some(host) => "1@" + host + ":5050"
        case None => "local[8]"
      }

    val sc = new SparkContext(scDest, "bruteForceAligner")
    val readsWithBfaResults = sc.objectFile[((Read, Read), AlignResult)](s3Path).collect
    
    readsWithBfaResults
  }

  def getReadPairAndBfaResult(entry: RichSamEntry, readPairsWithBfaResults: Array[((Read, Read), AlignResult)]): ((Read, Read), AlignResult) = {
    val (readPairs, bfaResults) = readPairsWithBfaResults.unzip
    val (reads1, reads2) = readPairs.unzip
    
    val read1 = reads1.filter(r => {
      // get read's ID prefix (not the 1/2 suffix that indicates whether it's the first or 2nd read)
      val rID = r.idStr
      val rPrefix = rID.split("/")(0)

      // get corresponding section of entry's ID (from Illumina, not the wgsim-style ID string)
      val ePrefix = entry.illuminaID

      rPrefix == ePrefix
    }).head

    // pull out corresponding bfa result
    val readNum = reads1.indexOf(read1)
    val read2 = reads2(readNum)
    val bfaResult = bfaResults(readNum)
    
    ((read1, read2), bfaResult)
  }
  
  def getReadPairAndSnapResult(entry: RichSamEntry, readsWithSnapResults: IndexedSeq[((Read, AlignResult), (Read, AlignResult))]): ((Read, AlignResult), (Read, AlignResult)) = {
    val (reads1WithSnapResults, reads2WithSnapResults) = readsWithSnapResults.unzip
    val (reads1, snapResults1) = reads1WithSnapResults.unzip
    val (reads2, snapResults2) = reads2WithSnapResults.unzip
    
    val read1 = reads1.filter(r => {
      // get read's ID prefix (not the 1/2 suffix that indicates whether it's the first or 2nd read)
      val rID = r.idStr
      val rPrefix = rID.split("/")(0)

      // get corresponding section of entry's ID (from Illumina, not the wgsim-style ID string)
      val ePrefix = entry.illuminaID

      rPrefix == ePrefix
    }).head

    // pull out corresponding bfa result
    val readNum = reads1.indexOf(read1)
    val read2 = reads2(readNum)
    val snapResult1 = snapResults1(readNum)
    val snapResult2 = snapResults2(readNum)
    
    ((read1, snapResult1), (read2, snapResult2))
  }

  def includeRead(entry: RichSamEntry, maxSeparation: Int, minSeparation: Int, clipMax: Int): Boolean = {
    // don't include read if:
    // bfa placed reads too close together
    // bfa placed reads too far apart
    val bfaPos = getBfaPos(entry.ID)
    val lowEnd = bfaPos.min
    val highEnd = bfaPos.max
    val numClipped = entry.read.qualityStr.indexWhere(c => c != '#')
    
    if (highEnd - lowEnd > maxSeparation || highEnd - lowEnd < minSeparation) false
    else if (numClipped > clipMax) false
    else true
  }

  def categorizeError(entry: RichSamEntry, maxDist: Int, maxSeparation: Int, minSeparation: Int): Int = {
    var errorNum = 0
    //val samPos = Set(entry.readPos, entry.pairPos)
    val bfaPos = getBfaPos(entry.ID)
    
    //val pos1Correct = (abs(samPos.min - bfaPos.min) <= maxDist || abs(samPos.min - bfaPos.max) <= maxDist)
    //val pos2Correct = (abs(samPos.max - bfaPos.min) <= maxDist || abs(samPos.max - bfaPos.max) <= maxDist)

    // count how many #'s the quality string starts with
    val numClipped = entry.read.qualityStr.indexWhere(c => c != '#')

    // TODO: cases are correct (0), error but don't verify (1), error but verify (2)
    if (entry.readPiece == getBfaPiece(entry.ID) && (abs(entry.readPos - bfaPos.min) <= maxDist || abs(entry.readPos - bfaPos.max) <= maxDist)) errorNum = 0  // aligner & bfa agree on both pos
    else if (entry.pairSeparation > maxSeparation) errorNum = 1 // aligner placed reads too far apart
    else if (entry.pairSeparation < minSeparation) errorNum = 1 // aligner placed reads too close together
    else errorNum = 2
    
    errorNum
  }
  
  def getBfaPos(ID: String): Set[Long] = {
    val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r
    val ID_REGEX(piece, start, end) = ID
    Set(start.toLong, end.toLong)
  }
  
  def getBfaPiece(ID: String): String = {
    val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r
    val ID_REGEX(piece, start, end) = ID
    piece
  }
  
  def verifyError(entry: RichSamEntry, readPairWithBfaResult: ((Read, Read), AlignResult), maxDist: Int): (Int, Int) = {
    val (readPair, bfaResult) = readPairWithBfaResult
    val (read1, read2) = readPair
    val readLen = read1.dataStr.length
    
    var errorScore = -1
    var bfaScore = -1
    
    // get error score
    // align read1 at its pos, with right orientation
    // align read2 at its pos, with right orientation
    // get combined score
    
    // I have one sam entry representing a pair
    // I have both reads in the read pair
    // I want to figure out which of the pair this read corresponds to (and use the read pair to get the other one)
    // Once I know whether read 1 or 2 is the one indicated in the pair, I know which position & orientation to use
    // If read1 is entry.read, align read1 against readPos & get orientation from readIsForward
    // If read2 is the pair, align read2 against pairPos & with opposite orientation
    // (vice versa if read2 is entry.read)
    val rc1 = new Array[Byte](readLen)
    DNA.getReverseComplement(read1.data, rc1)
    
    val rc2 = new Array[Byte](readLen)
    DNA.getReverseComplement(read2.data, rc2)
    
    val readPiece = GenomeLoader.genome.getPiece(entry.readPiece)
    val pairPiece = GenomeLoader.genome.getPiece(entry.pairPiece)
    
    if (entry.read.dataStr == read1.dataStr || entry.read.dataStr == new String(rc1)) {
      val score1 = alignRead(read1, entry.readIsRC, readPiece.startIndex + entry.readPos, maxDist)
      val score2 = alignRead(read2, !entry.readIsRC, pairPiece.startIndex + entry.pairPos, maxDist)
      errorScore = score1 + score2
    } else if (entry.read.dataStr == read2.dataStr || entry.read.dataStr == new String(rc2)) {
      val score1 = alignRead(read1, !entry.readIsRC, pairPiece.startIndex + entry.pairPos, maxDist)
      val score2 = alignRead(read2, entry.readIsRC, readPiece.startIndex + entry.readPos, maxDist)
      errorScore = score1 + score2
    } else {
      println("read1:    " + read1.dataStr)
      println("read2:    " + read2.dataStr)
      println("sam read: " + entry.read.dataStr)
      println("Mismatch between read pair & read in SAM entry.")
      return (-1, -1)
    }

    bfaResult match {
      case RichPairSingleHit(pos1, pos2, score, isRC) => {
        // get bfa score
        bfaScore = score
      }
      case _ => println("Didn't get single hit in BFA.")
    }

    val relativeStr = 
      if (errorScore >= 0 && errorScore < bfaScore) "worse"
      else "better"

    println("BFA score (" + bfaScore + ") is " + relativeStr + " than aligner score (" + errorScore + ")")

    (errorScore, bfaScore)
  }
  
  def alignRead(read: Read, isRC: Boolean, pos: Long, maxDist: Int): Int = {
    val lv = new LandauVishkin(maxDist)
    val readLen = read.data.length
    val ref = new Array[Byte](readLen)
    GenomeLoader.genome.getSubstring(pos, pos + readLen, ref)

    val readToCheck = 
      if (isRC) {
        val rc = new Array[Byte](readLen)
        DNA.getReverseComplement(read.data, rc)
        rc
      } else read.data
    
    lv.distance(ref, readLen, readToCheck, readLen, maxDist)
  }
  
  def testGetBfaPos = {
    val IDs = List("chr1_1794501_1794685_0:S2000-907_85:5:2108:17731:93158", "chrUn_gl000220_126043_126195_0:S2000-907_85:8:2305:14907:33559/12")
    val posSets = List(Set(1794501, 1794685), Set(126043, 126195))
    (0 until IDs.length).foreach(i => {
      val id = IDs(i)
      val posSet = posSets(i)
      assert(getBfaPos(id) == posSet)
    })
  }

  def testAlignRead = {
    val testPos = GenomeLoader.genome.getPiece("chr2").startIndex + 10000000
    val readLen = 100
    val ref = new Array[Byte](readLen)
    GenomeLoader.genome.getSubstring(testPos, testPos + readLen, ref)
    
    assert(alignRead(new Read("@".getBytes, ref, "2".getBytes), false, testPos, 10) == 0)
    
    val rc = new Array[Byte](readLen)
    DNA.getReverseComplement(ref, rc)
    
    assert(alignRead(new Read("@".getBytes, rc, "2".getBytes), true, testPos, 10) == 0)
    
  }

  def testCategorizeError = {
    var samLine = "chr2_165437074_165437289_0:S2000-907_85:6:2208:17294:119490	97	chr9	124020010	37	100M	chr2	165437190	0	GTGTTGGTGTGCACCTGTAATCTCAGCTACTCAGGAGGCTGAGGCAGGAGAATCACTTGAACCTGGGATGTGGAAGTTGCAGTGAGCCAAGATTGTGCCA	B@@DFFFFHHHHHJJJJJJJJJJJJJJIJJJJJJJHJJJJJIJIIIIGDGHJJJGIIJIJGIIFGHHFDB;BCEEEC;@C@C;@;@CDC<<<::>(4>@C	XT:A:U	NM:i:2	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:2	XO:i:0	XG:i:0	MD:Z:95C1T2"
    assert(categorizeError(new RichSamEntry(samLine), 10, 1000, 10) == 1)
    
    samLine = "chr1_121485067_121485329_0:S2000-907_85:6:1101:14798:3263/1	97	chr1	121485230	37	100M	=	121485068	-6002	GAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTGACAGAGTTCAACCTTTCTTTTCATAGAGCA	CCCFFFFFHHHHHJJJJJJJJIIJJJIJIJJJIJJJJJIJJJJJJJJJJJJGHCGIGIAGGHIJJJEHGJIJIIIHIHIHHHHHHFFFFFFECEDCCCDC	XT:A:U	NM:i:4	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:4	XO:i:0	XG:i:0	MD:Z:11G9A47C8T21"
    assert(categorizeError(new RichSamEntry(samLine), 10, 1000, 10) == 2)
    
    samLine = "chr1_121485067_121485329_0:S2000-907_85:6:1101:14798:3263/1	97	chr1	121485230	37	100M	=	121485068	-6	GAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTGACAGAGTTCAACCTTTCTTTTCATAGAGCA	CCCFFFFFHHHHHJJJJJJJJIIJJJIJIJJJIJJJJJIJJJJJJJJJJJJGHCGIGIAGGHIJJJEHGJIJIIIHIHIHHHHHHFFFFFFECEDCCCDC	XT:A:U	NM:i:4	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:4	XO:i:0	XG:i:0	MD:Z:11G9A47C8T21"
    assert(categorizeError(new RichSamEntry(samLine), 10, 1000, 10) == 3)
    
    samLine = "chr1_121485068_121485329_0:S2000-907_85:6:1101:14798:3263/1	97	chr1	121485230	37	100M	=	121485068	-62	GAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTGACAGAGTTCAACCTTTCTTTTCATAGAGCA	CCCFFFFFHHHHHJJJJJJJJIIJJJIJIJJJIJJJJJIJJJJJJJJJJJJGHCGIGIAGGHIJJJEHGJIJIIIHIHIHHHHHHFFFFFFECEDCCCDC	XT:A:U	NM:i:4	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:4	XO:i:0	XG:i:0	MD:Z:11G9A47C8T21"
    assert(categorizeError(new RichSamEntry(samLine), 10, 1000, 10) == 4)

    samLine = "chr1_121485067_121485329_0:S2000-907_85:6:1101:14798:3263/1	97	chr1	121485230	37	100M	=	121485068	-62	GAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTGACAGAGTTCAACCTTTCTTTTCATAGAGCA	CCCFFFFFHHHHHJJJJJJJJIIJJJIJIJJJIJJJJJIJJJJJJJJJJJJGHCGIGIAGGHIJJJEHGJIJIIIHIHIHHHHHHFFFFFFECEDCCCDC	XT:A:U	NM:i:4	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:4	XO:i:0	XG:i:0	MD:Z:11G9A47C8T21"
    assert(categorizeError(new RichSamEntry(samLine), 10, 1000, 10) == 5)
    
    samLine = "chr1_121485068_121485230_0:S2000-907_85:6:1101:14798:3263/1	97	chr1	121485230	37	100M	=	121485068	-62	GAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTGACAGAGTTCAACCTTTCTTTTCATAGAGCA	CCCFFFFFHHHHHJJJJJJJJIIJJJIJIJJJIJJJJJIJJJJJJJJJJJJGHCGIGIAGGHIJJJEHGJIJIIIHIHIHHHHHHFFFFFFECEDCCCDC	XT:A:U	NM:i:4	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:4	XO:i:0	XG:i:0	MD:Z:11G9A47C8T21"
    assert(categorizeError(new RichSamEntry(samLine), 10, 1000, 10) == 0)
  }
  
  def testVerifyError = {
    val readLen = 100

    // create a sam entry
    val samLine = "chr1_121485068_121485230_0:S2000-907_85:6:1101:14798:3263/1	97	chr1	121485230	37	100M	=	121485068	-62	GAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTGACAGAGTTCAACCTTTCTTTTCATAGAGCA	CCCFFFFFHHHHHJJJJJJJJIIJJJIJIJJJIJJJJJIJJJJJJJJJJJJGHCGIGIAGGHIJJJEHGJIJIIIHIHIHHHHHHFFFFFFECEDCCCDC	XT:A:U	NM:i:4	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:4	XO:i:0	XG:i:0	MD:Z:11G9A47C8T21"
    val entry = new RichSamEntry(samLine)

    // create a read pair & its corresponding bfaResult
    val read1 = entry.read

    val pos2 = GenomeLoader.genome.getPiece("chr1").startIndex + entry.pairPos - 1
    val ref = new Array[Byte](readLen)
    GenomeLoader.genome.getSubstring(pos2, pos2 + readLen, ref)
    val rc = new Array[Byte](readLen)
    DNA.getReverseComplement(ref, rc)

    val read2Data = 
      if (entry.readIsRC) ref
      else rc

    val read2 = new Read("@".getBytes, read2Data, "2".getBytes)
    val score = alignRead(read1, entry.readIsRC, entry.readPos, 10)
    val bfaResult = RichPairSingleHit(entry.readPos, pos2, score, entry.readIsRC)
    
    val (errorScore, bfaScore) = verifyError(entry, ((read1, read2), bfaResult), 20)
    assert(errorScore > bfaScore)
  }

  def testGetReadPairAndBfaResult(accessKey: String, secretKey: String) = {
    val readPairsWithBfaResults = loadReadPairsAndBfaResults(None, accessKey, secretKey, None, "NA12878_2000Reads_confDiff1_v2")
    
    val samLine = "chr1_121485067_121485329_0:S2000-907_85:6:1101:14798:3263/1	97	chr1	121485230	37	100M	=	121485068	-62	GAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTGACAGAGTTCAACCTTTCTTTTCATAGAGCA	CCCFFFFFHHHHHJJJJJJJJIIJJJIJIJJJIJJJJJIJJJJJJJJJJJJGHCGIGIAGGHIJJJEHGJIJIIIHIHIHHHHHHFFFFFFECEDCCCDC	XT:A:U	NM:i:4	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:4	XO:i:0	XG:i:0	MD:Z:11G9A47C8T21"
    val entry = new RichSamEntry(samLine)
    
    val ((read1, read2), bfaResult) = getReadPairAndBfaResult(entry, readPairsWithBfaResults)
    
    val id = "S2000-907_85:6:1101:14798:3263"
    assert(read1.idStr.split("/")(0) == id && read2.idStr.split("/")(0) == id)
  }
}