package snap
import java.util.Random
import it.unimi.dsi.fastutil.longs.LongArrayList
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import scala.math._

class SimilarRegionDetector (genome: Genome, index: Index, seedsToTry: Int, maxDist: Int, similarityThreshold: Int, readLen: Int = 100, confDiff: Int = 3, maxHits: Int = 40) {
  private val seedLen = index.seedLen
  private val rand = new Random(42)
  private val lv = new LandauVishkin(maxDist + confDiff - 1)
  private val fwdHits = new LongArrayList
  private val rcHits = new LongArrayList
  private val fwdSeen = new LongOpenHashSet
  private val rcSeen = new LongOpenHashSet

  // Byte arrays for the read's reverse complement and any genome substrings we compare with
  private val MAX_READ_LEN = 512 // Not a hard constraint; just for sizing the arrays
  private val rc = new Array[Byte](MAX_READ_LEN)
  private val ref = new Array[Byte](MAX_READ_LEN)

  private val callDistCounts = new Array[Int](maxDist + confDiff + 1);
  private val returnedDistCounts = new Array[Int](maxDist + confDiff + 1);
  private val bestScores = new Array[Int](maxDist + confDiff + 1);

  var correctOnFirstCount = 0
  var hasBigComponents = 0

  def checkForSimilarHits(read: Read) = {
    DNA.getReverseComplement(read.data, rc)

    var foundCorrectAns = false

    var offset = 0
    //var seedNum = 1
    //var offset = ((seedNum / (seedsToTry - 1.0)) * (read.data.length - seedLen)).round.toInt

    fwdHits.clear()
    rcHits.clear()
    var i = 0
    
    // for both forward & rc, get a seed & get hits from that seed

    // check forward
    if (index.get(DNA.substringToLong(read.data, offset, offset + seedLen), fwdHits, maxHits)) {
      // take first hit & see if it is correct
      if (fwdHits.size > 0) {
        val pos = fwdHits.getLong(i) - offset
        if (isCorrect(read.idStr, pos)) {
          correctOnFirstCount += 1
          foundCorrectAns = true
          println("Correct (forward)")
        }
      }
    }
    
    // check rc
    if (index.get(DNA.substringToLong(rc, offset, offset + seedLen), rcHits, maxHits)) {
      // take first hit & see if it is correct
      if (rcHits.size > 0) {
        val pos = rcHits.getLong(i) - offset
        if (isCorrect(read.idStr, pos)) {
          correctOnFirstCount += 1
          foundCorrectAns = true
          println("Correct (reverse)")
        }
      }
    }
    
    // check the rest of the hits - any similar regions?
    if (!foundCorrectAns) {
      var seedNum = 1
      while (seedNum < seedsToTry) {
        index.get(DNA.substringToLong(read.data, offset, offset + seedLen), fwdHits, maxHits)
        index.get(DNA.substringToLong(rc, offset, offset + seedLen), rcHits, maxHits)
        
        // Try next seed
        seedNum += 1
        offset = ((seedNum / (seedsToTry - 1.0)) * (read.data.length - seedLen)).round.toInt
      }

      // quick approach:  get consensus string & average # of diffs
      /*
      if (fwdHits.size > 0) {
        val avgNumDiffs = getAvgNumDiffs(fwdHits)
        println("Avg # diffs (f): " + avgNumDiffs)
      }
      
      if (rcHits.size > 0) {
        val avgNumDiffs = getAvgNumDiffs(rcHits)
        println("Avg # diffs (r): " + avgNumDiffs)
      }
      */
      
      // check for connected components among fwd, rc hits
      // TODO
      //List[(String, Long)]
      var foundComponents = false
      if (fwdHits.size > 0) {
        println("fwdHits sz: " + fwdHits.size)

        val vertices = getHitsAsStr(fwdHits).zip(getHitsAsLong(fwdHits))
        val components = AdjacencyGraphHandler.findConnectedComponents(vertices, similarityThreshold)
        AdjacencyGraphHandler.printConnectedComponentsSummary(vertices, components, 5, None, similarityThreshold)
        val numBigComponents = components.filter(c => c.size >= 5).size
        println("# big components: " + numBigComponents)
        
        if (numBigComponents > 0) {
          hasBigComponents += 1
          foundComponents = true
        }
      }
      
      if (rcHits.size > 0) {
        println("rcHits sz: " + rcHits.size)
        
        val vertices = getHitsAsStr(rcHits).zip(getHitsAsLong(rcHits))
        val components = AdjacencyGraphHandler.findConnectedComponents(vertices, similarityThreshold)
        AdjacencyGraphHandler.printConnectedComponentsSummary(vertices, components, 5, None, similarityThreshold)
        val numBigComponents = components.filter(c => c.size >= 5).size
        println("# big components: " + numBigComponents)
        
        if (numBigComponents > 0 && ! foundComponents /* avoid double-counting */) hasBigComponents += 1
      }
    }
  }
  
  def isCorrect(idStr: String, pos: Long): Boolean = {
    // get correct pos
    val ID_REGEX = """([^:]+)_(\d+)_(\d+)_\d+:.*""".r
    val ID_REGEX(piece, start, end) = idStr
    val lowEnd = min(start.toLong, end.toLong)
    val highEnd = max(start.toLong, end.toLong)
    
    // get piece, relative pos
    val (myPiece, myOffset) = genome.getLocation(pos)
    
    // check if equal
    //(piece == myPiece && (abs(start.toLong - myOffset) < maxDist || abs(end.toLong - myOffset) < maxDist))
    (myPiece == piece && myOffset >= lowEnd - maxDist && myOffset <= highEnd + maxDist)
  }
  
  def getAvgNumDiffs(hits: LongArrayList): Double = {
    val consensusStr = getConsensusStr(hits)
    
    val hitStrList = getHitsAsStr(hits)
    
    val numDiffs = hitStrList.map(hit => getNumDiffs(hit, consensusStr))
    println(numDiffs)
    
    numDiffs.sum.toDouble / numDiffs.length
  }
  
  def getNumDiffs(str1: String, str2: String): Int = {
    assert(str1.length == str2.length)
    
    var numDiffs = 0
    
    (0 until str1.length).foreach(i => {
      if (str1.charAt(i) != str2.charAt(i)) numDiffs += 1
    })
    
    numDiffs
  }

  def getHitsAsStr(hits: LongArrayList): List[String] = {
    var hitStrList: List[String] = Nil

    // get strings
    (0 until hits.size).foreach(i => {
      hitStrList = genome.substring(hits.getLong(i), hits.getLong(i) + readLen) :: hitStrList
    })
    
    hitStrList.reverse
  }
  
  def getHitsAsLong(hits: LongArrayList): List[Long] = {
    var hitLongList: List[Long] = Nil
    
    // get longs
    (0 until hits.size).foreach(i => {
      hitLongList = hits.getLong(i) :: hitLongList
    })
    
    hitLongList.reverse
  }
  
  def getConsensusStr(hits: LongArrayList): String = {
    val hitStrList = getHitsAsStr(hits)

    // get consensus string
    val consensusStr = Array.fill(readLen)(' ')
    
    val ct = Array.fill(readLen)(Array.fill(4)(0))
    (0 until readLen).foreach(r => {
      
      // get most popular char
      (0 until hitStrList.length).foreach(h => {
        ct(r)(DNA.BASE_TO_CODE(hitStrList(h).charAt(r))) += 1
      })
    })
    
    // figure out which char is most popular at each base
    (0 until readLen).foreach(r => {
      consensusStr(r) = DNA.CODE_TO_BASE(argMax(ct(r)))
    })

    new String(consensusStr)
  }
  
  def argMax(baseCts: Array[Int]): Int = {
    val max = baseCts.max
    val argMax = baseCts.indexOf(max)
    argMax
  }

  def testGetConsensusStr = {
    // create read
    val read = new Read("@chr22_23770102_23770595_1:0:0_3:0:0_1/1".getBytes, "TAACTTAACAATACAATATGGAGTGATTAATTCAAGCTAATGAACATACCTGTTACCTCACTTATTTGGC".getBytes, "2222222222222222222222222222222222222222222222222222222222222222222222".getBytes)
    
    // get hits
    fwdHits.clear
    var seedNum = 0
    var offset = 0
    
    while (seedNum < seedsToTry) {
      index.get(DNA.substringToLong(read.data, offset, offset + seedLen), fwdHits, maxHits)
      
      // Try next seed
      seedNum += 1
      offset = ((seedNum / (seedsToTry - 1.0)) * (read.data.length - seedLen)).round.toInt
    }
    
    println("hits:")
    val hitsList = getHitsAsStr(fwdHits)
    hitsList.foreach(println)
    
    println("consensus:")
    val consensusStr = getConsensusStr(fwdHits)
    println(consensusStr)
  }

  def testGetNumDiffs = {
    val str1 = "AACATAGA"
    val str2 = "AACATATA"
    val str3 = "ACCAGATA"
    
    assert(getNumDiffs(str1, str2) == 1)
    assert(getNumDiffs(str2, str3) == 2)
    assert(getNumDiffs(str1, str3) == 3)
  }

  
}