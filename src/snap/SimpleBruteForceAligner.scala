package snap

import spark.SparkContext
import SparkContext._

object SimpleBruteForceAligner {
  // Note:  for simplicity, I'm not yet using the snapFirstReadIsRC
  // optimization would be to check the direction snap suggests first
  def pairedAlign(read1: Array[Byte], read2: Array[Byte], params: PairedBruteForceParams, snapUpperBound: Int, snapFirstReadIsRC: Boolean, debug: Boolean = false): AlignResult = {
    val readLen = read1.length
    var ref = new Array[Byte](readLen)
    var rc1 = new Array[Byte](readLen)
    var rc2 = new Array[Byte](readLen)

    DNA.getReverseComplement(read1, rc1)
    DNA.getReverseComplement(read2, rc2)

    val score = new Score(params.maxDist, params.confDiff)
    val lv = new LandauVishkin(params.getUpperBound)
    
    var pos1 = 0L
    var pos2 = 0L
    var read2UpperBound = 0
    
    var distToCheck = 0
    var dist1 = 0
    var dist2 = 0
    
    for (p <- GenomeLoader.genome.pieces) {
      println("Aligning against " + p.name + "...")

      if (debug)
      	 println(p.startIndex)

      // check forward => use read1, rc2
      pos1 = p.startIndex
      while (pos1 < p.endIndex) {
        if (debug && pos1 % 10000000 == 0)
          println("Pos: " + pos1)

        // get reference
        GenomeLoader.genome.getSubstring(pos1, pos1 + readLen, ref)
      
        // get current dist to check (based on best score so far)
        distToCheck = score.getDistToCheck(snapUpperBound)
        
        // get dist for read1
        dist1 = lv.distance(ref, readLen, read1, readLen, distToCheck)
      
        if (dist1 != -1 && dist1 <= params.getUpperBound) {
          // found a pretty good hit for read1; now check around it for rc2

          pos2 = params.getRead2StartPos(p, pos1)
          read2UpperBound = distToCheck - dist1
          
          while (pos2 < params.getRead2EndPos(p, pos1)) {
            GenomeLoader.genome.getSubstring(pos2, pos2 + readLen, ref)
            dist2 = lv.distance(ref, readLen, rc2, readLen, read2UpperBound)
            if (dist2 != -1)  
              score.update(dist1 + dist2, pos1, pos2, false, debug)
            
            pos2 += 1
          }
        }
        
        pos1 += 1
      }
      
      // check reverse => use rc1, read2

      pos1 = p.startIndex
      while (pos1 < p.endIndex) {
        if (debug && pos1 % 10000000 == 0)
          println("Pos: " + pos1)
        
        // get reference
        GenomeLoader.genome.getSubstring(pos1, pos1 + readLen, ref)
        
        // get current dist to check (based on best score so far)
        distToCheck = score.getDistToCheck(snapUpperBound)

        // get dist for rc1
        dist1 = lv.distance(ref, readLen, rc1, readLen, distToCheck)
        
        if (dist1 != -1 && dist1 <= params.getUpperBound) {
          // found a pretty good hit for rc1; now check around it for read2
          
          pos2 = params.getRead2StartPos(p, pos1)
          read2UpperBound = distToCheck - dist1
          
          while (pos2 < params.getRead2EndPos(p, pos1)) {
            GenomeLoader.genome.getSubstring(pos2, pos2 + readLen, ref)
            dist2 = lv.distance(ref, readLen, read2, readLen, read2UpperBound)
	    
	    if (dist2 != -1) {
              score.update(dist1 + dist2, pos1, pos2, true, debug)
            }
            
            pos2 += 1
          }
        }
        
        pos1 += 1
      }
    }

    score.getResult
  }  
}