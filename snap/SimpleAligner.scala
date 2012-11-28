package snap
import java.util.Random
import it.unimi.dsi.fastutil.longs.LongArrayList
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import scala.math.min
import scala.collection.mutable.Set
import scala.collection.JavaConversions._

class SimpleAligner(genome: Genome, index: Index, seedsToTry: Int, maxDist: Int, maxHits: Int = 40) {
  private val seedLen = index.seedLen
  private val rand = new Random(42)
  private val lv = new LandauVishkin(maxDist)

  // given a read, return a list of all the positions to which it aligns
  def align(read: String): Set[Long] = {
    val hits = Set[Long]()
    var offset = 0
    (1 to seedsToTry).foreach(seedNum => {
      val fwdHits = new LongArrayList
      
      index.get(DNA.substringToLong(read, offset, offset + seedLen), fwdHits, maxHits)

      /*
      // debugging
      if (seedNum < 100) {
        println("seed" + seedNum + ", " + read.substring(offset, offset + seedLen) + ", has " + fwdHits.size + " hits")
        fwdHits.foreach(println)
      }
      */
            
      fwdHits.foreach(hit => hits += hit)

      offset = ((seedNum / (seedsToTry - 1.0)) * (read.length - seedLen)).round.toInt
    })
    
    // filter out any hits that are too far away
    hits.filter(hit => lv.distance(read, genome.substring(hit, hit + read.length - 1), maxDist) != -1)
    
    hits
  }
}
