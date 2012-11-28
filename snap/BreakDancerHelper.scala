package snap

class BreakDancerHelper(bdLine: String) {
  val fields = bdLine.split(" +")
  
  def chromosome = fields(0)
  
  def startPos = fields(1).toLong
  
  def endPos = fields(4).toLong  // assumes breakpoints are on same chromosome
  
  def size = fields(7).toInt
  
  def svType = fields(6)
  
  def getVcfLine = {
    svType match {
      // dummy format
      case "DEL" => List(chromosome, startPos, ".", Array.fill(size + 1)("n").mkString(""), "n", "200", "PASS", ".", "GT", "0|0").mkString("\t")
      case "INS" => List(chromosome, startPos, ".", "n", Array.fill(size + 1)("n").mkString(""), "200", "PASS", ".", "GT", "0|0").mkString("\t")
      case _ => ""
    }
  }
}

object BreakDancerClient {
  def toVcf(fname: String) = {
    val source = scala.io.Source.fromFile(fname)
    val lines = source.getLines.filter(l => !l.startsWith("#")).toList
    
    var helper:BreakDancerHelper = null
    
    lines.foreach(l => {
      helper = new BreakDancerHelper(l)
      println(helper.getVcfLine)
    })
  }
}