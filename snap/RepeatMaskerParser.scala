package snap

class RepeatMaskerParser(repeatLine: String) {
  val fields = repeatLine.split(" +")
  
  def querySequence = fields(5)
  
  def beginInQuery = fields(6).toLong
  
  def endInQuery = fields(7).toLong
  
  def repeatLength = endInQuery - beginInQuery + 1
  
  def isRC = fields(9) == "C"
  
  def repeatDescriptor = fields(10)
  
  def repeatClass = fields(11)
}