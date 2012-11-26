package snap

import snap._

sealed trait Point extends Serializable {
  def distance(that: Point): Int
}

case class SimplePoint(val pos: Int) extends Point {
  override def distance(that: Point): Int = {
    scala.math.abs(this.pos - that.asInstanceOf[SimplePoint].pos)
  }

  override def toString = pos.toString
}

case class ReadPoint(val read: String) extends Point {
  override def distance(that: Point): Int = {
    val l = new Levenshtein
    l.distance(this.read, that.asInstanceOf[ReadPoint].read)
  }
  
  override def toString = read
}
