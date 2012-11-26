package snap

object GenomeLoader {
  lazy val genome = FASTA.read("/disk/1/kcurtis/simFinder/hg19.fa")
  //lazy val genome = FASTA.read("/home/eecs/kcurtis/chr22.fa")
}