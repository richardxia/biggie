package snap

/**
 * Represents a two-allele combination at a given location, which may either match the
 * reference or be a heterozygous or homozygous SNP.
 */
case class SNP(val reference: Char, val allele1: Char, val allele2: Char) {
  override def toString = {
    if (allele1 == reference && allele2 == reference) {
      reference + "/" + reference
    } else if (allele1 == reference && allele2 != reference) {
      reference + "/" + allele2 + "*"
    } else if (allele1 != reference && allele2 == reference) {
      reference + "/" + allele1 + "*"
    } else {  // Both not equal to reference
      allele1 + "/" + allele2 + "!"
    }
  }

  def diffString: String = {
    if (reference == '?')
      "?"
    else if (reference == allele1 && reference == allele2)
      "-"
    else
      allele1 + "/" + allele2
  }

  def contains(allele: Char): Boolean = allele1 == allele || allele2 == allele

  def uncalled: Boolean = reference == '?'
}

object SNP {
  val UNCALLED = SNP('?', '?', '?')

  // Use only one copy of the "unchanged from reference" SNPs to save memory
  val MATCHING = Map(
    'A' -> SNP('A', 'A', 'A'),
    'C' -> SNP('C', 'C', 'C'),
    'G' -> SNP('G', 'G', 'G'),
    'T' -> SNP('T', 'T', 'T')
  )
}
