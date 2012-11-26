package snap

import scala.io.Source

import scala.math.{max, min}

object TVSim {
  // Read a VCF file produced by tvsim and return a list of SNPs as (chrom, pos, SNP) tuples
  def readSNPs(trueVCF: String): Array[(String, Int, SNP)] = {
    val entries = Source.fromFile(trueVCF).getLines
                        .filter(!_.startsWith("#"))
                        .map(line => Utils.split(line, '\t'))
    // Keep only the entries where the REF and ALT fields are one base, which represent SNPs
    val snps = entries.filter(fields => isOneBase(fields(3)) && isOneBase(fields(4)))
    snps.map { fields =>
      val chrom = new String(fields(0))
      val pos = fields(1).toInt
      val ref = fields(3).charAt(0)
      val alts = fields(4).toUpperCase
      val alt1 = alts.charAt(0)
      val alt2 = if (alts.length == 1) alts.charAt(0) else alts.charAt(2)
      val snp = if (fields(9) == "1|1") {
        SNP(ref, alt1, alt1)
      } else if (fields(9) == "0|1" || fields(9) == "1|0") {
        SNP(ref, ref, alt1)
      } else {
        SNP(ref, min(alt1, alt2).toChar, max(alt1, alt2).toChar)
      }
      (chrom, pos, snp)
    }.toArray
  }

  def isOneBase(s: String): Boolean = s.length == 1 || s.length == 3 && s(1) == ','
}
