package biggie

import java.io.File
import net.sf.samtools.SAMRecordIterator
import net.sf.samtools.SAMRecord
import net.sf.samtools.SAMFileReader

class SamRegionReader(samFile: String, var refSeq: String, var region: Range) {
  val indexFile = samFile + ".bai"
  val file = new SAMFileReader(new File(samFile), new File(indexFile))

  // Don't validate SAM file
  file.setValidationStringency(SAMFileReader.ValidationStringency.SILENT)

  def reads(): SAMRecordIterator = {
    file.queryOverlapping(refSeq, region.start, region.end)
  }
}
