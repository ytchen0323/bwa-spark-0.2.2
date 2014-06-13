package cs.ucla.edu.bwaspark.worker2

import cs.ucla.edu.bwaspark.datatype._
import cs.ucla.edu.bwaspark.worker2.MemMarkPrimarySe._
import cs.ucla.edu.bwaspark.worker2.MemRegToADAMSAM._

object BWAMemWorker2 {
  /**
    *  Main function of BWA-mem worker2
    *
    *  @param opt the input MemOptType object
    *  @param regs the alignment registers to be transformed
    *  @param bns the input BNSSeqType object
    *  @param pac the PAC array
    *  @param seq the read (NOTE: currently we use Array[Byte] first. may need to be changed!!!)
    */
  def bwaMemWorker2(opt: MemOptType, regs: Array[MemAlnRegType], bns: BNTSeqType, pac: Array[Byte], seq: FASTQSingleNode, numProcessed: Long) {
    var regsOut: Array[MemAlnRegType] = null
    if(regs != null)
      regsOut = memMarkPrimarySe(opt, regs, numProcessed)

    memRegToSAMSe(opt, bns, pac, seq, regsOut, 0, null)
  }

}

