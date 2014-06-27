package cs.ucla.edu.bwaspark.worker2

import cs.ucla.edu.bwaspark.datatype._

object MemSamPe {
  private def calSub(opt: MemOptType, regs: Array[MemAlnRegType]) : Int = {
    var j = 1
    var isBreak = false
    while(j < regs.length && !isBreak) { // choose unique alignment
      var bMax = regs(0).qBeg
      var eMin = regs(0).qEnd

      if(regs(j).qBeg > regs(0).qBeg) bMax = regs(j).qBeg
      if(regs(j).qEnd < regs(0).qEnd) eMin = regs(j).qEnd
      if(eMin > bMax) { // have overlap
        var minL = regs(0).qEnd - regs(0).qBeg
	if(regs(j).qEnd - regs(j).qBeg < minL) minL = regs(j).qEnd - regs(j).qBeg
	if(eMin - bMax >= minL * opt.maskLevel) { 
          isBreak = true
	  j -= 1
        }
      }

      j += 1    
    }

    if(j < regs.length) regs(j).score
    else opt.minSeedLen * opt.a
  }

  def memPeStat(opt: MemOptType, pacLen: Long, n: Int, regs: Array[MemAlnRegType], pes: Array[MemPeStat]) {
    var iSize: Array[Long] = new Array[Long](4)

    var i = 0
    while(i < (n>>>1)) {
      var r: Array[MemAlnRegType] = new Array[MemAlnRegType](2)
      r(0) = regs(i<<1|0)
      r(1) = regs(i<<1|1)

      //if(r(0).n != 0 && r(1).n != 0) {
        //if(calSub(opt, r(0)) 
      //}

      i += 1
    }
    
  }
}

