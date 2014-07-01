package cs.ucla.edu.bwaspark.worker2

import scala.math.sqrt

import cs.ucla.edu.bwaspark.datatype._

object MemSamPe {
  private val MIN_RATIO = 0.8
  private val MIN_DIR_CNT = 10
  private val MIN_DIR_RATIO = 0.05
  private val OUTLIER_BOUND = 2.0
  private val MAPPING_BOUND = 3.0
  private val MAX_STDDEV = 4.0

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

  // Not called in the current setting; Needs to be tested
  def memPeStat(opt: MemOptType, pacLen: Long, n: Int, regArray: Array[MemAlnRegArrayType], pes: Array[MemPeStat]) {
    var iSize: Array[Vector[Int]] = new Array[Vector[Int]](4)

    var i = 0
    while(i < (n>>>1)) {
      var r: Array[MemAlnRegArrayType] = new Array[MemAlnRegArrayType](2)
      var dir: Int = 0
      var is: Long = 0
      r(0) = regArray(i<<1|0)
      r(1) = regArray(i<<1|1)

      if(r(0).regs.size != 0 && r(1).regs.size != 0) {
        if(calSub(opt, r(0).regs) <= MIN_RATIO * r(0).regs(0).score) {
          if(calSub(opt, r(1).regs) <= MIN_RATIO * r(1).regs(0).score) {
            // Inline mem_infer_dir
            var r1: Boolean = false
            var r2: Boolean = false
            if(r(0).regs(0).rBeg >= pacLen) r1 = true
            if(r(1).regs(0).rBeg >= pacLen) r2 = true

            var rBegLarger = r(1).regs(0).rBeg
            // rBegLarger is the coordinate of read 2 on the read 1 strand
            if(r1 != r2) rBegLarger = (pacLen << 1) - 1 - r(1).regs(0).rBeg 
            var dist: Int = (r(0).regs(0).rBeg - rBegLarger).toInt
            if(rBegLarger > r(0).regs(0).rBeg) dist = (rBegLarger - r(0).regs(0).rBeg).toInt
            
            var cond1 = 1
            if(r1 == r2) cond1 = 0
            var cond2 = 3
            if(rBegLarger > r(0).regs(0).rBeg) cond2 = 0
            
            dir = cond1 ^ cond2
            if(dist > 0 && dist <= opt.maxIns) iSize(dir) = iSize(dir) :+ dist
          }
        }
      }

      i += 1
    }
   
    var d = 0
    while(d < 4) {
      val qInit: Vector[Int] = iSize(d)
      if(qInit.size < MIN_DIR_CNT) {
        println("skip orientation as there are not enough pairs")
        pes(d).failed = 1
        d += 1
      } 
      else {
        println("analyzing insert size distribution for orientation")
        val q = qInit.sortWith(_.compareTo(_) < 0) // ks_introsort_64
        var p25: Int = q((0.25 * iSize(d).size + 0.499).toInt)
        var p50: Int = q((0.50 * iSize(d).size + 0.499).toInt)
        var p75: Int = q((0.75 * iSize(d).size + 0.499).toInt)
        pes(d).low = (p25 - OUTLIER_BOUND * (p75 - p25) + 0.499).toInt
        if(pes(d).low < 1) pes(d).low = 1
        pes(d).high = (p75 - OUTLIER_BOUND * (p75 - p25) + 0.499).toInt
        println("(25, 50, 75) percentile: (" + p25 + ", " + p50 + ", " + p75 + ")")
        println("low and high boundaries for computing mean and std.dev: (" + pes(d).low + ", " + pes(d).high + ")")

        i = 0
        var x = 0
        pes(d).avg = 0
        while(i < q.size) {
          if(q(i) >= pes(d).low && q(i) <= pes(d).high) {
            pes(d).avg += q(i)
            x += 1
          }
          i += 1
        }
        pes(d).avg /= x
        
        i = 0
        pes(d).std = 0
        while(i < q.size) {
          if(q(i) >= pes(d).low && q(i) <= pes(d).high)
            pes(d).std += (q(i) - pes(d).avg) * (q(i) - pes(d).avg)
        }
        pes(d).std = sqrt(pes(d).std / x)
        println("mean and std.dev: (" + pes(d).avg + ", " + pes(d).std + ")")

        pes(d).low = (p25 - MAPPING_BOUND * (p75 - p25) + .499).toInt
        pes(d).high = (p75 + MAPPING_BOUND * (p75 - p25) + .499).toInt
        if(pes(d).low > pes(d).avg - MAX_STDDEV * pes(d).std) pes(d).low = (pes(d).avg - MAX_STDDEV * pes(d).std + .499).toInt
        if(pes(d).high < pes(d).avg - MAX_STDDEV * pes(d).std) pes(d).high = (pes(d).avg - MAX_STDDEV * pes(d).std + .499).toInt
        if(pes(d).low < 1) pes(d).low = 1
        println("low and high boundaries for proper pairs: (" + pes(d).low + ", " + pes(d).high + ")")

        d += 1
      }
    } 

    d = 0
    var max = 0
    while(d < 4) {
      if(max < iSize(d).size) max = iSize(d).size
      d += 1
    }

    d = 0
    while(d < 4) {
      if(pes(d).failed == 0 && iSize(d).size < max * MIN_DIR_RATIO) {
        pes(d).failed = 1
        println("skip orientation")
      }     
    }
  }
/*
  private def memMateSw(opt: MemOptType, pacLen: Long, pac: Array[Byte], pes: Array[MemPeStat], reg: MemAlnRegType, 
                        mateSeqLen: Int, mateSeq: Array[Byte], mateRegs: Array[MemAlnRegType]): Int = {
       
  }
*/
}

