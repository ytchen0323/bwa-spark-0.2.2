package cs.ucla.edu.bwaspark.worker2

import org.apache.commons.math3.special.Erf.erfc

import scala.math.sqrt
import scala.math.log
import scala.math.abs
import scala.collection.immutable.Vector

import cs.ucla.edu.bwaspark.datatype._
import cs.ucla.edu.bwaspark.worker1.MemSortAndDedup.memSortAndDedup
import cs.ucla.edu.bwaspark.util.BNTSeqUtil.bnsGetSeq
import cs.ucla.edu.bwaspark.util.SWUtil.SWAlign2
import cs.ucla.edu.bwaspark.worker2.MemMarkPrimarySe.{hash64, memMarkPrimarySe}
import cs.ucla.edu.bwaspark.worker2.MemRegToADAMSAM.{memApproxMapqSe, memRegToAln, memAlnToSAM, memRegToSAMSe}

import java.io.FileReader
import java.io.BufferedReader

// testing use
import java.math.BigInteger

object MemSamPe {
  private val MIN_RATIO = 0.8
  private val MIN_DIR_CNT = 10
  private val MIN_DIR_RATIO = 0.05
  private val OUTLIER_BOUND = 2.0
  private val MAPPING_BOUND = 3.0
  private val MAX_STDDEV = 4.0
  private val KSW_XBYTE = 0x10000
  private val KSW_XSTOP = 0x20000
  private val KSW_XSUBO = 0x40000
  private val KSW_XSTART = 0x80000
  private val M_SQRT1_2 = 7.0710678118654752440E-1  
  private val MEM_F_NO_RESCUE = 0x20
  private val MEM_F_NOPAIRING = 0x4

  //pre-process: transform A/C/G/T to 0,1,2,3
  private def locusEncode(locus: Char): Byte = {
    //transforming from A/C/G/T to 0,1,2,3
    locus match {
      case 'A' => 0
      case 'a' => 0
      case 'C' => 1
      case 'c' => 1
      case 'G' => 2
      case 'g' => 2
      case 'T' => 3
      case 't' => 3
      case '-' => 5
      case _ => 4
    }
  }

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
    var j = 0
    while(j < 4) {
      iSize(j) = scala.collection.immutable.Vector.empty
      j += 1
    }

    var i = 0
    while(i < (n>>>1)) {
      var r: Array[MemAlnRegArrayType] = new Array[MemAlnRegArrayType](2)
      var dir: Int = 0
      var is: Long = 0
      r(0) = regArray(i<<1|0)
      r(1) = regArray(i<<1|1)

      if(r(0) != null && r(1) != null) {
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
        pes(d).high = (p75 + OUTLIER_BOUND * (p75 - p25) + 0.499).toInt
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

          i += 1
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

      d += 1
    }
  }

  // verified
  private def memMateSw(opt: MemOptType, pacLen: Long, pac: Array[Byte], pes: Array[MemPeStat], reg: MemAlnRegType, 
                        mateSeqLen: Int, mateSeq: Array[Byte], mateRegs: Array[MemAlnRegType]): (Int, Array[MemAlnRegType]) = {
    var mateRegsUpdated: Vector[MemAlnRegType] = scala.collection.immutable.Vector.empty
    var skip: Array[Int] = new Array[Int](4)
    var n = 0
    var regArray = new MemAlnRegArrayType

    var r = 0
    while(r < 4) {
      if(pes(r).failed > 0) skip(r) = 1
      else skip(r) = 0
      r += 1     
    }
    //println("[0] " + skip(0) + " " + skip(1) + " " + skip(2) + " " + skip(3))

    var i = 0
    if(mateRegs != null) {
      while(i < mateRegs.size) { // consistent pair exist; no need to perform SW
        // Inline mem_infer_dir
        var r1: Boolean = false
        var r2: Boolean = false
        if(reg.rBeg >= pacLen) r1 = true
        if(mateRegs(i).rBeg >= pacLen) r2 = true

        var rBegLarger = mateRegs(i).rBeg
        // rBegLarger is the coordinate of read 2 on the read 1 strand
        if(r1 != r2) rBegLarger = (pacLen << 1) - 1 - mateRegs(i).rBeg
        var dist: Int = (reg.rBeg - rBegLarger).toInt
        if(rBegLarger > reg.rBeg) dist = (rBegLarger - reg.rBeg).toInt

        var cond1 = 1
        if(r1 == r2) cond1 = 0
        var cond2 = 3
        if(rBegLarger > reg.rBeg) cond2 = 0

        r = cond1 ^ cond2      
        if(dist >= pes(r).low && dist <= pes(r).high) skip(r) = 1     

        i += 1
      }
    }
 
    //println("[1] " + skip(0) + " " + skip(1) + " " + skip(2) + " " + skip(3))
    if(skip(0) + skip(1) + skip(2) + skip(3) == 4) (0, mateRegs)   // consistent pair exist; no need to perform SW

    if(mateRegs != null) {
      i = 0
      while(i < mateRegs.size) {
        mateRegsUpdated = mateRegsUpdated :+ mateRegs(i)
        i += 1
      }
    }

    r = 0
    while(r < 4) {
      if(skip(r) == 0) {
        var seq: Array[Byte] = mateSeq
        var rBeg: Long = -1
        var rEnd: Long = -1
        var len: Long = 0

        var isRev = 0
        if((r >> 1) != (r & 1)) isRev = 1   // whether to reverse complement the mate

        var isLarger = 1
        if((r >> 1) > 0) isLarger = 0   // whether the mate has larger coordinate

        if(isRev > 0) {
          var rev: Array[Byte] = new Array[Byte](mateSeqLen)
          i = 0
          while(i < mateSeqLen) {
            if(mateSeq(i) < 4) rev(mateSeqLen - 1 - i) = (3 - mateSeq(i)).toByte
            else rev(mateSeqLen - 1 - i) = 4
            i += 1
          }
          seq = rev
        }
        
        if(isRev == 0) {
          if(isLarger > 0) rBeg = reg.rBeg + pes(r).low
          else rBeg = reg.rBeg - pes(r).high

          if(isLarger > 0) rEnd = reg.rBeg + pes(r).high + mateSeqLen // if on the same strand, end position should be larger to make room for the seq length
          else rEnd = reg.rBeg - pes(r).low + mateSeqLen
        }
        else {
          if(isLarger > 0) rBeg = reg.rBeg + pes(r).low - mateSeqLen // similarly on opposite strands
          else rBeg = reg.rBeg - pes(r).high - mateSeqLen

          if(isLarger > 0) rEnd = reg.rBeg + pes(r).high
          else rEnd = reg.rBeg - pes(r).low
        }

        if(rBeg < 0) rBeg = 0
        if(rEnd > (pacLen << 1)) rEnd = pacLen << 1

        val ref = bnsGetSeq(pacLen, pac, rBeg, rEnd) 
        if(ref._2 == rEnd - rBeg) { // no funny things happening  ( ref._2 -> len)
          var xtraTmp = 0
          if(mateSeqLen * opt.a < 250) xtraTmp = KSW_XBYTE
          val xtra = KSW_XSUBO | KSW_XSTART | xtraTmp | (opt.minSeedLen * opt.a)
          val aln = SWAlign2(mateSeqLen, seq, ref._2.toInt, ref._1, 5, opt, xtra) // ref._1 -> ref
          
          //println("aln score: " + aln.score + ", xtra: " + xtra)  
          var alnTmp = new MemAlnRegType
          if(aln.score >= opt.minSeedLen && aln.qBeg >= 0) { // something goes wrong if aln.qBeg < 0
            if(isRev > 0) {
              alnTmp.qBeg = mateSeqLen - (aln.qEnd + 1)
              alnTmp.qEnd = mateSeqLen - aln.qBeg
              alnTmp.rBeg = (pacLen << 1) - (rBeg + aln.tEnd + 1)
              alnTmp.rEnd = (pacLen << 1) - (rBeg + aln.tBeg)
            }
            else {
              alnTmp.qBeg = aln.qBeg
              alnTmp.qEnd = aln.qEnd + 1
              alnTmp.rBeg = rBeg + aln.tEnd + 1
              alnTmp.rEnd = rBeg + aln.tEnd + 1
            }

            alnTmp.score = aln.score
            alnTmp.csub = aln.scoreSecond
            alnTmp.secondary = -1

            if(alnTmp.rEnd - alnTmp.rBeg < alnTmp.qEnd - alnTmp.qBeg) alnTmp.seedCov = ((alnTmp.rEnd - alnTmp.rBeg) >>> 1).toInt
            else alnTmp.seedCov = (alnTmp.qEnd - alnTmp.qBeg) >>> 1
            
            // move b s.t. ma is sorted
            mateRegsUpdated = mateRegsUpdated :+ alnTmp

            // original sorting in c codes
            /*
            i = 0
            var isBreak = false
            while(i < mateRegs.size && isBreak) {
              if(mateRegs(i).score < alnTmp.score) isBreak = true
              else i += 1
            }
            val breakIdx = i

            i = mateRegs.size - 1
            while(i > breakIdx) {
              mateRegs(i) = mateRegs(i - 1)
              i -= 1
            }
            mateRegs(breakIdx) = alnTmp
            */
          }
          
          n += 1
        }

        if(n > 0) {
          // Sort here!
          //println("size: " + mateRegsUpdated.size)
          //mateRegsUpdated.foreach(ele => println(ele.rBeg + " " + ele.rEnd + " " + ele.qBeg + " " + ele.qEnd + " " + ele.score + " " + ele.trueScore + " " + ele.sub + " " + ele.csub + " " + ele.subNum + " " + ele.width + " " + ele.seedCov + " " + ele.secondary + " " + ele.hash))
          mateRegsUpdated = mateRegsUpdated.sortBy(seq => (seq.score))
          regArray.regs = mateRegsUpdated.toArray
          regArray.curLength = mateRegsUpdated.size
          regArray.maxLength = mateRegsUpdated.size
          regArray = memSortAndDedup(regArray, opt.maskLevelRedun)
          
        }
      }

      r += 1
    }

    if(n > 0) (n, regArray.regs)
    else (n, mateRegs)
  }

  /**
    *  Data structure which keeps a pair of Long integer
    */
  class PairLong {
    var x: Long = 0
    var y: Long = 0
  } 

  /**
    *  Data structure which keeps a pair of Long integer
    */
  class QuadLong {
    var xMSB: Long = 0
    var xRest: Long = 0
    var x: Long = 0
    var y: Long = 0
  } 


  // Verification done with the c version!!!
  private def memPair(opt: MemOptType, pacLen: Long, pac: Array[Byte], pes: Array[MemPeStat], alnRegVec: Array[Array[MemAlnRegType]], id: Long): (Int, Int, Int, Array[Int]) = {
    var keyVec: Vector[PairLong] = scala.collection.immutable.Vector.empty
    //var keyUVec: Vector[QuadLong] = scala.collection.immutable.Vector.empty
    var keyUVec: Vector[PairLong] = scala.collection.immutable.Vector.empty
    var y: Array[Int] = new Array[Int](4)
    var r = 0
    var i = 0
    var z: Array[Int] = new Array[Int](2)

    z(0) = -1
    z(1) = -1

    while(r < 2) {
      i = 0
      while(i < alnRegVec(r).size) {
        var key: PairLong = new PairLong

        key.x = alnRegVec(r)(i).rBeg
        if(alnRegVec(r)(i).rBeg >= pacLen) key.x = (pacLen << 1) - 1 - alnRegVec(r)(i).rBeg   // forward position
        if(alnRegVec(r)(i).rBeg >= pacLen) key.y = (alnRegVec(r)(i).score.toLong << 32) | (i << 2) | 2 | r
        else key.y = (alnRegVec(r)(i).score.toLong << 32) | (i << 2) | 0 | r

        //println("x: " + key.x + ", y: " + key.y)

        keyVec = keyVec :+ key
        i += 1
      }

      r += 1
    }

    val sortedKeyVec = keyVec.sortBy(key => (key.x, key.y))
    y(0) = -1 
    y(1) = -1
    y(2) = -1
    y(3) = -1

    i = 0
    while(i < sortedKeyVec.size) {
      r = 0
      while(r < 2) { // loop through direction
        var dir = ((r << 1) | (sortedKeyVec(i).y >>> 1 & 1)).toInt
        var which: Int = 0
        if(pes(dir).failed > 0) { // invalid orientation
          r += 1
        }
        else {
          which = ((r << 1) | ((sortedKeyVec(i).y & 1) ^ 1)).toInt
          if(y(which) < 0) { // no previous hits
            r += 1
          }
          else {
            var k = y(which)
            var isBreak = false
            while(k >= 0 && !isBreak) { // TODO: this is a O(n^2) solution in the worst case; remember to check if this loop takes a lot of time (I doubt)
              if((sortedKeyVec(k).y & 3) != which) 
                k -= 1
              else {
                var dist: Long = sortedKeyVec(i).x - sortedKeyVec(k).x
                //println(k + ": " + dist)
                if(dist > pes(dir).high) isBreak = true
                else if(dist < pes(dir).low) k -= 1
                else {
                  var ns: Double = (dist - pes(dir).avg) / pes(dir).std
                  var q: Int = ((sortedKeyVec(i).y >>> 32) + (sortedKeyVec(k).y >>> 32) + 0.721 * log(2.0 * erfc(abs(ns) * M_SQRT1_2)) * opt.a + 0.499).toInt
                  if(q < 0) q = 0
                  
                  //var key: QuadLong = new QuadLong
                  var key: PairLong = new PairLong
                  key.y = k.toLong << 32 | i
                  key.x = q.toLong << 32 | ((hash64(key.y ^ id << 8) << 32) >>> 32) // modify from (hash64(key.y ^ id << 8) & 0xffffffff)
                  //key.xMSB = key.x >>> 63
                  //key.xRest = (key.x << 1) >>> 1
                  //println("[" + sortedKeyVec(k).x + "," + sortedKeyVec(i).x + "]\t" + q + "\tdist=" + dist)
                  //val tmp = (hash64(key.y ^ id << 8) << 32) >>> 32
                  //println("[0] " + key.y + " " + key.x + " " + key.xMSB + " " + key.xRest + " " + tmp);
                  keyUVec = keyUVec :+ key
                  k -= 1
                }
              }
            }
            
            r += 1
          }
        }
      }      

      y((sortedKeyVec(i).y & 3).toInt) = i
      i += 1
    }

    //println("keyU size: " + keyUVec.size)
    //keyUVec.foreach(key => println(key.y + " " + key.x + " " + key.xMSB + " " + key.xRest))
    var ret = 0
    var sub = 0
    var numSub = 0
    if(keyUVec.size > 0) { // found at least one proper pair
      var tmp: Int = opt.a + opt.b
      if(tmp < opt.oDel + opt.eDel) tmp = opt.oDel + opt.eDel
      if(tmp < opt.oIns + opt.eIns) tmp = opt.oIns + opt.eIns
      
      //val sortedKeyUVec = keyUVec.sortBy(key => (key.xMSB, key.xRest, key.y))
      val sortedKeyUVec = keyUVec.sortBy(key => (key.x, key.y))
      var i = (sortedKeyUVec(sortedKeyUVec.size - 1).y >>> 32).toInt
      var k = (sortedKeyUVec(sortedKeyUVec.size - 1).y << 32 >>> 32).toInt
      //println("i: " + i + ", k: " + k)
      z((sortedKeyVec(i).y & 1).toInt) = (sortedKeyVec(i).y << 32 >>> 34).toInt
      z((sortedKeyVec(k).y & 1).toInt) = (sortedKeyVec(k).y << 32 >>> 34).toInt
      ret = (sortedKeyUVec(sortedKeyUVec.size - 1).x >>> 32).toInt
      if(sortedKeyUVec.size > 1) sub = (sortedKeyUVec(sortedKeyUVec.size - 2).x >>> 32).toInt

      //println("rest: " + (sortedKeyUVec(sortedKeyUVec.size - 1).xRest >>> 32) + " " + (sortedKeyUVec(sortedKeyUVec.size - 2).xRest >>> 32))

      i = sortedKeyUVec.size - 2
      while(i >= 0) {
        if(sub - (sortedKeyUVec(i).x >>> 32).toInt <= tmp) numSub += 1
        i -= 1
      }
    } 
    else {
      ret = 0
      sub = 0
      numSub = 0
    }

    //println(ret + " " + sub + " " + numSub + " " + z(0) + " " + z(1))
    (ret, sub, numSub, z)
  }

  def memSamPe(opt: MemOptType, bns: BNTSeqType, pac: Array[Byte], pes: Array[MemPeStat], id: Long, 
               seqs: Array[FASTQSingleNode], alnRegVec: Array[Array[MemAlnRegType]]): Int = {
    var n: Int = 0
    var z: Array[Int] = new Array[Int](2)
    var subo: Int = 0
    var numSub: Int = 0
    var extraFlag: Int = 1
    var seqsTrans: Array[Array[Byte]] = new Array[Array[Byte]](2)
    var noPairingFlag: Boolean = false

    // for testMemSamPe use
    //seqsTrans(0) = seqs(0).seq.toCharArray.map(ele => (ele - 48).toByte)
    //seqsTrans(1) = seqs(1).seq.toCharArray.map(ele => (ele - 48).toByte)

    seqsTrans(0) = seqs(0).seq.toCharArray.map(ele => locusEncode(ele))
    seqsTrans(1) = seqs(1).seq.toCharArray.map(ele => locusEncode(ele))

    //print("seq0: ")
    //seqsTrans(0).foreach(ele => print(ele))
    //println
    //print("seq1: ")
    //seqsTrans(1).foreach(ele => print(ele))
    //println

    if((opt.flag & MEM_F_NO_RESCUE) == 0) { // then perform SW for the best alignment
      var i = 0
      var alnRegTmpVec = new Array[Vector[MemAlnRegType]](2)
      alnRegTmpVec(0) = scala.collection.immutable.Vector.empty
      alnRegTmpVec(1) = scala.collection.immutable.Vector.empty

      while(i < 2) {
        if(alnRegVec(i) != null) {
          var j = 0
          while(j < alnRegVec(i).size) {
            if(alnRegVec(i)(j).score >= alnRegVec(i)(0).score - opt.penUnpaired)
              alnRegTmpVec(i) = alnRegTmpVec(i) :+ alnRegVec(i)(j)
            j += 1
          }
        }

        i += 1
      }

      //println("Len: " + alnRegTmpVec(0).size + ", " + alnRegTmpVec(1).size)

      i = 0
      while(i < 2) {
        var j = 0
        while(j < alnRegTmpVec(i).size && j < opt.maxMatesw) {
          var iBar = 0
          if(i == 0) iBar = 1
          val ret = memMateSw(opt, bns.l_pac, pac, pes, alnRegTmpVec(i)(j), seqs(iBar).seqLen, seqsTrans(iBar), alnRegVec(iBar)) 
          n += ret._1
          alnRegVec(iBar) = ret._2
          j += 1
        }

        i += 1
      }  
    } 
    
    alnRegVec(0) = memMarkPrimarySe(opt, alnRegVec(0), id<<1|0)
    alnRegVec(1) = memMarkPrimarySe(opt, alnRegVec(1), id<<1|1)

    if((opt.flag & MEM_F_NOPAIRING) == 0) {
      // pairing single-end hits
      var a0Size = 0
      var a1Size = 0
      if(alnRegVec(0) != null) a0Size = alnRegVec(0).size
      if(alnRegVec(1) != null) a1Size = alnRegVec(1).size
      //println(a0Size + " " + a1Size)

      if(alnRegVec(0) != null && alnRegVec(1) != null) {
        val retVal = memPair(opt, bns.l_pac, pac, pes, alnRegVec, id)
        val ret = retVal._1
        subo = retVal._2
        numSub = retVal._3
        z = retVal._4

        //println(subo + " " + numSub + " " + ret + " " + z(0) + " " + z(1))

        if(ret > 0) {
          var scoreUn: Int = 0
          var isMulti: Array[Boolean] = new Array[Boolean](2)
          var qPe: Int = 0
          var qSe: Array[Int] = new Array[Int](2)

          var i = 0
          while(i < 2) {
            var j = 1
            var isBreak = false
            while(j < alnRegVec(i).size && !isBreak) {
              if(alnRegVec(i)(j).secondary < 0 && alnRegVec(i)(j).score >= opt.T) 
                isBreak = true
              else 
                j += 1
            }
              
            if(j < alnRegVec(i).size) isMulti(i) = true
            else isMulti(i) = false

            //println("j " + j + ", alnRegVec(i).size: " + alnRegVec(i).size)
            i += 1
          }

          //println(isMulti(0) + " " + isMulti(1))
          if(!isMulti(0) && !isMulti(1)) {  
            // compute mapQ for the best SE hit
            scoreUn = alnRegVec(0)(0).score + alnRegVec(1)(0).score - opt.penUnpaired
            if(subo < scoreUn) subo = scoreUn
            // Inline raw_mapq(ret - subo, opt.a)
            // #define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))
            qPe = (6.02 * (ret - subo) / opt.a + 0.499).toInt
            
            if(numSub > 0) qPe -= (4.343 * log(numSub + 1) + .499).toInt
            if(qPe < 0) qPe = 0
            if(qPe > 60) qPe = 60

            // the following assumes no split hits
            if(ret > scoreUn) { // paired alignment is preferred
              var tmpRegs = new Array[MemAlnRegType](2)
              //println(alnRegVec(0).size + " " + z(0) + " " + alnRegVec(1).size + " " + z(1))
              tmpRegs(0) = alnRegVec(0)(z(0))
              tmpRegs(1) = alnRegVec(1)(z(1))

              var i = 0
              while(i < 2) {
                if(tmpRegs(i).secondary >= 0) {
                  tmpRegs(i).sub = alnRegVec(i)(tmpRegs(i).secondary).score
                  tmpRegs(i).secondary = -1
                }

                qSe(i) = memApproxMapqSe(opt, tmpRegs(i))
                i += 1
              }

              if(qSe(0) < qPe) {
                if(qPe < qSe(0) + 40) qSe(0) = qPe
                else qSe(0) = qSe(0) + 40
              }

              if(qSe(1) < qPe) {
                if(qPe < qSe(1) + 40) qSe(1) = qPe
                else qSe(1) = qSe(1) + 40
              }
              
              extraFlag |= 2
                           
              // cap at the tandem repeat score
              // Inline raw_mapq(tmpRegs(0).score - tmpRegs(0).csub, opt.a)
              var tmp = (6.02 * (tmpRegs(0).score - tmpRegs(0).csub) / opt.a + 0.499).toInt
              if(qSe(0) > tmp) qSe(0) = tmp
              // Inline raw_mapq(tmpRegs(1).score - tmpRegs(1).csub, opt.a)
              tmp = (6.02 * (tmpRegs(1).score - tmpRegs(1).csub) / opt.a + 0.499).toInt
              if(qSe(1) > tmp) qSe(1) = tmp
            }
            else { // the unpaired alignment is preferred
              z(0) = 0
              z(1) = 0
              qSe(0) = memApproxMapqSe(opt, alnRegVec(0)(0))
              qSe(1) = memApproxMapqSe(opt, alnRegVec(1)(0))
            }
            
            // write SAM
            val aln0 = memRegToAln(opt, bns, pac, seqs(0).seqLen, seqsTrans(0), alnRegVec(0)(z(0))) 
            aln0.mapq = qSe(0).toShort
            aln0.flag |= 0x40 | extraFlag
            val aln1 = memRegToAln(opt, bns, pac, seqs(1).seqLen, seqsTrans(1), alnRegVec(1)(z(1)))
            aln1.mapq = qSe(1).toShort
            aln1.flag |= 0x80 | extraFlag

            var samStr0 = new SAMString
            var alnList0 = new Array[MemAlnType](1)
            alnList0(0) = aln0
            memAlnToSAM(bns, seqs(0), seqsTrans(0), alnList0, 0, aln1, samStr0)
            seqs(0).sam = samStr0.str.dropRight(samStr0.size - samStr0.idx).mkString
            var samStr1 = new SAMString
            var alnList1 = new Array[MemAlnType](1)
            alnList1(0) = aln1
            memAlnToSAM(bns, seqs(1), seqsTrans(1), alnList1, 0, aln0, samStr1)
            seqs(1).sam = samStr1.str.dropRight(samStr1.size - samStr1.idx).mkString

            if(seqs(0).name != seqs(1).name) println("[Error] paired reads have different names: " + seqs(0).name + ", " + seqs(1).name)
            
            //println("[0] " + n)
          }
          else // else: goto no_pairing (TODO: in rare cases, the true hit may be long but with low score)
            noPairingFlag = true
        }
        else // else: goto no_pairing
          noPairingFlag = true
      }
      else // else: goto no_pairing 
        noPairingFlag = true 
    }
    else // else: goto no_pairing
      noPairingFlag = true

    if(!noPairingFlag) n
    else { // no_pairing: start from here
      var alnVec: Array[MemAlnType] = new Array[MemAlnType](2)

      var i = 0
      while(i < 2) {
        if(alnRegVec(i) != null && alnRegVec(i).size > 0 && alnRegVec(i)(0).score >= opt.T) alnVec(i) = memRegToAln(opt, bns, pac, seqs(i).seqLen, seqsTrans(i), alnRegVec(i)(0)) 
        else alnVec(i) = memRegToAln(opt, bns, pac, seqs(i).seqLen, seqsTrans(i), null) 

        i += 1
      }

      // if the top hits from the two ends constitute a proper pair, flag it.
      if((opt.flag & MEM_F_NOPAIRING) == 0 && alnVec(0).rid == alnVec(1).rid && alnVec(0).rid >= 0) {
        // Inline mem_infer_dir
        var r1: Boolean = false
        var r2: Boolean = false
        if(alnRegVec(0)(0).rBeg >= bns.l_pac) r1 = true
        if(alnRegVec(1)(0).rBeg >= bns.l_pac) r2 = true

        var rBegLarger = alnRegVec(1)(0).rBeg
        // rBegLarger is the coordinate of read 2 on the read 1 strand
        if(r1 != r2) rBegLarger = (bns.l_pac << 1) - 1 - alnRegVec(1)(0).rBeg
        var dist: Int = (alnRegVec(0)(0).rBeg - rBegLarger).toInt
        if(rBegLarger > alnRegVec(0)(0).rBeg) dist = (rBegLarger - alnRegVec(0)(0).rBeg).toInt

        var cond1 = 1
        if(r1 == r2) cond1 = 0
        var cond2 = 3
        if(rBegLarger > alnRegVec(0)(0).rBeg) cond2 = 0

        var d = cond1 ^ cond2
     
        if(pes(d).failed == 0 && dist >= pes(d).low && dist <= pes(d).high) extraFlag |= 2 
      }

      memRegToSAMSe(opt, bns, pac, seqs(0), alnRegVec(0), 0x41 | extraFlag, alnVec(1))
      memRegToSAMSe(opt, bns, pac, seqs(1), alnRegVec(1), 0x81 | extraFlag, alnVec(0))

      if(seqs(0).name != seqs(1).name) println("[Error] paired reads have different names: " + seqs(0).name + ", " + seqs(1).name)

      //println("[1] " + n)
      n   // return value
    }

  }


  // test functions
  def testMemMateSw(fileName: String, opt: MemOptType, pacLen: Long, pac: Array[Byte], pes: Array[MemPeStat]) {

    val reader = new BufferedReader(new FileReader(fileName))
    var line = reader.readLine
    var i = 0

    while(line != null) {
      val lineFields = line.split(" ")
      var reg = new MemAlnRegType
      reg.rBeg = lineFields(0).toLong
      reg.rEnd = lineFields(1).toLong
      reg.qBeg = lineFields(2).toInt
      reg.qEnd = lineFields(3).toInt
      reg.score = lineFields(4).toInt
      reg.trueScore = lineFields(5).toInt
      reg.sub = lineFields(6).toInt
      reg.csub = lineFields(7).toInt
      reg.subNum = lineFields(8).toInt
      reg.width = lineFields(9).toInt
      reg.seedCov = lineFields(10).toInt
      reg.secondary = lineFields(11).toInt
      reg.hash = lineFields(12).toLong

      line = reader.readLine
      var seqLen = line.toInt
      line = reader.readLine
      val seq = line.toCharArray.map(ele => ((ele.toByte) - 48).toByte)
      //seq.foreach(e => print(e))
      //println
      
      line = reader.readLine
      var vecSize = line.toInt
      var regs: Array[MemAlnRegType] = new Array[MemAlnRegType](vecSize)
      var j = 0
      while(j < vecSize) {
        line = reader.readLine
        val lineFields_ = line.split(" ")
        regs(j) = new MemAlnRegType
        regs(j).rBeg = lineFields_(0).toLong
        regs(j).rEnd = lineFields_(1).toLong
        regs(j).qBeg = lineFields_(2).toInt
        regs(j).qEnd = lineFields_(3).toInt
        regs(j).score = lineFields_(4).toInt
        regs(j).trueScore = lineFields_(5).toInt
        regs(j).sub = lineFields_(6).toInt
        regs(j).csub = lineFields_(7).toInt
        regs(j).subNum = lineFields_(8).toInt
        regs(j).width = lineFields_(9).toInt
        regs(j).seedCov = lineFields_(10).toInt
        regs(j).secondary = lineFields_(11).toInt
        regs(j).hash = lineFields_(12).toLong
        
        j += 1
      }

      line = reader.readLine
      // next input
      line = reader.readLine

      val ret = memMateSw(opt, pacLen, pac, pes, reg, seqLen, seq, regs)
      println("ret: " + ret._1 + "; size: " + ret._2.size)

      i += 1    
    }
  }

  def testMemPair(fileName: String, opt: MemOptType, pacLen: Long, pac: Array[Byte], pes: Array[MemPeStat]) {

    val reader = new BufferedReader(new FileReader(fileName))
    var line = reader.readLine
    var i = 0

    while(line != null) {
      var vecSize = line.toInt
      var regs0: Array[MemAlnRegType] = new Array[MemAlnRegType](vecSize)
      var j = 0
      while(j < vecSize) {
        line = reader.readLine
        val lineFields_ = line.split(" ")
        regs0(j) = new MemAlnRegType
        regs0(j).rBeg = lineFields_(0).toLong
        regs0(j).rEnd = lineFields_(1).toLong
        regs0(j).qBeg = lineFields_(2).toInt
        regs0(j).qEnd = lineFields_(3).toInt
        regs0(j).score = lineFields_(4).toInt
        regs0(j).trueScore = lineFields_(5).toInt
        regs0(j).sub = lineFields_(6).toInt
        regs0(j).csub = lineFields_(7).toInt
        regs0(j).subNum = lineFields_(8).toInt
        regs0(j).width = lineFields_(9).toInt
        regs0(j).seedCov = lineFields_(10).toInt
        regs0(j).secondary = lineFields_(11).toInt
        //regs0(j).hash = 0
        //regs0(j).hash = lineFields_(12).toLong
        val bigInt = new BigInteger(lineFields_(12))
        regs0(j).hash = bigInt.longValue
        //println("BigInt: " + bigInt + ", Long: " + regs0(j).hash.toBinaryString)

        j += 1
      }

      line = reader.readLine
      vecSize = line.toInt
      var regs1: Array[MemAlnRegType] = new Array[MemAlnRegType](vecSize)
      j = 0
      while(j < vecSize) {
        line = reader.readLine
        val lineFields_ = line.split(" ")
        regs1(j) = new MemAlnRegType
        regs1(j).rBeg = lineFields_(0).toLong
        regs1(j).rEnd = lineFields_(1).toLong
        regs1(j).qBeg = lineFields_(2).toInt
        regs1(j).qEnd = lineFields_(3).toInt
        regs1(j).score = lineFields_(4).toInt
        regs1(j).trueScore = lineFields_(5).toInt
        regs1(j).sub = lineFields_(6).toInt
        regs1(j).csub = lineFields_(7).toInt
        regs1(j).subNum = lineFields_(8).toInt
        regs1(j).width = lineFields_(9).toInt
        regs1(j).seedCov = lineFields_(10).toInt
        regs1(j).secondary = lineFields_(11).toInt
        //regs1(j).hash = lineFields_(12).toLong
        //regs1(j).hash = 0 
        val bigInt = new BigInteger(lineFields_(12))
        regs1(j).hash = bigInt.longValue
        //println("BigInt: " + bigInt + ", Long: " + regs1(j).hash.toBinaryString)

        j += 1
      }

      line = reader.readLine
      var id: Long = line.toLong     

      var alnRegVec = new Array[Array[MemAlnRegType]](2)
      alnRegVec(0) = regs0
      alnRegVec(1) = regs1
      if(regs0.size > 0 && regs1.size > 0) {
        var ret = memPair(opt, pacLen, pac, pes, alnRegVec, id)
        if(ret._1 > 0) {
          val z = ret._4     
          //println(ret._1 + " " + ret._2 + " " + ret._3 + " " + z(0) + " " + z(1))
        }
      }

      line = reader.readLine
      i += 1
    }
  }


  def testMemSamPe(fileName: String, opt: MemOptType, bns: BNTSeqType, pac: Array[Byte], pes: Array[MemPeStat]) {

    val reader = new BufferedReader(new FileReader(fileName))
    var line = reader.readLine
    var i = 0

    while(line != null) {
      var id = line.toInt
      
      line = reader.readLine
      val seqLen0 = line.toInt
      line = reader.readLine
      val seq0 = line

      line = reader.readLine
      val seqLen1 = line.toInt
      line = reader.readLine
      val seq1 = line

      line = reader.readLine
      var vecSize = line.toInt
      var regs0: Array[MemAlnRegType] = new Array[MemAlnRegType](vecSize)
      var j = 0
      while(j < vecSize) {
        line = reader.readLine
        val lineFields_ = line.split(" ")
        regs0(j) = new MemAlnRegType
        regs0(j).rBeg = lineFields_(0).toLong
        regs0(j).rEnd = lineFields_(1).toLong
        regs0(j).qBeg = lineFields_(2).toInt
        regs0(j).qEnd = lineFields_(3).toInt
        regs0(j).score = lineFields_(4).toInt
        regs0(j).trueScore = lineFields_(5).toInt
        regs0(j).sub = lineFields_(6).toInt
        regs0(j).csub = lineFields_(7).toInt
        regs0(j).subNum = lineFields_(8).toInt
        regs0(j).width = lineFields_(9).toInt
        regs0(j).seedCov = lineFields_(10).toInt
        regs0(j).secondary = lineFields_(11).toInt
        //regs0(j).hash = 0
        //regs0(j).hash = lineFields_(12).toLong
        val bigInt = new BigInteger(lineFields_(12))
        regs0(j).hash = bigInt.longValue
        //println("BigInt: " + bigInt + ", Long: " + regs0(j).hash.toBinaryString)

        j += 1
      }

      line = reader.readLine
      vecSize = line.toInt
      var regs1: Array[MemAlnRegType] = new Array[MemAlnRegType](vecSize)
      j = 0
      while(j < vecSize) {
        line = reader.readLine
        val lineFields_ = line.split(" ")
        regs1(j) = new MemAlnRegType
        regs1(j).rBeg = lineFields_(0).toLong
        regs1(j).rEnd = lineFields_(1).toLong
        regs1(j).qBeg = lineFields_(2).toInt
        regs1(j).qEnd = lineFields_(3).toInt
        regs1(j).score = lineFields_(4).toInt
        regs1(j).trueScore = lineFields_(5).toInt
        regs1(j).sub = lineFields_(6).toInt
        regs1(j).csub = lineFields_(7).toInt
        regs1(j).subNum = lineFields_(8).toInt
        regs1(j).width = lineFields_(9).toInt
        regs1(j).seedCov = lineFields_(10).toInt
        regs1(j).secondary = lineFields_(11).toInt
        //regs1(j).hash = lineFields_(12).toLong
        //regs1(j).hash = 0 
        val bigInt = new BigInteger(lineFields_(12))
        regs1(j).hash = bigInt.longValue
        //println("BigInt: " + bigInt + ", Long: " + regs1(j).hash.toBinaryString)

        j += 1
      }

      var seqs = new Array[FASTQSingleNode](2)
      seqs(0) = new FASTQSingleNode
      seqs(0).seq = seq0
      seqs(0).seqLen = seqLen0
      seqs(1) = new FASTQSingleNode
      seqs(1).seq = seq1
      seqs(1).seqLen = seqLen1
      var alnRegVec = new Array[Array[MemAlnRegType]](2)
      alnRegVec(0) = regs0
      alnRegVec(1) = regs1
      memSamPe(opt, bns, pac, pes, id, seqs, alnRegVec)

      line = reader.readLine
      i += 1
    }
  }
}

