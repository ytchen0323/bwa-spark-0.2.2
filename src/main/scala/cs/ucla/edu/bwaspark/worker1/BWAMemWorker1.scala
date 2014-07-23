package cs.ucla.edu.bwaspark.worker1

import cs.ucla.edu.bwaspark.datatype._
import scala.collection.mutable.MutableList
import java.util.TreeSet
import java.util.Comparator
import cs.ucla.edu.bwaspark.worker1.MemChain._
import cs.ucla.edu.bwaspark.worker1.MemChainFilter._
import cs.ucla.edu.bwaspark.worker1.MemChainToAlign._
import cs.ucla.edu.bwaspark.worker1.MemSortAndDedup._

//this standalone object defines the main job of BWA MEM:
//1)for each read, generate all the possible seed chains
//2)using SW algorithm to extend each chain to all possible aligns
object BWAMemWorker1 {
  
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

  //the function which do the main task
  def bwaMemWorker1(opt: MemOptType, //BWA MEM options
                    bwt: BWTType, //BWT and Suffix Array
                    bns: BNTSeqType, //.ann, .amb files
                    pac: Array[Byte], //.pac file uint8_t
                    len: Int, //the length of the read
                    seq: String //a read
                    ): MemAlnRegArrayType = { //all possible alignment  

    //println(seq)
    val read: Array[Byte] = seq.toCharArray.map(ele => locusEncode(ele))

    //first step: generate all possible MEM chains for this read
    val chains = generateChains(opt, bwt, bns.l_pac, len, read) 

    //second step: filter chains
    val chainsFiltered = memChainFilter(opt, chains)

    if (chainsFiltered == null) 
      null
    else {
      // build the references of the seeds in each chain
      var totalSeedNum = 0
      chainsFiltered.foreach(chain => {
        totalSeedNum += chain.seeds.length
        } )

      //third step: for each chain, from chain to aligns
      var regArray = new MemAlnRegArrayType
      regArray.maxLength = totalSeedNum
      regArray.regs = new Array[MemAlnRegType](totalSeedNum)

        // test
/*
        var tmp = chainsFiltered(2)
        chainsFiltered(2) = chainsFiltered(1)
        chainsFiltered(1) = tmp
*/

      for (i <- 0 until chainsFiltered.length) {
        memChainToAln(opt, bns.l_pac, pac, len, read, chainsFiltered(i), regArray)
      }

      regArray.regs = regArray.regs.filter(r => (r != null))
      regArray.maxLength = regArray.regs.length
      assert(regArray.curLength == regArray.maxLength, "[Error] After filtering array elements")

/*
    println("#####")
    var i = 0
    regArray.regs.foreach(r => {
        print("Reg " + i + "(")
        print(r.rBeg + ", " + r.rEnd + ", " + r.qBeg + ", " + r.qEnd + ", " + r.score + ", " + r.trueScore + ", ")
        println(r.sub + ", "  + r.csub + ", " + r.subNum + ", " + r.width + ", " + r.seedCov + ", " + r.secondary + ")")
        i += 1
    } )
*/
      //last step: sorting and deduplication
      val pureRegArray = memSortAndDedup(regArray, opt.maskLevelRedun)

      pureRegArray
    }
  }

}
