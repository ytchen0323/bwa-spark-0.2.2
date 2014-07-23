package cs.ucla.edu.bwaspark

import cs.ucla.edu.bwaspark.datatype._
import cs.ucla.edu.bwaspark.worker1.BWAMemWorker1._
import cs.ucla.edu.bwaspark.worker2.BWAMemWorker2._
import cs.ucla.edu.bwaspark.worker2.MemSamPe._
import cs.ucla.edu.bwaspark.sam.SAMHeader._
import cs.ucla.edu.bwaspark.sam.SAMWriter
import cs.ucla.edu.bwaspark.debug.DebugFlag._
import cs.ucla.edu.bwaspark.fastq._
import cs.ucla.edu.bwaspark.util.SWUtil._

// for test use
import cs.ucla.edu.bwaspark.worker2.MemRegToADAMSAM._

import java.io.FileReader
import java.io.BufferedReader

object FastMap {
  private val MEM_F_PE: Int = 0x2
  private val MEM_F_ALL = 0x8
  private val MEM_F_NO_MULTI = 0x10

  // load reads from the FASTQ file (for testing use)
  def loadBatchFASTQSeqs(reader: BufferedReader, batchNum: Int): Array[FASTQSingleNode] = {
    
    var line = reader.readLine
    var i = 0
    var readIdx = 0    
    var seqs: Vector[FASTQSingleNode] = scala.collection.immutable.Vector.empty

    while(line != null && i < batchNum) {

      val lineFields = line.split(" ")
      var seq = new FASTQSingleNode
  
      if(lineFields.length == 1) {
        if(lineFields(0).charAt(0) == '@') seq.name = lineFields(0).substring(1).dropRight(2)
        else seq.name = lineFields(0).dropRight(2)
        seq.seq = reader.readLine
        seq.seqLen = seq.seq.size
        reader.readLine
        seq.qual = reader.readLine
        seq.comment = ""
        seqs = seqs :+ seq
        i += 1
      }
      else if(lineFields.length == 2) {
        if(lineFields(0).charAt(0) == '@') seq.name = lineFields(0).substring(1).dropRight(2)
        else seq.name = lineFields(0).dropRight(2)
        seq.comment = lineFields(1)
        seq.seq = reader.readLine
        seq.seqLen = seq.seq.size
        reader.readLine
        seq.qual = reader.readLine
        seqs = seqs :+ seq
        i += 1
      }
      else 
        println("Error: Input format not handled")

      if(i < batchNum)
        line = reader.readLine
    } 
    
    seqs.toArray
  }

  def loadBatchPairFASTQSeqs(reader1: BufferedReader, reader2: BufferedReader, batchNum: Int): Array[FASTQSingleNode] = {

    var line = reader1.readLine
    var i = 0
    var readIdx = 0    
    var seqs: Vector[FASTQSingleNode] = scala.collection.immutable.Vector.empty

    while(line != null && i < batchNum) {

      val lineFields = line.split(" ")
      var seq = new FASTQSingleNode
  
      if(lineFields.length == 1) {
        if(lineFields(0).charAt(0) == '@') seq.name = lineFields(0).substring(1).dropRight(2)
        else seq.name = lineFields(0).dropRight(2)
        seq.seq = reader1.readLine
        seq.seqLen = seq.seq.size
        reader1.readLine
        seq.qual = reader1.readLine
        seq.comment = ""
        seqs = seqs :+ seq
        i += 1
      }
      else if(lineFields.length == 2) {
        if(lineFields(0).charAt(0) == '@') seq.name = lineFields(0).substring(1).dropRight(2)
        else seq.name = lineFields(0).dropRight(2)
        seq.comment = lineFields(1)
        seq.seq = reader1.readLine
        seq.seqLen = seq.seq.size
        reader1.readLine
        seq.qual = reader1.readLine
        seqs = seqs :+ seq
        i += 1
      }
      else 
        println("Error: Input format not handled")

      line = reader2.readLine

      if(line == null) println("Error: the number of two FASTQ files are different")
      else {
        val lineFields = line.split(" ")
        var seq = new FASTQSingleNode

        if(lineFields.length == 1) {
          if(lineFields(0).charAt(0) == '@') seq.name = lineFields(0).substring(1).dropRight(2)
          else seq.name = lineFields(0).dropRight(2)
          seq.seq = reader2.readLine
          seq.seqLen = seq.seq.size
          reader2.readLine
          seq.qual = reader2.readLine
          seq.comment = ""
          seqs = seqs :+ seq
          i += 1
        }
        else if(lineFields.length == 2) {
          if(lineFields(0).charAt(0) == '@') seq.name = lineFields(0).substring(1).dropRight(2)
          else seq.name = lineFields(0).dropRight(2)
          seq.comment = lineFields(1)
          seq.seq = reader2.readLine
          seq.seqLen = seq.seq.size
          reader2.readLine
          seq.qual = reader2.readLine
          seqs = seqs :+ seq
          i += 1
        }
        else
          println("Error: Input format not handled")
      }

      if(i < batchNum)
        line = reader1.readLine
    } 
    
    seqs.toArray

  }

  class testRead {
    var seq: FASTQSingleNode = _
    var regs: MemAlnRegArrayType = _
  }

  def memMain() {

    if(bwaSetReadGroup("@RG\tID:HCC1954\tLB:HCC1954\tSM:HCC1954")) {
      println("Read line: " + readGroupLine)
      println("Read Group ID: " + bwaReadGroupID)
    }
    else println("Error on reading header")

    //loading index files
    println("Load Index Files")
    val bwaIdx = new BWAIdxType
    val prefix = "/home/hadoopmaster/genomics/ReferenceMetadata/human_g1k_v37.fasta"
    bwaIdx.load(prefix, 0)

    //loading BWA MEM options
    println("Load BWA-MEM options")
    val bwaMemOpt = new MemOptType
    bwaMemOpt.load

    //loading reads
    println("Load FASTQ files")

    // set the single/pair end mode
    bwaMemOpt.flag |= MEM_F_PE
    bwaMemOpt.flag |= MEM_F_ALL
    bwaMemOpt.flag |= MEM_F_NO_MULTI

    // different input sizes
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_20reads.fq", 80)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_100reads.fq", 400)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_10Mreads.fq", 10000000)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_5reads_err.fq")
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_1read_err.fq", 4)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_15000reads.fq", 60000)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/ERR013140_1.filt.fastq")

    var batchNum = 99010
    var n = 0
    var numProcessed = 0
    var readNum = 0
    var seqs: Array[FASTQSingleNode] = new Array[FASTQSingleNode](0)
    var reader1: BufferedReader = null
    var reader2: BufferedReader = null

    if((bwaMemOpt.flag & MEM_F_PE) > 0) {
      reader1 = new BufferedReader(new FileReader("/home/ytchen/genomics/data/HCC1954_1_10Mreads.fq"))
      reader2 = new BufferedReader(new FileReader("/home/ytchen/genomics/data/HCC1954_2_10Mreads.fq"))
    }
    else
      reader1 = new BufferedReader(new FileReader("/home/ytchen/genomics/data/HCC1954_1_20reads.fq"))      

    val samWriter = new SAMWriter("test.sam")
    samWriter.init
    samWriter.writeString(bwaGenSAMHeader(bwaIdx.bns))

    do {
      if((bwaMemOpt.flag & MEM_F_PE) > 0)
        seqs = loadBatchPairFASTQSeqs(reader1, reader2, batchNum)
      else      
        seqs = loadBatchFASTQSeqs(reader1, batchNum)

      // Debugging
      n = seqs.size
      
      var bpNum: Long = 0
      var i = 0   
      while(i < n) {
        bpNum += seqs(i).seqLen
        i += 1
      }
 
      println("read " + n + " sequences (" + bpNum + " bp)")
      //seqs.foreach(s => println(s.seq))

      memProcessSeqs(bwaMemOpt, bwaIdx.bwt, bwaIdx.bns, bwaIdx.pac, numProcessed, n, seqs, null, samWriter)
      numProcessed += n
      
      println("Num processed: " + numProcessed)

    } while(seqs.size == batchNum)
    
    if((bwaMemOpt.flag & MEM_F_PE) > 0) {
      reader1.close
      reader2.close
    }
    else
      reader1.close
    samWriter.close
  } 


  def memProcessSeqs(opt: MemOptType, bwt: BWTType, bns: BNTSeqType, pac: Array[Byte], numProcessed: Long, n: Int, seqs: Array[FASTQSingleNode], pes0: Array[MemPeStat], samWriter: SAMWriter) {

    var pes: Array[MemPeStat] = new Array[MemPeStat](4)   
    var j = 0
    while(j < 4) {
      pes(j) = new MemPeStat
      j += 1
    }

    //println("PE: " + (opt.flag & MEM_F_PE))

    // worker1
    // find mapping positions
    println("@Worker1")
    var i = 0
    var regsAllReads: Array[MemAlnRegArrayType] = new Array[MemAlnRegArrayType](seqs.length)

    if((opt.flag & MEM_F_PE) == 0) {
      while(i < seqs.length) {
        regsAllReads(i) = bwaMemWorker1(opt, bwt, bns, pac, seqs(i).seqLen, seqs(i).seq)
        i += 1
        if((i % 10000) == 0) println(i)
      }
    }
    else {
      while(i < seqs.length) {
        regsAllReads(i) = bwaMemWorker1(opt, bwt, bns, pac, seqs(i).seqLen, seqs(i).seq)
        i += 1
        regsAllReads(i) = bwaMemWorker1(opt, bwt, bns, pac, seqs(i).seqLen, seqs(i).seq)
        i += 1

        if((i % 10000) == 0) println(i)
      }
    }

    var testReads = new Array[testRead](seqs.length)
    i = 0
    while(i < seqs.length) {
      var read = new testRead
      read.seq = seqs(i)
      read.regs = regsAllReads(i)
      testReads(i) = read
      i += 1
    }

   
    println("@MemPeStat")
    if((opt.flag & MEM_F_PE) > 0) { // infer insert sizes if not provided
      if(pes0 != null) // if pes0 != NULL, set the insert-size distribution as pes0
        pes = pes0
      else // otherwise, infer the insert size distribution from data
        memPeStat(opt, bns.l_pac, n, regsAllReads, pes)
    }    

/*
    println("pes(0): " + pes(0).low + " " + pes(0).high + " " + pes(0).failed + " " + pes(0).avg + " " + pes(0).std)
    println("pes(1): " + pes(1).low + " " + pes(1).high + " " + pes(1).failed + " " + pes(1).avg + " " + pes(1).std)
    println("pes(2): " + pes(2).low + " " + pes(2).high + " " + pes(2).failed + " " + pes(2).avg + " " + pes(2).std)
    println("pes(3): " + pes(3).low + " " + pes(3).high + " " + pes(3).failed + " " + pes(3).avg + " " + pes(3).std)
*/
/*
    pes(0).low = 0 
    pes(0).high = 0 
    pes(0).failed = 1 
    pes(0).avg = 0.0 
    pes(0).std = 0.0
    pes(1).low = 1
    pes(1).high = 1066 
    pes(1).failed = 0
    pes(1).avg = 253.625 
    pes(1).std = 113.5113109268558
    pes(2).low = 0 
    pes(2).high = 0 
    pes(2).failed = 1 
    pes(2).avg = 0.0 
    pes(2).std = 0.0
    pes(3).low = 0 
    pes(3).high = 0 
    pes(3).failed = 1 
    pes(3).avg = 0.0 
    pes(3).std = 0.0
*/

    println("@Worker2")
    i = 0
    if((opt.flag & MEM_F_PE) == 0) {
      testReads.foreach(read => {
        bwaMemWorker2(opt, read.regs.regs, bns, pac, read.seq, 0) 
        i += 1
        if((i % 10000) == 0) println(i) } )
    }
    else {
      while(i < testReads.size) {
        var alnRegVec: Array[Array[MemAlnRegType]] = new Array[Array[MemAlnRegType]](2)
        var seqs: Array[FASTQSingleNode] = new Array[FASTQSingleNode](2)
      
        if(testReads(i).regs != null) 
          alnRegVec(0) = testReads(i).regs.regs
        seqs(0) = testReads(i).seq
        if(testReads(i + 1).regs != null) 
          alnRegVec(1) = testReads(i + 1).regs.regs
        seqs(1) = testReads(i + 1).seq
        bwaMemWorker2Pair(opt, alnRegVec, bns, pac, seqs, numProcessed + i, pes)      

        i += 2
        if((i % 1000) == 0) println(i)
      }
    }

    testReads.foreach(r => samWriter.writeString((r.seq.sam)))

/*    
    println("@testMemMateSw")
    testMemMateSw("/home/ytchen/bwa/bwa-0.7.8/matesw_many.input", opt, bns.l_pac, pac, pes)
*/
/*
    println("@testMemPair")
    testMemPair("/home/ytchen/bwa/bwa-0.7.8/mempair.input", opt, bns.l_pac, pac, pes)
*/
/*   
    println("@testBwaFixXref2")
    testBwaFixXref2("/home/ytchen/bwa/bwa-0.7.8/bwafixxref2.input", opt, bns, pac)
*/
/*
    println("@testMemSamPe")
    testMemSamPe("/home/ytchen/bwa/bwa-0.7.8/memsampe.input", opt, bns, pac, pes)
*/
  }
} 
