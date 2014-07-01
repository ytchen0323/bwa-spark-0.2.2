package cs.ucla.edu.bwaspark

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
//import org.apache.spark.storage.StorageLevel

import scala.collection.mutable.MutableList

import cs.ucla.edu.bwaspark.datatype._
import cs.ucla.edu.bwaspark.worker1.BWAMemWorker1._
import cs.ucla.edu.bwaspark.worker2.BWAMemWorker2._
import cs.ucla.edu.bwaspark.sam.SAMHeader._
import cs.ucla.edu.bwaspark.sam.SAMWriter
import cs.ucla.edu.bwaspark.debug.DebugFlag._
import cs.ucla.edu.bwaspark.fastq._

import java.io.FileReader
import java.io.BufferedReader

object BWAMEMSpark {
  // load reads from the FASTQ file (for testing use)
  private def loadFASTQSeqs(fileName: String, readNum: Int): Array[FASTQSingleNode] = {
    
    val reader = new BufferedReader(new FileReader(fileName))
    var line = reader.readLine
    var i = 0
    var readIdx = 0    
    var seqs = new Array[FASTQSingleNode](readNum / 4)

    while(line != null) {
      //if(i % 4 == 1) {
      //  seqs(readIdx) = line
      //  readIdx += 1
      //}
      //i += 1
      //line = reader.readLine

      val lineFields = line.split(" ")
      seqs(i) = new FASTQSingleNode
  
      if(lineFields.length == 1) {
        if(lineFields(0).charAt(0) == '@') seqs(i).name = lineFields(0).substring(1).dropRight(2)
        else seqs(i).name = lineFields(0).dropRight(2)
        seqs(i).seq = reader.readLine
        seqs(i).seqLen = seqs(i).seq.size
        reader.readLine
        seqs(i).qual = reader.readLine
        seqs(i).comment = ""
        i += 1
      }
      else if(lineFields.length == 2) {
        if(lineFields(0).charAt(0) == '@') seqs(i).name = lineFields(0).substring(1).dropRight(2)
        else seqs(i).name = lineFields(0).dropRight(2)
        seqs(i).comment = lineFields(1)
        seqs(i).seq = reader.readLine
        seqs(i).seqLen = seqs(i).seq.size
        reader.readLine
        seqs(i).qual = reader.readLine
        i += 1
      }
      else 
        println("Error: Input format not handled")

      line = reader.readLine
    } 
    
    //seqs.foreach(s => println(s.seq))
    seqs
  }

  class testRead {
    var seq: FASTQSingleNode = _
    var regs: Array[MemAlnRegType] = _
  }

  class testSWAlnInput {
    var qLen: Int = _
    var tLen: Int = _
    var oDel: Int = _
    var eDel: Int = _
    var oIns: Int = _
    var eIns: Int = _
    var xtra: Int = _
    var query: Array[Byte] = _
    var target: Array[Byte] = _
  }

  private def loadSWAlnInput(reader: BufferedReader, inputNum: Int): Array[testSWAlnInput] = {
    var line = reader.readLine    
    var i = 0
    var swInputArray = new Array[testSWAlnInput](inputNum)
 
    while(line != null) {
      val lineFields = line.split(" ")
      swInputArray(i) = new testSWAlnInput
      swInputArray(i).qLen = lineFields(0).toInt
      swInputArray(i).tLen = lineFields(1).toInt
      swInputArray(i).oDel = lineFields(2).toInt
      swInputArray(i).eDel = lineFields(3).toInt
      swInputArray(i).oIns = lineFields(4).toInt
      swInputArray(i).eIns = lineFields(5).toInt
      swInputArray(i).xtra = lineFields(6).toInt
      swInputArray(i).query = new Array[Byte](swInputArray(i).qLen)
      swInputArray(i).target = new Array[Byte](swInputArray(i).tLen)

      var j = 0
      while(j < swInputArray(i).qLen) {
        swInputArray(i).query(j) = (lineFields(7).charAt(j).toInt - 48).toByte
        j += 1
      }

      j = 0
      while(j < swInputArray(i).tLen) {
        swInputArray(i).target(j) = (lineFields(8).charAt(j).toInt - 48).toByte
        j += 1
      }

      line = reader.readLine
      i += 1
    }

    swInputArray
  }

  def main(args: Array[String]) {
/*
    //val sc = new SparkContext("local[12]", "BWA-mem Spark",
       //"/home/hadoopmaster/spark/spark-0.9.0-incubating-bin-hadoop2-prebuilt/", List("/home/ytchen/incubator/bwa-spark-0.2.0/target/bwa-spark-0.2.0.jar"))
    //val sc = new SparkContext("spark://Jc11:7077", "BWA-mem Spark",
       //"/home/hadoopmaster/spark/spark-0.9.0-incubating-bin-hadoop2-prebuilt/", List("/home/ytchen/incubator/bwa-spark-0.2.0/target/bwa-spark-0.2.0.jar"))

    //val fastqLoader = new FASTQLocalFileLoader()
    //fastqLoader.storeFASTQInHDFS(sc, args(0), args(1))
    //val fastqRDDLoader = new FASTQRDDLoader(sc, "hdfs://Jc11:9000/user/ytchen/ERR013140_2.filt.fastq.test4/", 13)
    //val fastqRDD = fastqRDDLoader.RDDLoadAll()
    //val fastqRDD = fastqRDDLoader.RDDLoad("hdfs://Jc11:9000/user/ytchen/ERR013140_2.filt.fastq.test4/2/")

    if(bwaSetReadGroup("@RG\tID:HCC1954\tLB:HCC1954\tSM:HCC1954")) {
      println("Read line: " + readGroupLine)
      println("Read Group ID: " + bwaReadGroupID)
    }
    else println("Error on reading header")

    //loading index files
    println("Load Index Files")
    val bwaIdx = new BWAIdxType
    val prefix = "/home/pengwei/genomics/ReferenceMetadata/human_g1k_v37.fasta"
    bwaIdx.load(prefix, 0)

    //loading BWA MEM options
    println("Load BWA-MEM options")
    val bwaMemOpt = new MemOptType
    bwaMemOpt.load

    //debugLevel = 1

    //loading reads
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_1read.fq")
    println("Load FASTQ files")

    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_20reads.fq", 80)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_100reads.fq", 400)
    var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_10Mreads.fq", 10000000)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_5reads_err.fq")
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_1read_err.fq", 4)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_15000reads.fq", 60000)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/ERR013140_1.filt.fastq")

    //val regsAllReads = seqs.map( seq => bwaMemWorker1(bwaMemOpt, bwaIdx.bwt, bwaIdx.bns, bwaIdx.pac, null, seq.length, seq) )
    //println("Processing")
    
    println("@Worker1")
    var i = 0
    var regsAllReads: Array[Array[MemAlnRegType]] = new Array[Array[MemAlnRegType]](seqs.length)
    //val regsAllReads = seqs.map( {
    //seqs.foreach( {
      //seq => regsAllReads(i) = bwaMemWorker1(bwaMemOpt, bwaIdx.bwt, bwaIdx.bns, bwaIdx.pac, null, seq.seqLen, seq.seq) 
      //if(i >= 14700) debugLevel = 1
      //if(i >= 14700) println(i)
      //println("Read: " + i)
      
      //i += 1
      //if((i % 10000) == 0) println(i)
      //} )

    while(i < seqs.length) {
      regsAllReads(i) = bwaMemWorker1(bwaMemOpt, bwaIdx.bwt, bwaIdx.bns, bwaIdx.pac, null, seqs(i).seqLen, seqs(i).seq)
      i += 1
      if((i % 10000) == 0) println(i)
    }

    var testReads = new Array[testRead](seqs.length)
    for(i <- 0 to (seqs.length - 1)) {
      var read = new testRead
      read.seq = seqs(i)
      read.regs = regsAllReads(i)
      testReads(i) = read
    }

    println("@Worker2")
    i = 0
    testReads.foreach(read => {
      //println
      //println("Read " + i)
      bwaMemWorker2(bwaMemOpt, read.regs, bwaIdx.bns, bwaIdx.pac, read.seq, 0) 
      i += 1
      if((i % 10000) == 0) println(i)
      } )

    val samWriter = new SAMWriter("test.sam")
    samWriter.init
    samWriter.writeString(bwaGenSAMHeader(bwaIdx.bns))
    testReads.foreach(r => samWriter.writeString((r.seq.sam)))
    samWriter.close
*/

// Testing: Pair-End SW
    val reader = new BufferedReader(new FileReader("/home/ytchen/bwa/bwa-0.7.8/sw.input"))
    val inputSW = loadSWAlnInput(reader,100)
    inputSW.foreach(sw => {
      print(sw.qLen + " " + sw.tLen + " " + sw.oDel + " " + sw.eDel + " " + sw.oIns + " " + sw.eIns + " " + sw.xtra + " ")
      sw.query.foreach(print(_))
      print(" ")
      sw.target.foreach(print(_))
      println } )

// Testing
/*
    var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_20reads.fq", 80)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_100reads.fq", 400)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_1read_No3.fq", 4)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_1read_No14.fq", 4)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_1read_No33.fq", 4)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_1read_No32.fq", 4)
    //var seqs = loadFASTQSeqs("/home/ytchen/genomics/data/HCC1954_1_1read_No12.fq", 4)
    val regsAllReads = seqs.map(seq => bwaMemWorker1(bwaMemOpt, bwaIdx.bwt, bwaIdx.bns, bwaIdx.pac, null, seq.seqLen, seq.seq))


    // print regs for all reads
//    var readNum = 0
//    regsAllReads.foreach(read => {
//      var i = 0
//     println("#####")
//      println("Read " + readNum)
//      read.foreach(r => {
//        print("Reg " + i + "(")
//        print(r.rBeg + ", " + r.rEnd + ", " + r.qBeg + ", " + r.qEnd + ", " + r.score + ", " + r.trueScore + ", ")
//        println(r.sub + ", "  + r.csub + ", " + r.subNum + ", " + r.width + ", " + r.seedCov + ", " + r.secondary + ")")
//        i += 1
//      } )
//      readNum += 1
//    } )



    var testReads = new MutableList[testRead]
    for(i <- 0 to (seqs.length - 1)) {
      var read = new testRead
      read.seq = seqs(i)
      read.regs = regsAllReads(i)
      testReads += read
    }

    var i = 0
    testReads.foreach(read => {
      //println
      //println("Read " + i)
      bwaMemWorker2(bwaMemOpt, read.regs, bwaIdx.bns, bwaIdx.pac, read.seq, 0) 
      i += 1
      } )

    val samWriter = new SAMWriter("test.sam")
    samWriter.init
    samWriter.writeString(bwaGenSAMHeader(bwaIdx.bns))
    testReads.foreach(r => samWriter.writeString((r.seq.sam)))
    samWriter.close
*/
  } 
}
