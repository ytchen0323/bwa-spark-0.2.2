package cs.ucla.edu.bwaspark.datatype

import java.io.{FileInputStream, IOException}
import java.nio.channels.FileChannel
import cs.ucla.edu.bwaspark.datatype.BinaryFileReadUtil._

//BWAIdxType: maintaining all the information of BWA Index generated from FastA Reference
class BWAIdxType {  

  //1st: BWTType(".bwt", ".sa" files)
  var bwt: BWTType = _

  //2nd: BNTSeqType(".ann", ".amb" files)
  var bns: BNTSeqType = _

  //3rd: PACType(".pac" file)
  var pac: Array[Byte] = _  //uint8_t[]

  //loading files into fields
  //prefix: prefix of filenames
  //which: !!!to add!!!
  def load(prefix: String, which: Int) {
    //There is a function called "bwa_idx_infer_prefix" in original BWA,
    //but it seems to be useless

    //loading bwt
    //!!!In the beginning, set all as true
    //if (which & BWA_IDX_BWT) {
    if (true) {
      bwt = new BWTType
      bwt.load(prefix)
    }

    //loading bns
    //!!!In the beginning, set all as true
    //if (which & BWA_IDX_BNS) {
    if (true) {
      bns = new BNTSeqType
      bns.load(prefix)
      
      //loading pac
      //!!!In the beginning, set all as true
      //if (which & BWA_IDX_PAC) {
      if (true) {
        def pacLoader(filename: String, length: Long): Array[Byte] = {
          //to add: reading binary file
          val reader = new FileInputStream(filename).getChannel
          var pac = readByteArray(reader, (length/4+1).toInt, 0)          
          pac
        }
        pac = pacLoader(prefix+".pac", bns.l_pac)
      }
    }
  }
}
