package cs.ucla.edu.bwaspark.datatype

class BNTAnnType(offset_l: Long,
                 len_i: Int,
                 n_ambs_i: Int,
                 gi_i: Int,  //uint32_t
                 name_s: String,
                 anno_s: String) {
  var offset = offset_l
  var len = len_i
  var n_ambs = n_ambs_i
  var gi = gi_i
  var name = name_s
  var anno = anno_s
}
