package cs.ucla.edu.bwaspark.datatype

class MemAlnRegType {
  var rBeg: Long = 0       // [rBeg,rEnd): reference sequence in the alignment
  var rEnd: Long = 0       
  var qBeg: Int = 0        // [qBeg,qEnd): query sequence in the alignment
  var qEnd: Int = 0
  var score: Int = 0       // best local SW score
  var trueScore: Int = 0   // actual score corresponding to the aligned region; possibly smaller than $score
  var sub: Int = 0         // 2nd best SW score
  var csub: Int = 0        // SW score of a tandem hit
  var subNum: Int = 0      // approximate number of suboptimal hits
  var width: Int = 0       // actual band width used in extension
  var seedCov: Int = 0     // length of regions coverged by seeds
  var secondary: Int = 0   // index of the parent hit shadowing the current hit; <0 if primary
  var hash: Long = 0
}

