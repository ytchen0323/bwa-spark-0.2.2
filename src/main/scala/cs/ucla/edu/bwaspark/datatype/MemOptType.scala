package cs.ucla.edu.bwaspark.datatype 

import scala.math.log
//import scalaz._

class MemOptType
{
	var a : Int = 1;
	var b : Int = 4;
	var oDel : Int = 6;
	var eDel : Int = 1;
	var oIns : Int = 6;
	var eIns : Int = 1;
	var penUnpaired : Int  = 17;
	var penClip5 : Int  = 5;
	var penClip3 : Int  = 5;
	var w : Int = 100;
	var zdrop : Int = 100;

	var T : Int = 30;
	var flag: Int  = 0;
	var minSeedLen: Int  = 19;
	var splitFactor: Float  = 1.5f;
	var splitWidth: Int  = 10;
	var maxOcc: Int  = 10000;
	var maxChainGap: Int  = 10000;

	var chunkSize: Int  = 10000000;
	var maskLevel: Float  = 0.50f;
	var chainDropRatio: Float  = 0.50f;
	var maskLevelRedun: Float  = 0.95f;
	var mapQCoefLen: Float  = 50f;
	var mapQCoefFac: Int  = log(mapQCoefLen).asInstanceOf[Int]
	var maxIns: Int  = 10000;
	var maxMatesw: Int  = 100;
	var mat: Array[Byte]  = new Array[Byte](25);//all initialized to 0

	private def bwaFillScmat(){
	  var k = 0
	  for(i <- 0 to 3){
	  	for(j <- 0 to 3){
		var temp = ( if (i == j)  a else -b)
		
	  	  mat(k) = temp.asInstanceOf[Byte]
	  	  k = k + 1
	  	}
	  	mat(k) = -1
	  	k = k + 1
	  }
	  for(j <- 0 to 4){
	  	mat(k) = -1
	  	k = k + 1
	  }

	}
	
	def load() {
		bwaFillScmat()
	}

}
