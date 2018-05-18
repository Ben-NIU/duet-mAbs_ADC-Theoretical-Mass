## Return the chain formula depending on PTM input.(Not considering Glycan and disulfide bonds yet)

Calc.chain<-function(chain, PyroE, Lys){
 if(chain!=""){
  if(is.null(PyroE)){
    pyroe<-0
  } else {
    if(PyroE=="No"){
      pyroe<-0}
    if(PyroE=="Yes"){
     if(substr(chain, 1,1)=="E"){  
      pyroe<-c(C=0, H=-2, N=0, O=-1, S=0)} ## -18 Da
     else if(substr(chain, 1,1)=="Q") {
      pyroe<-c(C=0, H=-3, N=-1, O=0, S=0)} ## -17 Da
    }
  }
  
 if(is.null(Lys)){
   lys<-0
 } else {
   if(Lys=="No"){
     lys<-0
     } else if(Lys=="Yes") {
       lys<-c(C=-6, H=-12, N=-2, O=-1, S=0)}
 }
  
  total<-unlist(ConvertPeptide(chain, IAA=FALSE)) + pyroe + lys
 } else { total<-c("C"=0,"H"=0,"N"=0, "O"=0, "S"=0)}
  return(total)
}