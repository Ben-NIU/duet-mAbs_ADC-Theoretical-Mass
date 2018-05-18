read.duet.fasta<-function(x){
  r<-suppressWarnings(readLines(x))
  LH<-which(substr(r,1,1)==">")
  NM<-r[LH]
  name<-substr(NM, 2, 100)
  SQ<-r[LH+1]
  
  ct.lth<-function(x){
    if(x!=""){
    length(strsplit(as.character(x), split="")[[1]])
    } else {0}
  }
  
  L<-do.call("c",lapply(SQ, ct.lth))
  
  Seq<-data.frame("Name"=name, "SQ"=SQ, "LTH"=L)
  return(Seq)
}