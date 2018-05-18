show.table.Red<-function(Seq.info, ch1.PyroE, ch1.Lys, ch2.PyroE, ch2.Lys, ch3.PyroE, ch3.Lys, ch4.PyroE, ch4.Lys, Glycan, un1=3,un2=3,un3=3,un4=3, which.chain, Carbon, Hydrogen, Nitrogen, Oxygen, Sulfer, how.many1, how.many2, how.many3, how.many4, w.warhead){
  source("Calc.chain.R")

  ## compute the theoretical formula depending on PTM input.
 GLYCAN<-list("none"=c(C=0,H=0,N=0,O=0,S=0),"Deglycosylated"=c(C=0, H=-1, N=-1, O=1, S=0),'G0'=c("C"=50, "H"=82, "N"=4, "O"=35, "S"=0),"G0F"=c(C=56, H=92, N=4, O=39, S=0),"G1"=c(C=56, H=92, N=4, O=40, S=0),"G1F"=c(C=62, H=102, N=4, O=44, S=0),"G2F"=c(C=68, H=112, N=4, O=49, S=0),"G0-GN"=c(C=42, H=69, N=3, O=30, S=0))

 test.glycan<-function(x){
   if(is.null(x)){
     "none"
   } else {
     x}
 }
 
 HOWMANY<-function(x){
   if(is.null(x)){
     0
   } else{
     x}
 }

 
if(w.warhead=="Yep"){
 
    chain1.form.thr<-Calc.chain(as.character(Seq.info$SQ[1]), ch1.PyroE, ch1.Lys) + GLYCAN[[test.glycan(Glycan)]] + c(C=0, H=un1*(-2), N=0, O=0, S=0) + (1 %in% which.chain)*HOWMANY(how.many1)*c(C=Carbon, H=Hydrogen, N=Nitrogen, O=Oxygen, S=Sulfer)

 
    chain2.form.thr<-Calc.chain(as.character(Seq.info$SQ[2]), ch2.PyroE, ch2.Lys) + GLYCAN[[test.glycan(Glycan)]] + c(C=0, H=un2*(-2), N=0, O=0, S=0) + (2 %in% which.chain)*HOWMANY(how.many2)*c(C=Carbon, H=Hydrogen, N=Nitrogen, O=Oxygen, S=Sulfer)

  
    chain3.form.thr<-Calc.chain(as.character(Seq.info$SQ[3]), ch3.PyroE, ch3.Lys) + GLYCAN[[test.glycan(Glycan)]] + c(C=0, H=un3*(-2), N=0, O=0, S=0) + (3 %in% which.chain)*HOWMANY(how.many3)*c(C=Carbon, H=Hydrogen, N=Nitrogen, O=Oxygen, S=Sulfer)

  
    chain4.form.thr<-Calc.chain(as.character(Seq.info$SQ[4]), ch4.PyroE, ch4.Lys) + GLYCAN[[test.glycan(Glycan)]] + c(C=0, H=un4*(-2), N=0, O=0, S=0) + (4 %in% which.chain)*HOWMANY(how.many4)*c(C=Carbon, H=Hydrogen, N=Nitrogen, O=Oxygen, S=Sulfer)
} else {
  
  chain1.form.thr<-Calc.chain(as.character(Seq.info$SQ[1]), ch1.PyroE, ch1.Lys) + GLYCAN[[test.glycan(Glycan)]] + c(C=0, H=un1*(-2), N=0, O=0, S=0) 
  
  
  chain2.form.thr<-Calc.chain(as.character(Seq.info$SQ[2]), ch2.PyroE, ch2.Lys) + GLYCAN[[test.glycan(Glycan)]] + c(C=0, H=un2*(-2), N=0, O=0, S=0) 
  
  
  chain3.form.thr<-Calc.chain(as.character(Seq.info$SQ[3]), ch3.PyroE, ch3.Lys) + GLYCAN[[test.glycan(Glycan)]] + c(C=0, H=un3*(-2), N=0, O=0, S=0)
  
  
  chain4.form.thr<-Calc.chain(as.character(Seq.info$SQ[4]), ch4.PyroE, ch4.Lys) + GLYCAN[[test.glycan(Glycan)]] + c(C=0, H=un4*(-2), N=0, O=0, S=0)
}

  
 #===========================================================#
  Standard<-c("C"=12.01078, "H"=1.007947, "N"=14.00672, "O"=15.99943, "S"=32.0655)
  ## compute the AVERAGE mass of theoretical formula depending on PTM input.
  
  chain1.mass<-format(round(sum(chain1.form.thr*Standard),1), nsmall = 1)  ## REPORT
  chain2.mass<-format(round(sum(chain2.form.thr*Standard),1), nsmall = 1)  ## REPORT
  chain3.mass<-format(round(sum(chain3.form.thr*Standard),1), nsmall = 1)  ## REPORT
  chain4.mass<-format(round(sum(chain4.form.thr*Standard),1), nsmall = 1)  ## REPORT
  
  #===============================================================================#
  chain1.info<-data.frame( "C"=as.integer(chain1.form.thr[["C"]]),"H"=as.integer(chain1.form.thr[["H"]]),"N"=as.integer(chain1.form.thr[["N"]]),"O"=as.integer(chain1.form.thr[["O"]]),"S"=as.integer(chain1.form.thr[["S"]]), "Mass (Da)"=chain1.mass, "Description"="Expected", check.names = FALSE)
                           
  chain2.info<-data.frame( "C"=as.integer(chain2.form.thr[["C"]]),"H"=as.integer(chain2.form.thr[["H"]]),"N"=as.integer(chain2.form.thr[["N"]]),"O"=as.integer(chain2.form.thr[["O"]]),"S"=as.integer(chain2.form.thr[["S"]]), "Mass (Da)"=chain2.mass, "Description"="Expected", check.names = FALSE)
  chain3.info<-data.frame( "C"=as.integer(chain3.form.thr[["C"]]),"H"=as.integer(chain3.form.thr[["H"]]),"N"=as.integer(chain3.form.thr[["N"]]),"O"=as.integer(chain3.form.thr[["O"]]),"S"=as.integer(chain3.form.thr[["S"]]), "Mass (Da)"=chain3.mass, "Description"="Expected", check.names = FALSE)
  chain4.info<-data.frame( "C"=as.integer(chain4.form.thr[["C"]]),"H"=as.integer(chain4.form.thr[["H"]]),"N"=as.integer(chain4.form.thr[["N"]]),"O"=as.integer(chain4.form.thr[["O"]]),"S"=as.integer(chain4.form.thr[["S"]]), "Mass (Da)"=chain4.mass, "Description"="Expected", check.names = FALSE)
  
  
  Overall<-list(chain1=chain1.info, chain2=chain2.info, chain3=chain3.info, chain4=chain4.info)
  return(Overall)
}
