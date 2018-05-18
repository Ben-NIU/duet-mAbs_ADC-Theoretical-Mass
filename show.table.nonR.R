show.table.nonR<-function(Seq.info, ch1.PyroE, ch1.Lys, ch2.PyroE, ch2.Lys, ch3.PyroE, ch3.Lys, ch4.PyroE, ch4.Lys, Glycan1,Glycan2, Glycan3, Glycan4, total.DS, cys, warhead,Carb=0, Hydro=0, Nitro=0, Oxy=0, Sulf=0, w.warhead){
  source("Calc.chain.R")

  ## compute the theoretical formula depending on PTM input.
  chain1.form.thr<-Calc.chain(as.character(Seq.info$SQ[1]), ch1.PyroE, ch1.Lys)
  chain2.form.thr<-Calc.chain(as.character(Seq.info$SQ[2]), ch2.PyroE, ch2.Lys)
  chain3.form.thr<-Calc.chain(as.character(Seq.info$SQ[3]), ch3.PyroE, ch3.Lys)
  chain4.form.thr<-Calc.chain(as.character(Seq.info$SQ[4]), ch4.PyroE, ch4.Lys) 
  
  GLYCAN<-list("none"=c(C=0,H=0,N=0,O=0,S=0),"Deglycosylated"=c(C=0, H=-1, N=-1, O=1, S=0),'G0'=c("C"=50, "H"=82, "N"=4, "O"=35, "S"=0),"G0F"=c(C=56, H=92, N=4, O=39, S=0),"G1"=c(C=56, H=92, N=4, O=40, S=0),"G1F"=c(C=62, H=102, N=4, O=44, S=0),"G2F"=c(C=68, H=112, N=4, O=49, S=0),"G0-GN"=c(C=42, H=69, N=3, O=30, S=0))
  
  test.glycan<-function(x){
    if(is.null(x)){
      "none"
    } else {
      x}
  }
  
  CV<-function(x){
    if(is.null(x)){
      0
    }else{x}
  }

  
  if(w.warhead=="Yep"){
  duetmAb.form.thr<-chain1.form.thr+chain2.form.thr +chain3.form.thr +chain4.form.thr + GLYCAN[[test.glycan(Glycan1)]] + GLYCAN[[test.glycan(Glycan2)]] + GLYCAN[[test.glycan(Glycan3)]] + GLYCAN[[test.glycan(Glycan4)]] + c(C=0, H=as.numeric(total.DS)*(-2), N=0, O=0, S=0) + warhead*c(C=CV(Carb), H=CV(Hydro), N=CV(Nitro), O=CV(Oxy), S=CV(Sulf)) + cys*c(C=3, H=5, N=1, O=2, S=1)   ## REPORT
  } 
  if(w.warhead=="Nope"){
  
  duetmAb.form.thr<-chain1.form.thr+chain2.form.thr +chain3.form.thr +chain4.form.thr + GLYCAN[[test.glycan(Glycan1)]] + GLYCAN[[test.glycan(Glycan2)]] + GLYCAN[[test.glycan(Glycan3)]] + GLYCAN[[test.glycan(Glycan4)]] + c(C=0, H=as.numeric(total.DS)*(-2), N=0, O=0, S=0) + cys*c(C=3, H=5, N=1, O=2, S=1)   ## REPORT
  }
  
  ## compute the theoretical formula of polypeptide chain.
  chain1.poly<- if(as.character(Seq.info$SQ[1])=="") {c("C"=0,"H"=0, "N"=0, "O"=0, "S"=0)} else{unlist(ConvertPeptide(as.character(Seq.info$SQ[1]), IAA=FALSE))}
  chain2.poly<- if(as.character(Seq.info$SQ[2])=="") {c("C"=0,"H"=0, "N"=0, "O"=0, "S"=0)} else{unlist(ConvertPeptide(as.character(Seq.info$SQ[2]), IAA=FALSE))}
  chain3.poly<- if(as.character(Seq.info$SQ[3])=="") {c("C"=0,"H"=0, "N"=0, "O"=0, "S"=0)} else{unlist(ConvertPeptide(as.character(Seq.info$SQ[3]), IAA=FALSE))}
  chain4.poly<- if(as.character(Seq.info$SQ[4])=="") {c("C"=0,"H"=0, "N"=0, "O"=0, "S"=0)} else{unlist(ConvertPeptide(as.character(Seq.info$SQ[4]), IAA=FALSE))}
  duetmAb.poly<-chain1.poly+chain2.poly+chain3.poly+chain4.poly
  #===========================================================#
  Standard<-c("C"=12.01078, "H"=1.007947, "N"=14.00672, "O"=15.99943, "S"=32.0655)
  ## compute the AVERAGE mass of theoretical formula depending on PTM input.
  duetmAb.mass.thr<-format(round(sum(duetmAb.form.thr*Standard),1), nsmall = 1)  ## REPORT
  
  chain1.mass.poly<-format(round(sum(chain1.poly*Standard),1), nsmall = 1)  ## REPORT
  chain2.mass.poly<-format(round(sum(chain2.poly*Standard),1), nsmall = 1)  ## REPORT
  chain3.mass.poly<-format(round(sum(chain3.poly*Standard),1), nsmall = 1)  ## REPORT
  chain4.mass.poly<-format(round(sum(chain4.poly*Standard),1), nsmall = 1)  ## REPORT
  
  
  duetmAb.mass.poly<-format(round(sum(duetmAb.poly*Standard),1), nsmall = 1)  ## REPORT
  
  #===============================================================================#
  chain1.info<-data.frame( "C"=chain1.poly[["C"]],"H"=chain1.poly[["H"]],"N"=chain1.poly[["N"]],"O"=chain1.poly[["O"]],"S"=chain1.poly[["S"]], "Mass (Da)"=chain1.mass.poly, "Description"="Polypeptide", check.names = FALSE)
  chain2.info<-data.frame( "C"=chain2.poly[["C"]],"H"=chain2.poly[["H"]],"N"=chain2.poly[["N"]],"O"=chain2.poly[["O"]],"S"=chain2.poly[["S"]], "Mass (Da)"=chain2.mass.poly, "Description"="Polypeptide", check.names = FALSE)
  chain3.info<-data.frame( "C"=chain3.poly[["C"]],"H"=chain3.poly[["H"]],"N"=chain3.poly[["N"]],"O"=chain3.poly[["O"]],"S"=chain3.poly[["S"]], "Mass (Da)"=chain3.mass.poly, "Description"="Polypeptide", check.names = FALSE)
  chain4.info<-data.frame( "C"=chain4.poly[["C"]],"H"=chain4.poly[["H"]],"N"=chain4.poly[["N"]],"O"=chain4.poly[["O"]],"S"=chain4.poly[["S"]], "Mass (Da)"=chain4.mass.poly, "Description"="Polypeptide", check.names = FALSE)
  
  
  duetmAb.info<-data.frame( "C"=c(duetmAb.form.thr[["C"]], duetmAb.poly[["C"]]),"H"=c(duetmAb.form.thr[["H"]], duetmAb.poly[["H"]]), "N"=c(duetmAb.form.thr[["N"]], duetmAb.poly[["N"]]),"O"=c(duetmAb.form.thr[["O"]], duetmAb.poly[["O"]]),"S"=c(duetmAb.form.thr[["S"]], duetmAb.poly[["S"]]), "Mass (Da)"=c(duetmAb.mass.thr, duetmAb.mass.poly), "Description"=c("Expected","Polypeptide"), check.names = FALSE )
  Overall<-list(chain1=chain1.info, chain2=chain2.info, chain3=chain3.info, chain4=chain4.info, duetmAb=duetmAb.info)
  return(Overall)
}
