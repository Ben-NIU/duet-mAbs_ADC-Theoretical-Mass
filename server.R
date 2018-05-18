source("read.duet.fasta.R")
source("show.table.nonR.R")
source("show.table.Red.R")
library(OrgMassSpecR)

shinyServer(function(input, output) {
  
  PRO<-reactive({ 
    validate(
      need(input$input.file$datapath !="", "")
    )
    read.duet.fasta(input$input.file[['datapath']])
     })
  
  observeEvent(input$input.file,{
   output$seq1<-renderText({paste(as.character(PRO()$Name[1]),"\n",as.character(PRO()$SQ[1]), sep="")})
  })
  
  observeEvent(input$input.file,{
    output$seq2<-renderText({paste(as.character(PRO()$Name[2]),"\n",as.character(PRO()$SQ[2]), sep="")})
  })
  
  observeEvent(input$input.file,{
    output$seq3<-renderText({paste(as.character(PRO()$Name[3]),"\n",as.character(PRO()$SQ[3]), sep="")})
  })

  observeEvent(input$input.file,{
    output$seq4<-renderText({paste(as.character(PRO()$Name[4]),"\n",as.character(PRO()$SQ[4]), sep="")})
  })
  
  observeEvent(input$input.file, {
     output$seq1.count<-renderValueBox({
       valueBox(value=tags$p(paste(as.character(PRO()$LTH[1]), "  amino acids", sep=""), style="font-size: 50%"), subtitle = tags$p("Chain 1", style="font-size: 120%"), icon=icon("angle-double-left"), color = "green")})
   })
  
   
   observeEvent(input$input.file, {
     output$seq2.count<-renderValueBox({
       valueBox( value=tags$p(paste(as.character(PRO()$LTH[2]), "  amino acids", sep=""),style="font-size: 50%"), subtitle = tags$p("Chain 2", style="font-size: 120%"), icon=icon("angle-left"), color = "light-blue")})
   })
   
   observeEvent(input$input.file, {
     output$seq3.count<-renderValueBox({
       valueBox( value=tags$p(paste(as.character(PRO()$LTH[3]), "  amino acids", sep=""),style="font-size: 50%"), subtitle = tags$p("Chain 3", style="font-size: 120%"), icon=icon("angle-left"), color = "teal")})
   })   
   
   observeEvent(input$input.file, {
     output$seq4.count<-renderValueBox({
       valueBox( value=tags$p(paste(as.character(PRO()$LTH[4]), "  amino acids", sep=""),style="font-size: 50%"), subtitle = tags$p("Chain 4", style="font-size: 120%"), icon=icon("angle-left"), color = "yellow")})
   })  
   
  
   output$warhead.f<-renderUI({
     if(input$with.warhead=="Yep" ){
       fluidRow(
         column(width=6,
                numericInput("carbon#","C", value=0, min=0, max=500, step=1),
                numericInput("nitrogen#", "N", value=0, min=0, max=500, step=1),
                numericInput("sulfer#", "S", value=0, min=0, max=100, step=1)),
         column(width=6,
                numericInput("hydrogen#","H", value=0, min=0, max=500, step=1),
                numericInput("oxygen#", "O", value=0, min=0, max=500, step=1))
       )
       
       
         } else {
       NULL}
   })
   
   
   
   output$CH1<-renderText({
     if((substr(PRO()$SQ[1],1,1) %in% c("Q","E")) |(substr(PRO()$SQ[1], PRO()$LTH[1]-5, PRO()$LTH[1]) =="SLSPGK") ){
       paste(">" ,"Chain 1",sep=" ")     } else {
       NULL}
       
   })
 
   
  output$nterm.ch1<-renderUI({
    if(substr(PRO()$SQ[1],1,1) %in% c("Q","E")){
          radioButtons("nterm.ch1", label="N-terminal PyroE", choices=list("Yes","No"), selected = "Yes", inline=TRUE)
    } else {
      NULL}
  })
  
  output$cterm.ch1<-renderUI({
    if(substr(PRO()$SQ[1], PRO()$LTH[1]-5, PRO()$LTH[1]) =="SLSPGK") {
      radioButtons("cterm.ch1", label="C-terminal Lys Loss", choices=list("Yes","No"), selected="Yes", inline=TRUE)
    } else{
      NULL}
  })
  
  
  output$CH2<-renderText({
    if((substr(PRO()$SQ[2],1,1) %in% c("Q","E")) |(substr(PRO()$SQ[2], PRO()$LTH[2]-5, PRO()$LTH[2]) =="SLSPGK") ){
      paste(">" ,"Chain 2",sep=" ")     } else {
        NULL}
    
  })
  
  output$nterm.ch2<-renderUI({
    if(substr(PRO()$SQ[2],1,1) %in% c("Q","E")){
      radioButtons("nterm.ch2", label="N-terminal PyroE", choices=list("Yes","No"), selected = "Yes", inline=TRUE)
    } else {
      NULL}
  })
  
  output$cterm.ch2<-renderUI({
    if(substr(PRO()$SQ[2], PRO()$LTH[2]-5, PRO()$LTH[2]) =="SLSPGK") {
      radioButtons("cterm.ch2", label="C-terminal Lys Loss", choices=list("Yes","No"), selected="Yes", inline=TRUE)
    } else{
      NULL}
  })
  
  output$CH3<-renderText({
    if((substr(PRO()$SQ[3],1,1) %in% c("Q","E")) |(substr(PRO()$SQ[3], PRO()$LTH[3]-5, PRO()$LTH[3]) =="SLSPGK") ){
      paste(">" ,"Chain 3",sep=" ")     } else {
        NULL}
    
  })
  
  output$nterm.ch3<-renderUI({
    if(substr(PRO()$SQ[3],1,1) %in% c("Q","E")){
      radioButtons("nterm.ch3", label="N-terminal PyroE", choices=list("Yes","No"), selected = "Yes", inline=TRUE)
    } else {
      NULL}
  })
  
  output$cterm.ch3<-renderUI({
    if(substr(PRO()$SQ[3], PRO()$LTH[3]-5, PRO()$LTH[3]) =="SLSPGK") {
      radioButtons("cterm.ch3", label="C-terminal Lys Loss", choices=list("Yes","No"), selected="Yes", inline=TRUE)
    } else{
      NULL}
  })
  
  output$CH4<-renderText({
    if((substr(PRO()$SQ[4],1,1) %in% c("Q","E")) |(substr(PRO()$SQ[4], PRO()$LTH[4]-5, PRO()$LTH[4]) =="SLSPGK") ){
      paste(">" ,"Chain 4",sep=" ")     } else {
        NULL}
    
  })
  output$nterm.ch4<-renderUI({
    if(substr(PRO()$SQ[4],1,1) %in% c("Q","E")){
      radioButtons("nterm.ch4", label="N-terminal PyroE", choices=list("Yes","No"), selected = "Yes", inline=TRUE)
    } else {
      NULL}
  })
  
  output$cterm.ch4<-renderUI({
    if(substr(PRO()$SQ[4], PRO()$LTH[4]-5, PRO()$LTH[4]) =="SLSPGK") {
      radioButtons("cterm.ch4", label="C-terminal Lys Loss", choices=list("Yes","No"), selected="Yes", inline=TRUE)
    } else{
      NULL}
  })
      
  
  observeEvent(input$input.file,{
    output$poly.note<-renderText('The masses and elemental compositions shown below correspond to each polypeptide chain (as indicated in each box), with no consideration of PTMs and disulfide bonds.')
  })
  
  
  
  
## ===================================== non-reduced, N-glycan UI==================================== 
  output$nglycan1<-renderUI({
    if(PRO()$SQ[1]!="") {
      selectInput("nglycan1", label="N-Glycan of Chain 1", choices = list("none","Deglycosylated",'G0', "G0F", "G1", "G1F","G2F", "G0-GN"), selected = "none", width="100%")
    } else{
      NULL}
  })
  output$nglycan2<-renderUI({
    if(PRO()$SQ[2]!="") {
      selectInput("nglycan2", label="N-Glycan of Chain 2", choices = list("none","Deglycosylated",'G0', "G0F", "G1", "G1F","G2F", "G0-GN"), selected = "none", width="100%")
    } else{
      NULL}
  })
  output$nglycan3<-renderUI({
    if(PRO()$SQ[3]!="") {
      selectInput("nglycan3", label="N-Glycan of Chain 3", choices = list("none","Deglycosylated",'G0', "G0F", "G1", "G1F","G2F", "G0-GN"), selected = "none", width="100%")
    } else{
      NULL}
  })
  output$nglycan4<-renderUI({
    if(PRO()$SQ[4]!="") {
      selectInput("nglycan4", label="N-Glycan of Chain 4", choices = list("none","Deglycosylated",'G0', "G0F", "G1", "G1F","G2F", "G0-GN"), selected = "none", width="100%")
    } else{
      NULL}
  })
  
  output$'#Ch1'<-renderUI({
    if(1 %in% input$which.chain) {
      numericInput("#Ch1", label="Chain 1", value=1, min=0, max=8, step=1, width="100%")
    } else{
      NULL}
  })
  
  output$'#Ch2'<-renderUI({
    if(2 %in% input$which.chain) {
      numericInput("#Ch2", label="Chain 2", value=1, min=0, max=8, step=1, width="100%")
    } else{
      NULL}
  })
  
  output$'#Ch3'<-renderUI({
    if(3 %in% input$which.chain) {
      numericInput("#Ch3", label="Chain 3", value=1, min=0, max=8, step=1, width="100%")
    } else{
      NULL}
  })
  
  output$'#Ch4'<-renderUI({
    if(4 %in% input$which.chain) {
      numericInput("#Ch4", label="Chain 4", value=1, min=0, max=8, step=1, width="100%")
    } else{
      NULL}
  })
  
  
  Tables<-reactive({ show.table.nonR(PRO(), input$nterm.ch1, input$cterm.ch1, input$nterm.ch2, input$cterm.ch2,input$nterm.ch3, input$cterm.ch3, input$nterm.ch4, input$cterm.ch4, input$nglycan1, input$nglycan2, input$nglycan3, input$nglycan4, input$tot.DS, input$Cys, input$Warhead,input$'carbon#', input$'hydrogen#', input$'nitrogen#', input$'oxygen#', input$'sulfer#', input$with.warhead)}) 
  
  
  Tables.Red<-reactive({ show.table.Red(PRO(), input$nterm.ch1, input$cterm.ch1, input$nterm.ch2, input$cterm.ch2,input$nterm.ch3, input$cterm.ch3, input$nterm.ch4, input$cterm.ch4, input$glyco.Red, input$unR.ch1, input$unR.ch2, input$unR.ch3, input$unR.ch4, input$which.chain, input$'carbon#', input$'hydrogen#', input$'nitrogen#', input$'oxygen#', input$'sulfer#', input$'#Ch1',input$'#Ch2', input$'#Ch3', input$'#Ch4', input$with.warhead) })
  
  observeEvent(input$input.file, {
    output$poly.mass.chain1<-renderValueBox({
      valueBox( value=tags$p(paste(Tables()$chain1$'Mass (Da)', "  Da", sep=""), style="font-size: 70%"), subtitle = tags$p("Chain1", style="font-size: 130%"), icon=icon("star-half-empty"), color = "green")})
  })

  observeEvent(input$input.file,{
    output$poly.mass.chain2<-renderValueBox({
      valueBox( value=tags$p(paste(Tables()$chain2$'Mass (Da)', "  Da", sep=""), style="font-size: 70%"), subtitle = tags$p("Chain2", style="font-size: 130%"), icon=icon("star-half-o"), color = "light-blue")})
  })
  
  observeEvent(input$input.file,{
    output$poly.mass.chain3<-renderValueBox({
      valueBox( value=tags$p(paste(Tables()$chain3$'Mass (Da)', "  Da", sep=""), style="font-size: 70%"), subtitle = tags$p("Chain3", style="font-size: 130%"), icon=icon("star-half-o"), color = "teal")})
  })
  
  observeEvent(input$input.file,{
    output$poly.mass.chain4<-renderValueBox({
      valueBox( value=tags$p(paste(Tables()$chain4$'Mass (Da)', "  Da", sep=""), style="font-size: 70%"), subtitle = tags$p("Chain4", style="font-size: 130%"), icon=icon("star-half-o"), color = "yellow")})
  })
  
  observeEvent(input$input.file,{
    output$poly.mass.duetmAb<-renderValueBox({
      valueBox( value=tags$p(paste(Tables()$duetmAb$'Mass (Da)'[2], "  Da", sep=""), style="font-size: 70%"), subtitle = tags$p("duet-mAb", style="font-size: 130%"), icon=icon("paper-plane-o"), color = "purple")})
  })
  
  output$tb1<-renderTable({
    Tables()$chain1[,1:5]}, caption=NULL,digits=0, align="c")
  
  output$tb2<-renderTable({
    Tables()$chain2[,1:5]}, caption=NULL,digits=0, align="c")
  
  output$tb3<-renderTable({
    Tables()$chain3[,1:5]}, caption=NULL,digits=0, align="c")
  
  output$tb4<-renderTable({
    Tables()$chain4[,1:5]}, caption=NULL,digits=0, align="c")
  
  output$tb5<-renderTable({
    Tables()$duetmAb[,1:5][2,]}, caption=NULL,digits=0, align="c")

  output$tb22<-renderTable({
    Tables()$duetmAb[,1:5][1,]}, caption=NULL,digits=0, align="c")
  
  
  
  observeEvent(input$input.file,{
    output$NR.mass<-renderValueBox({
      valueBox( value=tags$p(paste(Tables()$duetmAb$'Mass (Da)'[1], "  Da", sep=""), style="font-size: 70%"), subtitle = tags$p("duet-mAb", style="font-size: 130%"), icon=icon("paper-plane-o"), color = "purple")})
  })
  
  output$tbch1<-renderTable({
    Tables.Red()$chain1}, caption="Chain1", caption.placement = getOption("xtable.caption.placement", "top"), digits = 1,align="c")
  output$tbch2<-renderTable({
    Tables.Red()$chain2}, caption="Chain2", caption.placement = getOption("xtable.caption.placement", "top"), digits = 1, align="c")
  output$tbch3<-renderTable({
    Tables.Red()$chain3}, caption="Chain3", caption.placement = getOption("xtable.caption.placement", "top"), digits = 1,align="c")
  output$tbch4<-renderTable({
    Tables.Red()$chain4}, caption="Chain4", caption.placement = getOption("xtable.caption.placement", "top"), digits = 1, align="c")
  
  
  observeEvent(input$input.file,{
    output$note<-renderText('The "Expected" compositions/masses above, by default, correspond to chains with 3 un-reduced intra-disulfide bonds. Change the number of un-reduced intra-disulfide bonds as needed.')
  })
    
  
}
)
