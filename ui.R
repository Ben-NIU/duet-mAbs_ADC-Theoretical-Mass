library(shinydashboard)

dashboardPage(skin="red",
  dashboardHeader(title="duet-mAbs/ADC Intact Mass Calculator", titleWidth = 400),
  dashboardSidebar(
    sidebarMenu(
      fileInput("input.file", label=span("Input FASTA here", style="font-family:'calibri'; font-size:12pt")),
      hr(),
      menuItem("Polypeptide", icon=icon("paw"), tabName="poly"),
      menuItem("non-Reduced MS", icon=icon("sitemap"), tabName = "nonRed"),
      menuItem("Reduced MS", icon=icon("cubes"), tabName = "Red"),
      hr(),
      span(textOutput("CH1"), style="font-family:'calibri';font-size: 14pt;color:#FFFF00; text-decoration: underline"),
      uiOutput("nterm.ch1"),
      uiOutput("cterm.ch1"),
      span(textOutput("CH2"), style="font-family:'calibri';font-size: 14pt;color:#FFFF00; text-decoration: underline"),
      uiOutput("nterm.ch2"),
      uiOutput("cterm.ch2"),     
      span(textOutput("CH3"), style="font-family:'calibri';font-size: 14pt;color:#FFFF00; text-decoration: underline"),
      uiOutput("nterm.ch3"),
      uiOutput("cterm.ch3"), 
      span(textOutput("CH4"), style="font-family:'calibri';font-size: 14pt;color:#FFFF00; text-decoration: underline"),
      uiOutput("nterm.ch4"),
      uiOutput("cterm.ch4"),
      hr(),
      radioButtons("with.warhead", label="Contains warhead?", choices = list("Yep","Nope"), selected = "Nope", inline = TRUE),
      
      uiOutput("warhead.f")
      )
    ),
    
  
  dashboardBody(
    tags$style(HTML("

                              
                              .box.box-solid.box-primary>.box-header {
                              color:#fff;
                              background:#CDAA7D
                              }
                              
                              .box.box-solid.box-primary{
                              border-bottom-color:##CDAA7D;
                              border-left-color:##CDAA7D;
                              border-right-color:##CDAA7D;
                              border-top-color:##CDAA7D;
                              }
                              
                              ")),
    
    tabItems(
      tabItem(tabName = "poly",
      
    fluidRow(
   
      box(title="Show Sequences", status="primary", solidHeader = TRUE, collapsible = TRUE, width=12,
          fluidRow(
            column(width=9,
            p(verbatimTextOutput("seq1"))),
            column(width=3,
            fluidRow(
              br(),
            valueBoxOutput("seq1.count", width=12)))),
          fluidRow(
            column(width=9,
            p(verbatimTextOutput("seq2"))),
            column(width=3,
             fluidRow(
               br(),
            valueBoxOutput("seq2.count", width=12)))),
          fluidRow(
            column(width=9,
            p(verbatimTextOutput("seq3"))),
            column(width=3,
               fluidRow(
                br(),
              valueBoxOutput("seq3.count", width=12)))),
          fluidRow(
            column(width=9,
              p(verbatimTextOutput("seq4"))),
            column(width=3,
                fluidRow(
                 br(),
                valueBoxOutput("seq4.count", width=12))))
          
      )),
  
    
    fluidRow(
      box(title="Polypeptide Stats", status="primary", solidHeader = TRUE, collapsible = TRUE, width=9,
         
          fluidRow(
            column(width=12,
            verbatimTextOutput("poly.note")
                   )
          ),
          
          fluidRow(
           column(width=4,
            tableOutput("tb1")),
           column(width=4, 
          valueBoxOutput("poly.mass.chain1", width=10))),
         hr(),
         fluidRow(
           column(width=4,
             tableOutput("tb2")),
           column(width=4,
          valueBoxOutput("poly.mass.chain2", width=10))),
         hr(),
         fluidRow(
           column(width=4,
                  tableOutput("tb3")),
           column(width=4,
                  valueBoxOutput("poly.mass.chain3", width=10))),
         hr(),
         fluidRow(
           column(width=4,
                  tableOutput("tb4")),
           column(width=4,
                  valueBoxOutput("poly.mass.chain4", width=10))), 
         hr(),
         fluidRow(
           column(width=4,
              tableOutput("tb5")),
           column(width=4,
            valueBoxOutput("poly.mass.duetmAb", width=10)))
      ))
      ),
    
    tabItem(tabName = "nonRed",
    
    fluidRow(          
      column(width=5,
      box(title="Disulfide Bonds and N-glycans", status="primary", solidHeader = TRUE, collapsible = TRUE, width = "100%",
        fluidRow(
          column(width=6,
              numericInput("tot.DS", label="# total disulfide bonds", value = 16, width="100%"))),
          fluidRow(
          column(width=6,
            uiOutput("nglycan1")),
          column(width=6,
            uiOutput("nglycan2"))),
          fluidRow(
            column(width=6,
              uiOutput("nglycan3")),
            column(width=6,
              uiOutput("nglycan4")))
          
          
        )),
      column(width=7,
      box(title="Formula & Mass", status="primary", solidHeader = TRUE, collapsible = TRUE,width="100%",
          fluidRow(
            column(width=12,
            tableOutput("tb22"))),
          fluidRow(
            column(width=12,
            valueBoxOutput("NR.mass")))
      ))),
    fluidRow(
      column(width=5,
      box(title="ADC Modifications", status="primary", solidHeader = TRUE, collapsible = TRUE, width="100%",
          fluidRow(
          column(width=12,  
          sliderInput("Cys", "Cysteinylation #", min=0, max=10, value=2, step=1))),
          br(),
          fluidRow(
          column(width=12,  
          sliderInput("Warhead", "How many warheads in total?", min=0, max=8, value=2, step=1))
          
         )
         
          ) 
      )
    )
 
    ),
    
    tabItem(tabName = "Red",
      fluidRow(      
      column(width=4,
             
       box(title="# UnReduced intra-chain Disuldife bonds", status="primary", solidHeader = TRUE, collapsible = TRUE,width="100%",
          fluidRow(
            column(width=12,
                 numericInput("unR.ch1", "chain1", value=3, min=0, max=10, step=1))),
          fluidRow(
            column(width=12,
      
                 numericInput("unR.ch2", "chain2", value=3, min=0, max=10, step=1))),
          fluidRow(
            column(width=12,
                   
                numericInput("unR.ch3", "chain3", value=3, min=0, max=10, step=1))),       
          fluidRow(
            column(width=12,
                   
                numericInput("unR.ch4", "chain4", value=3, min=0, max=10, step=1)))),       
          
       
          box(title="N-glycan Identity", status="primary", solidHeader = TRUE, collapsible = TRUE, width = "100%",
            
              column(width=12,
                  span("Note: The selected N-glycan will be applied to all 4 chains.", style="font-family:'calibri'; font-size: 12pt; color:#CC6600")),
              column(width=6,
     
                 selectInput("glyco.Red", label=NULL, choices = list("none","Deglycosylated",'G0', "G0F", "G1", "G1F", "G2F","G0-GN"), selected = "G0"))),
          box(title="Locate Warheads", status="primary", solidHeader = TRUE, collapsible = TRUE, width="100%",
              column(width=12,
                    checkboxGroupInput("which.chain", "Select the chains that carry the warheads", choices = list("Chain1"=1, "Chain2"=2, "Chain3"=3, "Chain4"=4), selected=c(2,4)) 
                     )
             ),
         box(title="# of Warheads for each chain", status="primary", solidHeader = TRUE, collapsible = TRUE, width="100%",
             column(width=12,
                    span("Note: Specify the # of warheads for each chain below", style="font-family:'calibri'; font-size: 12pt;color:#CC6600")),
             column(width=8,
                    uiOutput("#Ch1"),
                    uiOutput("#Ch2"),
                    uiOutput("#Ch3"),
                    uiOutput("#Ch4")
                    ))
       
       
       ),
      
      column(width=6,
       box(title="Formula & Mass", status="primary", solidHeader = TRUE, collapsible = TRUE, width="100%",
            fluidRow(
              column(width=12,
              tableOutput("tbch1"))
              ),       
            fluidRow(
              column(width=12,
              tableOutput("tbch2"))
              ),
           fluidRow(
             column(width=12,
              tableOutput("tbch3"))
           ),       
           fluidRow(
             column(width=12,
              tableOutput("tbch4"))
           ),        
           
           
            fluidRow(
              column(width=12,
             verbatimTextOutput("note")))
             
          )) )       
          
    )
  
    
    
  )
)
)
          
          
          
          
          