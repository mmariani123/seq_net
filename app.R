#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinyBS)
library(markdown)
library(shiny)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library("rDNAse")
require(rDNAse)
library("ShortRead") #for converting fastq ONP to fasta
#C:\Users\Mike\Desktop\shiny_apps\seqnet\www\data
library("h2o")
library("data.table")

# Define UI for application that draws a histogram
ui <- dashboardPage(skin="green",
                    
  dashboardHeader(title = 'Frietze Lab Genomics'),

  dashboardSidebar(
    width=350,
    sidebarMenu(
    actionButton("do", 
                 "Run prediction",
                 icon("refresh"),
                 style="width:300px; 
                height:50px;
                margin-top:50px;
                font-size:20px;")
    ,
    downloadButton('downloadData', 
                   'Download Results', 
                   style="width:300px; 
                height:50px; 
                font-size:20px; 
                text-align:center;
                margin-left:16px;")
    #,
    #downloadButton('downloadExample', 
    #               'Download Example Input', 
    #               style="width:300px; 
    #            height:30px; 
    #            font-size:12px; 
    #            text-align:center;
    #            margin-left:-300px;
    #               margin-top:200px;")
  )),

  dashboardBody(
    fluidRow(align="left",
             tags$h1(
               style="font-size:16px;margin:40px;width:600px",
               "Welcome to ",tags$b("SeqNet"), ".  This APP uses a deep learning model 
                - a type of artificial neural network (ANN) - to predict whether long ONP reads
              are human or viral in nature.  Copy and paste some sequences (FASTA format) below 
               and use the buttons on the sidebar to run SeqNet and download the results once complete
              (Please be patient if inputting a large number of sequences).",
               br(),
               downloadLink("downloadExample", label = "Download example input")
               ),
           #htmlOutput("contents"),
           #plotOutput("distPlot"),
           div(style="margin:40px;",
             textAreaInput(
              "caption", 
              "SeqNet", 
              "Paste Fasta Sequences Here",
              width="600px",
              height="200px",
              resize="none")
           #tags$head(tags$script(src = "message-handler.js")),
           #verbatimTextOutput("value",style="height:400px;width:400px;")
           ),
           div(
             DT::dataTableOutput(
              "value",
              width="600px",
              height="300px"
              ),
             style="margin-left:40px;"
           )
  )
))

# Define server logic required to draw a histogram
server <- function(input, output) {

rv <- reactiveValues(m=data.frame())

observeEvent(input$do, {

withProgress(message = 'Running SeqNet', value = 0, {
#session$sendCustomMessage(
#  type = 'testmessage',
#  message = 'Thank you for clicking'
#)
  
#onp_fasta<-readFASTA(file="/slipstream/home/mmariani/hhv6_ann/reads/round_2_pass_merged.fasta")

write.table(
  input$caption, 
  file="/slipstream/home/mmariani/ShinyApps/seq_net/www/data/test.fasta",
  col.names=FALSE,
  row.names=FALSE,
  append = FALSE, 
  quote = FALSE, 
  sep = "\n")

fasta_seq <- readFASTA(file="/slipstream/home/mmariani/ShinyApps/seq_net/www/data/test.fasta")
#fasta_seq<-""

print("in")

x1<-t(sapply(fasta_seq,kmer)) #- Basic kmer and Reverse compliment kmer
incProgress(1/12, detail = paste("Finding features ", 1))
#x2<-t(sapply(fasta_seq,make_idkmer_vec)) #- Increment of diversity (ID)
#Autocorrelation
x3<-t(sapply(fasta_seq,extrDAC)) #- Dinucleotide-based auto covariance
incProgress(3/12, detail = paste("Finding features ", 3))
x4<-t(sapply(fasta_seq,extrDCC)) #- Dinucleotide-based cross covariance
incProgress(4/12, detail = paste("Finding features ", 4))
x5<-t(sapply(fasta_seq,extrDACC)) #- Dinucleotide-based auto-cross covariance
incProgress(5/12, detail = paste("Finding features ", 5))
x6<-t(sapply(fasta_seq,extrTAC)) #- Trinucleotide-based auto covariance
incProgress(6/12, detail = paste("Finding features ", 6))
x7<-t(sapply(fasta_seq,extrTCC)) #- Trinucleotide-based cross covariance
incProgress(7/12, detail = paste("Finding features ", 7))
x8<-t(sapply(fasta_seq,extrDACC)) #- Trinucleotide-based auto-cross covariance
incProgress(8/12, detail = paste("Finding features ", 8))
#Pseudo nucleotide composition
x9<-t(sapply(fasta_seq,extrPseDNC)) #- Pseudo dinucleotide composition
incProgress(9/12, detail = paste("Finding features ", 9))
x10<-t(sapply(fasta_seq,extrPseKNC)) #- Pseudo k-tupler nucleotide composition
incProgress(10/12, detail = paste("Finding features ", 10))
#x12<-t(sapply(fasta_seq,twoSeqSim))
incProgress(12/12, detail = paste("Finding features ", 12))
#x13<-t(sapply(fasta_seq,twoGOSim)) 

test_frame<-data.frame(unlist(cbind(
  x1,
  x3,
  x4,
  x5,
  x6,
  x7,
  x8,
  x9,
  x10
)))

test_frame$category=rep("virus",times=nrow(test_frame))

#There appears to be at least 1 duplicate column, so let's remove the duplicates:
test_frame <- test_frame[, !duplicated(colnames(test_frame))]
h2oServer <- h2o.init(ip="localhost", port=54321, max_mem_size="16g", nthreads=-1)
test_hex <- as.h2o(test_frame)
test_hex[,133] <- as.factor(test_hex[,133])
#model_path<-"/slipstream/home/mmariani/hhv6_ann/output/DeepLearning_model_R_1535305455367_4"

incProgress(12/12, detail = paste("Predicting sequences ... "))

model_path<-"/slipstream/home/mmariani/ShinyApps/seq_net/www/data/DeepLearning_model_R_1535305455367_4"
saved_model <- h2o.loadModel(model_path)
predictions <- h2o.predict(saved_model, test_hex)
pred_frame<-as.data.frame(predictions)
pred_frame$read<-rownames(test_frame)
pred_frame<-pred_frame[,c(4,2,3,1)]
colnames(pred_frame)<-c("read","p(human)", "p(virus)", "predict")

rv$m<-pred_frame

output$value <- DT::renderDataTable(
  {setDT(pred_frame,keep.rownames=FALSE)},
  options = list(scrollX = TRUE, scrollY=TRUE)
  )

#write.csv(as.data.frame(pred_frame), 
#          file="C:\\Users\\Mike\\Desktop\\shiny_apps\\seqnet\\www\\data\\shiny_test_predictions.csv",
#          row.names=FALSE)

})
  
})

output$downloadData <- downloadHandler(
  filename = function() {
    #paste(input$dataset, ".csv", sep = "")
    paste("virus_predictions", ".csv", sep = "")
  },
  content = function(file) {
    write.csv(rv$m, file, row.names = FALSE)
  }
)

output$downloadExample <- downloadHandler(
  filename = function() {
    paste("seqnet_example_input",".txt",sep="")
  },
  content = function(file) {
    file.copy(
      "/slipstream/home/mmariani/hhv6_ann/reads/seqnet_example_input.txt",
      file)
    contentType = "text/plain"
  }
)

#output$downloadData <- downloadHandler(
#  filename <- function() {
#    paste("output", "zip", sep=".")
#  },
#  
#  content <- function(file) {
#    file.copy("out.zip", file)
#  },
#  contentType = "application/zip"
#)

}

# Run the application 
shinyApp(ui = ui, server = server)

