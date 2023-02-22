#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(Familias)
library(Rmpfr)

# data("NorwegianFrequencies")



# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("ShinyPedigree"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("strResultFile", "STR results"),
      selectInput("dadID", "Father", ""),
      selectInput("momID", "Mother", ""),
      selectInput("childID", "Child", ""),
      tableOutput("inputTbl")
      # ,
      # checkboxInput("passUnknownAlleles",
      #               "Pass Unkonown Alleles",
      #               value = TRUE)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Result", 
                           h3("Excluded Alleles"),
                           tableOutput("excludedAlleles"),
                           h3("Result"),
                           htmlOutput("resultTxt"),
                           h3("RAW Result"),
                           verbatimTextOutput("rawResultTxt")),
                  # tabPanel("Input Table", 
                  # tableOutput("inputTbl")
                  # ),
                  tabPanel("Loci Table",
                           verbatimTextOutput("lociDataOutput"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  lociData <- reactive({
    # lociTableRAW <- read_tsv("cordis plus.csv")
    # loci <- colnames(lociTableRAW)[seq(1, ncol(lociTableRAW), by = 2)]
    # lociData <- lapply(seq(1, ncol(lociTableRAW), by = 2), function(locusI) {
    #   freqs <- as.numeric(lociTableRAW[[ locusI + 1]][-1])
    #   names(freqs) <- lociTableRAW[[ locusI]][-1]
    #   freqs <- freqs[!is.na(freqs)]
    #   nZero <- length(freqs[freqs == 0])
    #   freqs[freqs == 0] <- 0.0001
    #   maxFreq <- which(freqs == max(freqs))[1]
    #   freqs[maxFreq] <- freqs[maxFreq] - 0.0001 * nZero
    #   freqs
    # })
    # names(lociData) <- loci
    frinp <- readxl::read_excel("all_freq.xlsx")
    nColumns <- ncol(frinp)
    Alleles <- colnames(frinp)[7:nColumns]
    lociData <- map(
      1:nrow(frinp),
      \(locusI) {
        locusFreq <- unlist(frinp[locusI, 7:nColumns])
        locusFreq <- locusFreq[!is.na(locusFreq)]
        sumFreq <- sum(locusFreq)
        if (sumFreq < 1)
          locusFreq["correctionTo_1"] <- 1 - sumFreq
        else if (sumFreq > 1)
          locusFreq[1] <- locusFreq[1] - (sumFreq - 1)
        locusFreq
      })
    names(lociData) <- frinp$Locus
    lociData
  })
  
  strResult <- reactive({
    strFile <- input$strResultFile
    if (is.null(strFile)) {
      return()
    }
    strResult <- readODS::read_ods(strFile$datapath,
                                   col_names = FALSE,
                                   col_types = NA)
    colnames(strResult) <- c("smpl", "marker", "allele1", "allele2")
    strResult %>%
      pivot_longer(cols = c("allele1", "allele2"), names_to = NULL,
                   values_to = "allele") %>%
      mutate(allele = str_extract(allele, "[0-9XY]+\\.?[1-9]?+")) %>% 
      # filter(marker != "AMEL") %>% 
      group_by(marker) %>% 
      mutate(excludeUnknownAllele = !(allele %in% names(lociData()[[marker[1]]])),
             excludeMarker = any(excludeUnknownAllele))
    
  })
  
  output$excludedAlleles <- renderTable({
    if (is.null(strResult())) return()
    strResult() %>% 
      filter(excludeUnknownAllele & marker != "AMEL") %>% 
      select(-c(excludeUnknownAllele, excludeMarker))
  })
  
  output$inputTbl <- renderTable({
    if (is.null(strResult())) return()
    strResult() %>% 
      select(-c(excludeUnknownAllele, excludeMarker)) %>%
      group_by(smpl, marker) %>% 
      summarise(allele = paste(allele, collapse = " / ")) %>% 
      pivot_wider(names_from = smpl, values_from = allele)
  })
  
  observe({
    if (is.null(strResult())) return()
    isolate({
      smplNames <- c("", unique(strResult()$smpl))
      updateSelectInput(session, "dadID",
                        choices = smplNames)
      updateSelectInput(session, "momID",
                        choices = smplNames)
      updateSelectInput(session, "childID",
                        choices = smplNames)
      if (length(smplNames) == 3) {
        updateSelectInput(session, "dadID", 
                          selected = smplNames[2])
        updateSelectInput(session, "childID", 
                          selected = smplNames[3])
      }
      else if (length(smplNames) == 4) {
        updateSelectInput(session, "dadID", 
                          selected = smplNames[2])
        if ("AMEL" %in% strResult()$marker && 
            "Y" %in% (strResult() %>% 
                      filter(smpl == smplNames[3] & marker == "AMEL") %>% 
                      pull(allele))
        ) { 
          updateSelectInput(session, "momID", 
                            selected = smplNames[4])
          updateSelectInput(session, "childID", 
                            selected = smplNames[3])
        } else {
          updateSelectInput(session, "momID", 
                            selected = smplNames[3])
          updateSelectInput(session, "childID", 
                            selected = smplNames[4])
        }
      }
    })
  })
  
  resultRelationship <- reactive({
    if (input$dadID == "" || input$childID == "") 
      return()
    # if (any(strResult()[["excludeUnknownAllele"]] & !input$passUnknownAlleles))
    #   return()
    markers <- unique(strResult() %>% 
                        filter(!excludeMarker) %>%
                        pull(marker)) #%>% toupper()
    
    
    loci <- lapply(markers, function(marker)
    { 
      FamiliasLocus(
        lociData()[[marker]],
        name = marker)
    })
    
    strResult <- 
      strResult() %>% 
      filter(!excludeMarker) %>% 
      select(-c(excludeUnknownAllele, excludeMarker))
    # strResult <- strResult %>% 
    #   select(-excludeUnknownAllele)
    persons <- c(input$childID, input$dadID, 
                 {
                   if (input$momID != "") 
                     input$momID else "unkn"
                 })
    sex <- c("male", "male", "female")
    ped1 <- FamiliasPedigree(id = persons, dadid = c(NA, NA, NA),
                             momid = c(
                               {
                                 if (input$momID != "") 
                                   input$momID else NA
                               }, 
                               NA, NA), 
                             sex = sex)
    ped2 <- FamiliasPedigree(id = persons, dadid = c(input$dadID, NA, NA),
                             momid = c(
                               {
                                 if (input$momID != "") 
                                   input$momID else NA
                               },
                               NA, NA), 
                             sex = sex)
    pedigrees <- list(notFather = ped1, isFather = ped2)
    
    p1 <- strResult %>% 
      filter(smpl == input$childID) %>% 
      pull(allele)
    p2 <- strResult %>% 
      filter(smpl == input$dadID) %>% 
      pull(allele)
    p3 <- if (input$momID != "") {
      strResult %>% 
        filter(smpl == input$momID) %>% 
        pull(allele) 
    } else NULL
    datamatrix <- rbind(p1, p2, p3)
    rownames(datamatrix) <- 
      if (input$momID != "") 
        persons else persons[-3]
    
    # result <<-
    FamiliasPosterior(pedigrees, 
                      loci
                      , datamatrix)
    
    # result
  })
  
  
  
  output$resultTxt <- renderUI({
    if (is.null(resultRelationship())) return()
    LR <- as.numeric(resultRelationship()$LR["isFather"])
    calcr <- mpfr((LR / (LR + 1)) * 100, 64)
    sprintf("<b>Probability of paternity: %s %%<br>
            Likelihood Ratio (LR): %s</b>", 
            format(calcr , digits = NULL),
            LR) %>% 
      HTML()
  })
  
  output$rawResultTxt <- renderPrint({
    if (is.null(resultRelationship())) return()
    print(resultRelationship())
  })
  
  output$lociDataOutput <- renderPrint({
    if (is.null(lociData())) return()
    print(lociData())
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
