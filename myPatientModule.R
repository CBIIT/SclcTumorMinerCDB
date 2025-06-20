##UI CODE -------------------------------------------------------------------------
myPatientModuleUI <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "MyPatient",
    fluidPage(
      loadingModal(), ## new 04 17 2024
      title = "MyPatient",
      tabPanel(
        "MyPatient",
        fluidPage(
          useShinyjs(), 
          tags$head(
            tags$style(
              HTML("
                .hidden-tab {
                  display: none;
                }
              ")
            )
          ),
         
          tabsetPanel(
            id = ns("innerTabs1"),
            # tabPanel(
            #   'Find a patient',
            #   mainPanel(
            #     h3("Please Enter patient ID"), 
            #     width = 9,
            #     textInput("patID", "Patient Identifier:", ""),
            #     actionButton(ns("FindMyPatient"), "Find"),
            #   )
            # ),
            tabPanel(
              'Patient Look Up',
              mainPanel(
               ##  h3("Please select a patient by clicking on a row in the table below to view genomic data"), 
                h3("Please type patient Identifier:"), 
                ## selectInput(ns("patID"), "",choices = unique(srcContent$patient$sampleData$patientID), multiple = F, selected = ""),
                ## textInput(ns("patID"), "", ""),
                uiOutput(ns("patIDUi")),
                ## actionButton(ns("FindMyPatient"), "Find"),
                width = 9,
                h3("or select a patient by clicking on a row in the table below"),
                
                ## new hide table ....032425..............................
                
                checkboxInput(ns("showTable"), "Show Table", value=FALSE),
                uiOutput(ns("patTableUi"))
                ## actionButton(ns("MyPatientTable"), "Show patient table"),
                
                ## end new .................................................
                
                # DT::dataTableOutput(ns("myPatientSearchData"))
              )
            ),
            tabPanel(
              "Clinical Data",
              fluidRow(
                column(
                  width = 12,
                  # tags$div(
                  #   style = "font-size: 25px; color: #0000FF;", 
                  #   textOutput(ns("message"))
                  # ),
                 
                  ## h3("Demographic Information"),
                  uiOutput(ns("patient_demographic_info")),
                  br()
                ),
                column(
                  width = 12,
                  #h3("Sample Data"),
                  tags$p(style = "font-size: 24px; color: blue;","Sample Data"),
                  tags$div(
                  style = "font-size: 21px; color: #333333; border: 2px solid black;",
                  DT::dataTableOutput(ns("patient_sample_data"))
                  ), 
                  h4("Please select a row on the Sample Data table and select the Find Similar Patient or Cell line Button to find similar samples"), 
                ),
                column(
                  width = 12,
                  style = "display: flex; justify-content: flex-start; align-items: center;",
                  br(),
                    disabled(actionButton(ns("find_similar_patient"), "Find Similar Patient or Cell Line", width = "25%", class = "btn-primary btn-lg", style = "background-color: #e3e3e3; color: black; margin-right: 10px;")),
                   ##  disabled(actionButton(ns("find_similar_cell_line"), "Find Similar Cell Line", width = "15%", class = "btn-primary btn-lg", style = "background-color: #e3e3e3; color: black;"))
                )
              ),class = "hidden-tab"
            ),
            tabPanel(
              "Gene Expression",
              fluidRow(
                sidebarLayout(
                  sidebarPanel(
                    width = 2,
                    uiOutput(ns("gene_expression_info")),
                    uiOutput(ns("gene_filter_label")),
                    br(),
                    radioButtons(
                      ns("gene_filter"),
                      "Gene Filter:",
                      choices = c("All Genes", "Protein Coding", "Biomarkers"),
                      selected = "Protein Coding"
                    ),
                    br(),
                    radioButtons(
                      ns("slider_type"),
                      "Slider Type:",
                      ## choices = c("Gene Expression Range For At Least 1 Sample", "Top and Bottom Genes"),
                      choices = c("Gene Expression Range", "Top and Bottom Genes"),
                      selected = "Gene Expression Range"
                    ),
                    br(),
                    uiOutput(ns("gene_slider"))
                  ),
                  mainPanel(
                    column(
                      12,
                      br(), 
                      # selectInput(ns("pageLengthInputGene"), "Entries Displayed: ", choices = c(10, 15, 25), selected = 10, width = "165px"),
                      h3("Gene Expression Data"),
                      HTML(paste0("<b>","Note:","</b>", " The genes values are in log2(FPKM + 1).")),
                      br(),br(),
                      DT::dataTableOutput(ns("gene_expression_data"))
                    )
                  )
                )
              ), class = "hidden-tab"
            ),
            tabPanel(
              "Mutations",
              fluidRow(
                sidebarLayout(
                  sidebarPanel(
                    width = 2,
                    uiOutput(ns("mutation_info")),
                    uiOutput(ns("mutation_filter_label")),
                    br(),
                    radioButtons(
                      ns("mutation_filter"),
                      "Mutation Filter:",
                      choices = c("All Mutations", "Biomarkers"),
                      selected = "All Mutations"
                    ),
                    br(),
                    radioButtons(
                      ns("slider_type_mutation"),
                      "Slider Type:",
                      choices = c("Gene Mutation Range For At Least 1 Sample", "Top and Bottom Mutations"),
                      selected = "Top and Bottom Mutations"
                    ),
                    br(),
                    uiOutput(ns("mutation_slider"))
                  ),
                  mainPanel(
                    column(
                      12,
                      br(), 
                      # selectInput(ns("pageLengthInputMutation"), "Entries Displayed: ", choices = c(10, 15, 25), selected = 10, width = "165px"),
                      h3("Mutation Data"),
                      HTML(paste0("<b>","Note:","</b>", " The genes values are between 1 and 100 and represent the probability of homozygous, function-impacting somatic variants (see details in help section).")),
                      br(),br(),
                      DT::dataTableOutput(ns("mutation_data"))
                    )
                  )
                )
              ), class = "hidden-tab"
            ),
            ## methylation 
            
            tabPanel(
              "Methylation",
              fluidRow(
                sidebarLayout(
                  sidebarPanel(
                    width = 2,
                    uiOutput(ns("methylation_info")),
                    uiOutput(ns("methylation_filter_label")),
                    br(),
                    radioButtons(
                      ns("methylation_filter"),
                      "Methylation Filter:",
                      choices = c("All Genes", "Biomarkers"),
                      selected = "All Genes"
                    ),
                    br(),
                    radioButtons(
                      ns("slider_type_methylation"),
                      "Slider Type:",
                      choices = c("Methylation Range For At Least 1 Sample", "Top and Bottom Methylation Scores"),
                      selected = "Top and Bottom Methylation Scores"
                    ),
                    br(),
                    uiOutput(ns("methylation_slider"))
                  ),
                  mainPanel(
                    column(
                      12,
                      br(),
                      # selectInput(ns("pageLengthInputMethylation"), "Entries Displayed: ", choices = c(10, 15, 25), selected = 10, width = "165px"),
                      h3("Methylation Data"),
                      HTML(paste0("<b>","Note:","</b>", " The genes values are between 0 and 1 and represent the promoter methylation average beta value based on pre-selected probes (see details in help section).")),
                      br(),br(),
                      DT::dataTableOutput(ns("methylation_data"))
                    )
                  )
                )
              ),class = "hidden-tab"
            ),
            ###----------------------- cnv 
            tabPanel(
              "Copy number",
              fluidRow(
                sidebarLayout(
                  sidebarPanel(
                    width = 2,
                    uiOutput(ns("copy_info")),
                    uiOutput(ns("copy_filter_label")),
                    br(),
                    radioButtons(
                      ns("copy_filter"),
                      "copy Filter:",
                      choices = c("All Genes", "Biomarkers"),
                      selected = "All Genes"
                    ),
                    br(),
                    radioButtons(
                      ns("slider_type_copy"),
                      "Slider Type:",
                      choices = c("Copy Range For At Least 1 Sample", "Top and Bottom copy Scores"),
                      selected = "Top and Bottom copy Scores"
                    ),
                    br(),
                    uiOutput(ns("copy_slider"))
                  ),
                  mainPanel(
                    column(
                      12,
                      br(),
                      # selectInput(ns("pageLengthInputMethylation"), "Entries Displayed: ", choices = c(10, 15, 25), selected = 10, width = "165px"),
                      h3("Copy Number Data"),
                      HTML(paste0("<b>","Note:","</b>", " The genes values are the  average log2 probe intensity ratios based on methylation data (see details in help section).")),
                      br(),br(),
                      DT::dataTableOutput(ns("copy_data"))
                    )
                  )
                )
              ),class = "hidden-tab"
            ),
            ###----------------------- end cnv
            ###----------------------- specific genes
            tabPanel(
              "Specific genes",
              fluidRow(
                sidebarLayout(
                  sidebarPanel(
                    width = 2,
                    uiOutput(ns("spg_info")),
                    uiOutput(ns("spg_filter_label")),
                    br(),
                    textInput(ns("spg_genes"), "Gene Identifiers: (separated by a space e.g. TOP2A MKI67)", "TOP2A MKI67")
                    
                    
                    
                  ),
                  mainPanel(
                    column(
                      12,
                      br(), 
                      # selectInput(ns("pageLengthInputMethylation"), "Entries Displayed: ", choices = c(10, 15, 25), selected = 10, width = "165px"),
                      h3("Looking at specific genes values across all data types"),
                      HTML(paste0("<b>","Note:","</b>", " You can enter one or multiple genes (in the left panel); look at their values in expression, mutation, methylation or copy number; and compare to the median, minimum, maximun values for all patients with same disease")),
                      br(),br(),
                      DT::dataTableOutput(ns("spg_data"))
                    )
                  )
                )
              ),class = "hidden-tab"
            ),
            ###----------------------- end specific genes
            ###------------------------Immune Therapy Treatment 
            tabPanel(
             #  "Potential ADC treatment",
             # "Immune Therapy Treatment",
              "Predictive Biomarkers Details",
              fluidRow(
                sidebarLayout(
                  sidebarPanel(
                    width = 2,
                    uiOutput(ns("trt_info")),
                    uiOutput(ns("trt_filter_label")),
                    br(),
                  ##   textInput(ns("spg_genes"), "Gene Identifiers: (separated by a space e.g. TOP2A MKI67)", "TOP2A MKI67")
                    
                    
                    
                  ),
                  mainPanel(
                    column(
                      12,
                      br(), 
                      # selectInput(ns("pageLengthInputMethylation"), "Entries Displayed: ", choices = c(10, 15, 25), selected = 10, width = "165px"),
                      ## h3("Antibody-Drug Conjugate treatment options"),
                      h3("Immune therapy treatment options"),
                      HTML(paste0("<b>","Note:","</b>", " Drugs were selected based on target gene having expression values (in log2 FPKM) greater or equal than 6")),
                      br(),br(),
                      DT::dataTableOutput(ns("trt_data_v2"))
                    )
                  )
                )
              ),class = "hidden-tab"
            ),
            
            
            ###-------------------------end tretment
            tabPanel(
              "Find Similar Patient or Cell line",
              fluidRow(
                br(), # Add vertical space
                sidebarLayout(
                  sidebarPanel(
                    width = 2,
                    uiOutput(ns("selected_sample_info")),
                    uiOutput(ns("patient_filter_label")),
                    br(),
                    radioButtons(
                      ns("patient_filter"),
                      "Select the Genes to Compute Similarity",
                      choices = c("Top 1000 and Bottom 1000 Genes", "Biomarkers"),
                      selected = "Top 1000 and Bottom 1000 Genes"
                    ),
                    br(),
                    radioButtons(
                      ns("sim_filter"),
                      "Select target samples",
                      choices = c("Patients", "Cell lines"),
                      selected = "Patients"
                    )
                    
                  ),
                  mainPanel(
                    column(12, # Width of the data table column
                           #h3("Patient To Patient Comparison"), 
                           tags$div(
                             style = "font-size: 27px; font-weight: bold;",
                           uiOutput(ns("dynamicPatientHeading"))
                          ),
                          plotOutput(ns("balloonPlot"), height = 600, width = 900), 
                          br(), 
                           withSpinner(DT::dataTableOutput(ns("find_similar_patient_data")))
                    )
                  )
                )
              ),
              class = "hidden-tab"
            )
          )
        )
      )
    )
  )
}

myModule <- function(input, output, session, srcContent) {
  hmethyl = srcContent$patient$molPharmData$met    ## methylation null or not
  hcopy   = srcContent$patient$molPharmData$cop  ## copy number null or not
  currentpatient = ""
  
  output$patIDUi <- renderUI({
    mychoices = c("",unique(srcContent$patient$sampleData$patientID))
    selectizeInput(session$ns("patID"), "",choices = mychoices, multiple = F, selected = myselpat2())
    ## selectizeInput(session$ns("patID"), "",choices = mychoices, multiple = F , selected = "") ## new 
  })
  
  output$patTableUi <- renderUI({
    if (input$showTable) DT::dataTableOutput(session$ns("myPatientSearchData"))
  })
  
  
  biomarkers <- fromJSON("biomarkers.json")

    # results$Annotation <- geneSetPathwayAnalysis::geneAnnotTab[match(results$Gene,rownames(geneSetPathwayAnalysis::geneAnnotTab)), "SHORT_ANNOT"]
    # results$Annotation[is.na(results$Annotation)] <- ""
  
  #GENE EXPRESSION DATA
  #Combine Expression & Annotation Data
  combined_gene_expression <- data.frame(cbind(srcContent$patient$molPharmData$xsqA, srcContent$patient$molPharmData$xsq), check.names = FALSE)
  #Convert From Matrix to Table
  rownames(combined_gene_expression) <- NULL
  
  ## add gene annotation ??? future improvements !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ## combined_gene_expression <- cbind(Annotation = geneSetPathwayAnalysis::geneAnnotTab[match(combined_gene_expression$Symbol,rownames(geneSetPathwayAnalysis::geneAnnotTab)), "SHORT_ANNOT"], combined_gene_expression)  
  
  #MUTATION DATA
  #Combine Mutation & Annotation Data
  combined_mutation <- data.frame(cbind(srcContent$patient$molPharmData$mutA, srcContent$patient$molPharmData$mut), check.names = FALSE)
  
  #Convert From Matrix to Table
  rownames(combined_mutation) <- NULL
  
  #METHYLATION DATA
  #Combine Methylation & Annotation Data
  combined_methylation <- data.frame(cbind(srcContent$patient$molPharmData$metA, srcContent$patient$molPharmData$met), check.names = FALSE)
  #Convert From Matrix to Table
  rownames(combined_methylation) <- NULL
 
  #COPY NUMBER DATA
  #Combine Methylation & Annotation Data
  combined_copy <- data.frame(cbind(srcContent$patient$molPharmData$copA, srcContent$patient$molPharmData$cop), check.names = FALSE)
  #Convert From Matrix to Table
  rownames(combined_copy) <- NULL
   
  #Accessing Patient Search Data
  PatientSearchData <- data.frame(srcContent$patient$sampleData, check.names = FALSE) %>%
    mutate(Number_of_Samples = ave(Name, patientID, FUN = function(x) length(unique(x)))) %>%
    mutate(Treatment_Available = "No") %>%
    select(patientID, gender, age, race, vitalStatus, osMonths, smokingStatus, packYearsSmoking, disease, dataSet, Number_of_Samples, Treatment_Available)
  
  
  #Cleaning Up Patient Search Data
  PatientSearchData$gender <- gsub("Gender:", "", PatientSearchData$gender)  # new gender !!! nov 2023
  PatientSearchData$dataSet <- gsub("DataSet:", "", PatientSearchData$dataSet)
  
  #Ensure No Duplicate Patients & Convert From Matrix to Table
  PatientSearchData <- distinct(PatientSearchData, patientID, .keep_all = TRUE)
  rownames(PatientSearchData) <- NULL
  
  
  
  #Accessing Patient Sample Data
  # PatientSampleData <- data.frame(srcContent$patient$sampleData) %>%
  #   select(Name, patientID, TissueType, OncoTree1, OncoTree2, OncoTree3, OncoTree4, EMT, priorTreatment,
  #          sampleType, dataSource, biopsySite, tumorStage, MSI.TSO500, molecularSubtype, disease) %>%
  #   rename(Sample_Name = Name, Stage = tumorStage) ## 0313
  
  PatientSampleData <- data.frame(srcContent$patient$sampleData) %>%
    select(Name, patientID, TissueType, OncoTree1, OncoTree2, OncoTree3, OncoTree4, EMT, priorTreatment,
           sampleType, sampleDate, dataSource, biopsySite, tumorStage, MSI.TSO500, molecularSubtype, disease) %>%
    rename(Sample_Name = Name, Stage = tumorStage) ## 0313
  
  # Create a new column for genomic information in PatientSampleData
  xsq_information = apply(srcContent$patient$molPharmData$xsq, 2, function(x) { any(!is.na(x))}) 
  xsq_information[which(xsq_information==TRUE)] = "xsq"
  xsq_information[which(xsq_information==FALSE)] = "N/A"
  
  mut_information = apply(srcContent$patient$molPharmData$mut, 2, function(x) { any(!is.na(x))}) 
  mut_information[which(mut_information==TRUE)] = "mut"
  mut_information[which(mut_information==FALSE)] = "N/A"
 
  if (is.null(hmethyl))  
    met_information = replicate(nrow(srcContent$patient$sampleData), "N/A")
  else {
  met_information = apply(srcContent$patient$molPharmData$met, 2, function(x) { any(!is.na(x))}) ## FE
  met_information[which(met_information==TRUE)] = "met"
  met_information[which(met_information==FALSE)] = "N/A"
  }
  
  if (is.null(hcopy))  
    cop_information = replicate(nrow(srcContent$patient$sampleData), "N/A")
  else {
  cop_information = apply(srcContent$patient$molPharmData$cop, 2, function(x) { any(!is.na(x))}) ## FE
  cop_information[which(cop_information==TRUE)] = "cop"
  cop_information[which(cop_information==FALSE)] = "N/A"
  }
  
  genomic_informations = paste(xsq_information, mut_information, met_information, cop_information, sep = ", ")
  
  
  PatientSampleData$"Genomic Information (xsq, mut, met, cop)" = genomic_informations
  
  
  #CleanUp Patient Sample Data (gstring)
  PatientSampleData$EMT <- gsub("EMT:", "", PatientSampleData$EMT)
  PatientSampleData$priorTreatment <- gsub("PriorTreatment:", "", PatientSampleData$priorTreatment)
  PatientSampleData$sampleType <- gsub("SampleType:", "", PatientSampleData$sampleType)
  PatientSampleData$biopsySite <- gsub("BiopsySite:", "", PatientSampleData$biopsySite)
  PatientSampleData$molecularSubtype <- gsub("MolecularSubtype:", "", PatientSampleData$molecularSubtype)
  
  # Extract drug response data from molPharmData$act
  drug_response_data <- data.frame(t(srcContent$patient$molPharmData$act), check.names = FALSE)
  colnames(drug_response_data) <- gsub("^act", "", colnames(drug_response_data))
  
  trtname = apply(drug_response_data, 1, function(x) {i= which(!is.na(x)); if(length(i)>0) return(names(x[i])) else return(NA); })
  drugresponse = apply(drug_response_data, 1, function(x) {i= which(!is.na(x)); if(length(i)>0) return(x[i]) else return(NA); })
  drugresponse[which(drugresponse==1)] = "Good Responder"
  drugresponse[which(drugresponse==0)] = "Bad Responder"
  # PatientSampleData$Treatment_Name <- trtname
  # PatientSampleData$Drug_Response <- drugresponse
  PatientSampleData$Treatment <- trtname
  PatientSampleData$Response_to_treatment <- drugresponse
  
  patient_id <- unique(PatientSampleData$patientID[which(!is.na(trtname))])
  PatientSearchData$Treatment_Available[PatientSearchData$patientID %in% patient_id] <- "Yes"  
  
  
  
  #HORMONE DATA
  #Create Hormone Data
  hormone_data <- srcContent$patient$molPharmData$mda
  ## selected_hormone_data <- hormone_data[c("mdaTMB", "mdaKi67%", "mdaTMB.TSO500", "mdaHormone.prod", "mdaCgA.stain", "mdaSyP.positive", "mdaInsm1"), ]
  selected_hormone_data <- hormone_data["mdaTMB", ,drop=FALSE ] ## just for SCLC website
  
  selected_hormone_data_transpose <- data.frame(t(selected_hormone_data))
  
  #Merge Sample And Hormone Data
  PatientSampleData <- merge(PatientSampleData, selected_hormone_data_transpose,by.x=0, by.y=0) %>%
  rename(TMB = mdaTMB)
  # rename(TMB = mdaTMB,
  #        KI67perc = mdaKi67.,
  #        TMB.TSO500 = mdaTMB.TSO500,
  #        Hormone.prod = mdaHormone.prod,
  #        CgA.stain = mdaCgA.stain,
  #        SyP.positive = mdaSyP.positive,
  #        Insm1 = mdaInsm1)
  
  ## new April 2025 -------------------------------------------------------------
  ## NMF ****
    nmf_data <- srcContent$patient$molPharmData$nmf
  th = 3.8
  nmf_results <- apply(nmf_data,2, function(x) {
    # cat(x," ", is.numeric(x))
    res = ""
    # cat(x[1])
    if (x[1] >= th ) {
      res = "NMF1u,"
    } 
    if (x[2] >= th) {
      res = paste(res,"NMF2u,",sep="")
    }
    if (x[3] >= th) {
      res = paste(res,"NMF3u,",sep="")
    }
    if (x[4] >= th) {
      res = paste(res,"NMF4u,",sep="")
    }
    if (x[5] >= th) {
      res = paste(res,"NMF1r,",sep="")
    }
    if (x[6] >= th) {
      res = paste(res,"NMF2r,",sep="")
    }
    if (x[7] >= th) {
      res = paste(res,"NMF3r,",sep="")
    }
    if (x[8] >= th) {
      res = paste(res,"NMF4r",sep="")
    }
   
    if(res=="") res = "NoCall"
    return(res)
  }
  )
  
  nmf_results = gsub(",$","", nmf_results)
  
 ## PatientSampleData$NMF_clusters <- nmf_results
 
  ## NAPY ***
  
  napy_data <- srcContent$patient$molPharmData$xsq[c("xsqNEUROD1","xsqASCL1","xsqPOU2F3","xsqYAP1"),]
  thre =  5
  napy_results <- apply(napy_data,2, function(x) {
    # cat(x," ", is.numeric(x))
    res = ""
    # cat(x[1])
    if (x[1] >= thre ) {
      res = "NEUROD1,"
    } 
    if (x[2] >= thre) {
      res = paste(res,"ASCL1,",sep="")
    }
    if (x[3] >= thre) {
      res = paste(res,"POU2F3,",sep="")
    }
    if (x[4] >= thre) {
      res = paste(res,"YAP1",sep="")
    }
    
    if(res=="") res = "NoCall"
    return(res)
  }
  )
  
napy_results = gsub(",$","", napy_results)
## PatientSampleData$NAPY_classes <- napy_results

### ADC -------------------------

uniadc = paste0("xsq",unique(adcgenes$gene))
adc_data <- srcContent$patient$molPharmData$xsq[uniadc,]
rownames(adc_data) = substr(rownames(adc_data), 4, nchar(rownames(adc_data)))
thre2 =  6
adc_results <- apply(adc_data,2, function(x) {
  res = ""
   for (k in 1: length(x)) {
     if (x[k] >= thre2) res = paste(res,names(x)[k],",",sep="")
   }
  
  if(res=="") res = "NoCall"
  return(res)
}
)

adc_results = gsub(",$","", adc_results)

groupsnew = data.frame(NMF_clusters  =  nmf_results, NAPY_classes = napy_results,
                       Predictive_biomarkers = adc_results) 
colnames(groupsnew) = c("NMF_clusters (exp>=3.8)", "NAPY_clusters (exp>=5)", "Predictive_biomarkers (exp>=6)")
# PatientSampleData$Predictive_biomarkers <- adc_results
PatientSampleData <- merge(PatientSampleData, groupsnew,by.x=1, by.y=0)

PatientSampleData$MSI.TSO500 <- NULL
PatientSampleData$molecularSubtype <- NULL

## end New --------------------------------------------------------------------

  PatientSampleData$Row.names <- NULL
  
  
  
  ##DISPLAY PATIENT SEARCH DATA, MAKE ROWS CLICKABLE
  
  
  #Display Patient Search Data
  output$myPatientSearchData <- renderDT(
    PatientSearchData,
    filter = "top", style="bootstrap" ,
    rownames = FALSE,
    server = FALSE,
    ## selection = "single",
    selection = list(mode = "single", selected = myselpat()),
    options = list(pageLength = 10)
  )
  
  #output$message <- renderText({
  "Please select a patient" 
  #})
  ## hideTab(inputId = "innerTabs1", target = "Patient Look Up") ## new
  hideTab(inputId = "innerTabs1", target = "Clinical Data")
  hideTab(inputId = "innerTabs1", target = "Gene Expression")
  hideTab(inputId = "innerTabs1", target = "Mutations")
  hideTab(inputId = "innerTabs1", target = "Methylation")
  hideTab(inputId = "innerTabs1", target = "Copy number")
  hideTab(inputId = "innerTabs1", target = "Specific genes")
  hideTab(inputId = "innerTabs1", target = "Predictive Biomarkers Details")
  hideTab(inputId = "innerTabs1", target = "Find Similar Patient or Cell line")
  
  ### new command button ----------------------------------
  
  observeEvent(input$MyPatientTable, {
    
    DT::dataTableOutput("myPatientSearchData")
 
  })
  
  # observeEvent(input$FindMyPatient, {
  #   
  #   ind  = which(PatientSearchData$patientID == trimws(input$patID))
  #   
  #   shiny::validate(need( length(ind) > 0 ,
  #                        "Patient ID not found."))
  #   
  #   ## otherwise go to clinical Data 
  #   
  #   
  # })
  
  ### end command button ------------------------------------------
  
  selected_patient_sample_data <- reactiveVal(NULL)
  
  myselpat  <- reactive({
    ind  = which(PatientSearchData$patientID == trimws(input$patID))
    if (length(ind) > 0  ) return(ind) else return(NULL) ## new &
  })
  
  # myselpat  <- reactive({
  #   ind  = which(PatientSearchData$patientID == trimws(input$patID))
  #   selected_row <- input$myPatientSearchData_rows_selected
  #   if (length(ind) == 0 ) {
  #     return(NULL) 
  #     } 
  #   else {
  #     if (is.null(selected_row)) return(ind)
  #      else {
  #       if (ind!= selected_row) return(ind) else return(NULL) 
  #      }
  #     } 
  # })
  
  myselpat2  <- reactive({
    selected_row <- input$myPatientSearchData_rows_selected
    
    if (!is.null(selected_row)) return(PatientSearchData[selected_row, "patientID"]) else return("")
  })
  
  selrow <- reactive({ return(input$myPatientSearchData_rows_selected) })
  # observeEvent(input$myPatientSearchData_rows_selected | input$FindMyPatient, {
  #   cat("OK checking\n")
  #   selected_row <- input$myPatientSearchData_rows_selected
  #   OK <- F
  #   if (input$FindMyPatient != 0 & is.null(selected_row))  { 
  #     print("processing") 
  #     ind  = which(PatientSearchData$patientID == trimws(input$patID))
  #        if (length(ind) <=0 ) {
  #       shinyalert::shinyalert(
  #         title = paste("PatientID not found"), 
  #         text = "", 
  #         type = "error"
  #       )
  #         } else
  #         { patientID <- trimws(input$patID)
  #          OK = T
  #         }
  #     }
  #   else {
  #     ## selected_row <- input$myPatientSearchData_rows_selected
  #     if (!is.null(selected_row)) {
  #     patientID <- PatientSearchData[selected_row, "patientID"]
  #     OK = T
  #     }
  #   }
  #   
  #   ## selected_row <- input$myPatientSearchData_rows_selected
  #   ## if (!is.null(selected_row)) {
  #   if (OK) {
      
  observeEvent(input$patID, {
    ind  = which(PatientSearchData$patientID == trimws(input$patID))
    if (length(ind) == 0 ) {
      hideTab(inputId = "innerTabs1", target = "Clinical Data")
      hideTab(inputId = "innerTabs1", target = "Gene Expression")
      hideTab(inputId = "innerTabs1", target = "Mutations")
      hideTab(inputId = "innerTabs1", target = "Methylation")
      hideTab(inputId = "innerTabs1", target = "Copy number")
      hideTab(inputId = "innerTabs1", target = "Specific genes")
      hideTab(inputId = "innerTabs1", target = "Predictive Biomarkers Details")
      hideTab(inputId = "innerTabs1", target = "Find Similar Patient or Cell line")
    }
  })
  
  ####
 
  ####
   ####  observeEvent(input$myPatientSearchData_rows_selected,ignoreNULL = FALSE, {
    observeEvent(
      {
       input$myPatientSearchData_rows_selected
       input$patID
       } 
         ,ignoreNULL = FALSE, {
        selected_row <- input$myPatientSearchData_rows_selected
         cat("input list:", input$patID,is.null(input$patID), " input table:", selected_row,is.null(selected_row),"\n")
        ## if (is.null(selected_row) ) {
        if (is.null(selected_row) & is.null(input$patID)) { ## replace & by |
         
          cat("noselection \n")
          hideTab(inputId = "innerTabs1", target = "Clinical Data")
          hideTab(inputId = "innerTabs1", target = "Gene Expression")
          hideTab(inputId = "innerTabs1", target = "Mutations")
          hideTab(inputId = "innerTabs1", target = "Methylation")
          hideTab(inputId = "innerTabs1", target = "Copy number")
          hideTab(inputId = "innerTabs1", target = "Specific genes")
          hideTab(inputId = "innerTabs1", target = "Predictive Biomarkers Details")
          hideTab(inputId = "innerTabs1", target = "Find Similar Patient or Cell line")
        } else
          ## new
        { 
          ### new2
          if (trimws(input$patID)=="")
          {
            cat("noselection2 \n")
            currentpatient <<- "" # new
            hideTab(inputId = "innerTabs1", target = "Clinical Data")
            hideTab(inputId = "innerTabs1", target = "Gene Expression")
            hideTab(inputId = "innerTabs1", target = "Mutations")
            hideTab(inputId = "innerTabs1", target = "Methylation")
            hideTab(inputId = "innerTabs1", target = "Copy number")
            hideTab(inputId = "innerTabs1", target = "Specific genes")
            hideTab(inputId = "innerTabs1", target = "Predictive Biomarkers Details")
            hideTab(inputId = "innerTabs1", target = "Find Similar Patient or Cell line")
          }
          ###
          else
          {
            patientID <- input$patID
          # if (!is.null(selected_row)) patientID <- PatientSearchData[selected_row, "patientID"]
          #  else  patientID <- input$patID
          
          cat(currentpatient," : current pat \n")
          if( currentpatient != patientID)
          ##
          {
            
            currentpatient <<- patientID
            cat(currentpatient," : current pat new\n")
            cat(selected_row," : selected row \n")
      shinyjs::disable("find_similar_patient")
      ## patientID <- PatientSearchData[selected_row, "patientID"]
      showTab(inputId = "innerTabs1", target = "Clinical Data", select = T)
      showTab(inputId = "innerTabs1", target = "Gene Expression")
      showTab(inputId = "innerTabs1", target = "Mutations")
      ## showTab(inputId = "innerTabs1", target = "Methylation")
      ## showTab(inputId = "innerTabs1", target = "Copy number")
      showTab(inputId = "innerTabs1", target = "Specific genes")
      showTab(inputId = "innerTabs1", target = "Predictive Biomarkers Details")
      hideTab(inputId = "innerTabs1", target = "Find Similar Patient or Cell line")
      
      #Clincal Data <- Demographic & Sample Data
      selected_patient_demographic_data <- PatientSearchData[PatientSearchData$patientID == patientID, ]
      selected_patient_sample_data <- PatientSampleData[PatientSampleData$patientID == patientID, ]
      # cat(dim(selected_patient_sample_data), "dim sample data\n")
      # write.csv(selected_patient_sample_data,"test_sample_table.csv")
      selected_patient_sample_data(selected_patient_sample_data)
      
      
      #Update content of the "Clinical Data" tab with the selected patient's data
      output$patient_demographic_info <- renderUI({
        tags$div(
          style = "font-size: 21px; color: #333333;",
          tags$p(
            style = "font-size: 24px; color: blue;",
            #Display PatientID
            paste("Patient ID: ", selected_patient_demographic_data$patientID)
            ##paste("Patient ID: ", selected_patient_demographic_data$patientID, "DataSet: ",selected_patient_demographic_data$dataSet)
          ),
        
        tags$div(
          style = "font-size: 21px; color: #333333; border: 2px solid black; background-color: lightblue;",
          # tags$p(
          #   style = "font-size: 24px; color: blue;",
          #   #Display PatientID
          #   paste("Patient ID: ", selected_patient_demographic_data$patientID)
          #   ##paste("Patient ID: ", selected_patient_demographic_data$patientID, "DataSet: ",selected_patient_demographic_data$dataSet)
          # ),
          tags$div(
            style = "display: flex; justify-content: flex-start; flex-wrap: wrap;",
            tags$span(tags$strong("Gender:"), style = "margin-right: 10px;"),
            selected_patient_demographic_data$gender,
            tags$span(tags$strong("Age:"), style = "margin-right: 3px; margin-left: 25px;"),
            selected_patient_demographic_data$age,
            tags$span(tags$strong("Race:"), style = "margin-right: 3px; margin-left: 25px;"),
            selected_patient_demographic_data$race,
            tags$span(tags$strong("Vital Status:"), style = "margin-right: 3px; margin-left: 25px;"),
            selected_patient_demographic_data$vitalStatus,
            tags$span(tags$strong("DataSet:"), style = "margin-right: 3px; margin-left: 25px;"),
            selected_patient_demographic_data$dataSet
          ),
          tags$br(),
          tags$div(
            style = "display: flex; justify-content: flex-start; flex-wrap: wrap;",
            tags$span(tags$strong("OS Months:"), style = "margin-right: 5px;"),
            selected_patient_demographic_data$osMonths,
            tags$span(tags$strong("Smoking Status:"), style = "margin-right: 3px; margin-left: 25px;"),
            selected_patient_demographic_data$smokingStatus,
            tags$span(tags$strong("Pack Years Smoking:"), style = "margin-right: 3px; margin-left: 25px;"),
            selected_patient_demographic_data$packYearsSmoking,
            tags$span(tags$strong("Disease:"), style = "margin-right: 3px; margin-left: 25px;"),
            selected_patient_demographic_data$disease
          )
        )
        ) # new
      })
      
      #Display Patient Sample Data
      output$patient_sample_data <- renderDT (
        selected_patient_sample_data,
        select = "single",
        extensions = 'Buttons',
        filter='top', style='bootstrap' ,
        options = list(
          dom = 'Bfrtip',
         scrollX = TRUE,
          # dom = 'lipBt' ,
         buttons = list('colvis'),
          # buttons = list(
          #   list(
          #     extend = 'colvis',
          #     text = 'Addtional Information',
          #     postfixButtons = 'colvisRestore'
          # #    # className = 'btn btn-primary'
          #   )
          # ),
          columnDefs = list(
            # list(
            #   targets = c(0, 3, 4, 8, 9, 10, 11, 12, 14, 15),  # Adjusted column indices
            #   visible = TRUE
            # )
            list(
              # targets = c(1, 2, 3,5, 6, 7, 8, 14, 16, 22, 23, 24, 25, 26),  # index first column is ZERO !!!
              ## targets = c(1, 2, 3, 4, 5, 6, 7, 8, 14, 16),
              ## targets = c(1, 2, 3, 4, 5, 6, 7, 8, 12:18),
              targets = c(1:6,8:9,11,13:18),
              visible = FALSE
            )
          )
        ),
        rownames = FALSE
      )
      
      prev_selection <- reactiveVal(NULL)
      observe({
        selected_rows <- input$patient_sample_data_rows_selected
        if (is.null(selected_rows)) {
          if (!is.null(prev_selection())) {
            shinyjs::disable("find_similar_patient")
            hideTab(inputId = "innerTabs1", target = "Find Similar Patient or Cell line")
            
          }
        }
        
        # Update the previous selection with the current selection
        prev_selection(selected_rows)
      })
      
      
      #Switch from Patient Search Data to Clinical Data
      updateTabsetPanel(
        session = session,
        inputId = "innerTabs1",
        selected = "Clinical Data"
      )
          }
        ### end new2 
          }
        ## end new
        }
        
        
        # else
        # { 
        #   cat("no selection \n")
        #   hideTab(inputId = "innerTabs1", target = "Clinical Data")
        #   hideTab(inputId = "innerTabs1", target = "Gene Expression")
        #   hideTab(inputId = "innerTabs1", target = "Mutations")
        #   hideTab(inputId = "innerTabs1", target = "Methylation")
        #   hideTab(inputId = "innerTabs1", target = "Find Similar Patient or Cell line")
        # }
    
  }) ## end double observe event

  #Display Find Similar Patient or Cell line Tab
  observeEvent(input$find_similar_patient, {
    cat(input$patient_sample_data_rows_selected, input$find_similar_patient, "line 1\n")
    shiny::validate(need(!is.null(input$patient_sample_data_rows_selected), "Please select a patient sample"))
    # Check if "Clinical Data" tab is active and a row is selected
    if (input$innerTabs1 == "Clinical Data" && !is.null(input$patient_sample_data_rows_selected)) {
      # Show the "Find Similar Patient or Cell line" tab
      showTab(inputId = "innerTabs1", target = "Find Similar Patient or Cell line")
      updateTabsetPanel(
        session = session,
        inputId = "innerTabs1",
        selected = "Find Similar Patient or Cell line"
      )
    }
    else{
      cat("you need to select a patient sample \n")
    }
  })
  

selected_sample <- reactiveVal(NULL)


observeEvent(input$patient_sample_data_rows_selected, {
  selected_data <- selected_patient_sample_data()
  selected_row <- input$patient_sample_data_rows_selected
  #cat(input$patient_sample_data_rows_selected, "value of input\n")
  if (!is.null(selected_row)) {
    shinyjs::enable("find_similar_patient")
    # Step 2: Update selected_sample with the selected data
    selected_sample_value <- selected_data[selected_row, "Sample_Name"]
    selected_sample(selected_sample_value)
  }
  else {
    #cat("no sample selection...\n")
    shinyjs::disable("find_similar_patient")
  }
  updateRadioButtons(session, "patient_filter", selected = "Top 1000 and Bottom 1000 Genes")
})


output$dynamicPatientHeading <- renderText({
  data_set <- patient_to_patient()
  if(nrow(data_set) > 10) {    ## FE pat 1
  selected_data <- selected_patient_sample_data()
  selected_row <- input$patient_sample_data_rows_selected
  selected_sample_value <- selected_data[selected_row, "Sample_Name"]
  paste0("Table Similarity for Sample ", selected_sample_value, ", Ranked By Sum of Scores (profile, expression, mutation or methylation)")
  }
})



output$selected_sample_info <- renderUI({
  selected_data <- selected_patient_sample_data() # Ensure that a row is selected
  selected_row <- input$patient_sample_data_rows_selected
  selected_sample_value <- selected_data[selected_row, "Sample_Name"]
  tags$div(
    style = "font-size: 16px; color: #333333;",
    tags$p(
      style = "font-size: 20px; color: blue;",
      paste("Sample Name: ",selected_sample_value)
    )
  )
})

output$patient_filter_label <- renderUI({
  req(patientDisease <- PatientSearchData[selected_row(), "disease"])
  req(disease_label <- paste("Disease: ", patientDisease))
  tagList(
    tags$label(disease_label, style = "font-weight: bold;"),
    br()
  )
})

#Load Data For Reactive
mutation_a <- data.frame(srcContent$patient$molPharmData$mut, check.names = FALSE)
current_row_names_mut <- rownames(mutation_a)
new_row_names_mut <- gsub("^mut", "", current_row_names_mut)
rownames(mutation_a) <- new_row_names_mut
mutation_a_ids <- colnames(mutation_a)

expression_a <- cbind(data.frame(srcContent$patient$molPharmData$xsq, check.names = FALSE), data.frame(srcContent$patient$molPharmData$xsqA, check.names = FALSE)) %>%
  filter(biotype == "protein_coding") %>%
  select(-ID, -Symbol, -chr, -start, -end, -strand, -EnsID, -biotype)
current_row_names_xsq <- rownames(expression_a)
new_row_names_xsq <- gsub("^xsq", "", current_row_names_xsq)
rownames(expression_a) <- new_row_names_xsq

xsq_a <- data.frame(srcContent$patient$molPharmData$xsq, check.names = FALSE)
rownames(xsq_a) <- gsub("^xsq", "", rownames(xsq_a))

methylation_a <- data.frame(srcContent$patient$molPharmData$met, check.names = FALSE)
current_row_names_met <- rownames(methylation_a)
new_row_names_met <- gsub("^met", "", current_row_names_met)
rownames(methylation_a) <- new_row_names_met

copy_a <- data.frame(srcContent$patient$molPharmData$cop, check.names = FALSE)
current_row_names_cop <- rownames(copy_a)
new_row_names_cop <- gsub("^cop", "", current_row_names_cop)
rownames(copy_a) <- new_row_names_cop


## CCLE data ------------------------------------------------------------------
mutation_ccle <- data.frame(srcContent_cell$ccle$molPharmData$mut, check.names = FALSE)
rownames(mutation_ccle) <- gsub("^mut", "", rownames(mutation_ccle))
 
expression_ccle <- data.frame(srcContent_cell$ccle$molPharmData$xsq, check.names = FALSE)
rownames(expression_ccle) <- gsub("^xsq", "", rownames(expression_ccle))

methylation_ccle <- data.frame(srcContent_cell$ccle$molPharmData$rrb, check.names = FALSE)
rownames(methylation_ccle) <- gsub("^rrb", "", rownames(methylation_ccle))

# length(intersect(colnames(expression_ccle), ccleinfo$Name)) # 1089
rownames(ccleinfo) = ccleinfo$Name
ccleinfo = ccleinfo[colnames(expression_ccle),]
# stopifnot(identical(colnames(expression_ccle), ccleinfo$Name))

## end CCLE data --------------------------------------------------------------

#Jaccard Similarity Function
jaccard_similarity <- function(A, B) {
  A <- toupper(trimws(A))
  B <- toupper(trimws(B))
  A <- as.vector(na.omit(A))
  B <- as.vector(na.omit(B))
  # A <- toupper(trimws(A))
  # B <- toupper(trimws(B))
  intersection <- length(intersect(A, B))
  ## union <- length(A) + length(B) - intersection
  union <- length(union(A, B))
  return(intersection / union)
}

#Jaccard Profile Function
jaccard_profile <- function(A, B) {
  # we assume A and B have same length
  stopifnot(identical(length(A), length(B)))
  A = toupper(trimws(A))
  B = toupper(trimws(B))
  # if(identical(A,B)){
  #   return (1)
  # }
  #else {
  index  = which(!is.na(A) & !is.na(B))
  ## if (length(index>0)) {
  if (length(index) > 0) {
    compos = length(which(A[index] == B[index]))
    # intersection = length(intersect(A[index], B[index]))
    # return (intersection/length(A))

        ## return (compos/length(A)) old
    return (compos/length(which( !is.na(A) ))) # new
  }
  else
    return (0)
  }
#}

#Load Biomarkers
biomarkers <- fromJSON("biomarkers.json")

#Load Sample Data
sample_data <- data.frame(srcContent$patient$sampleData, check.names = FALSE) %>%
  select("Name", "disease", "sampleType", "gender", "age", "race", "biopsySite", "OncoTree1", "OncoTree2", "tumorStage") %>%
  mutate(age_group = case_when(
    age >= 0 & age <= 25 ~ "Young",
    age > 25 & age <= 65 ~ "Middle-aged",
    age > 65 ~ "Elderly",
    TRUE ~ "Unknown"
  )) %>%
  as.data.frame() %>%
  select(-age)

sample_data$sampleType[which(sample_data$sampleType=="SampleType:not reported")]= NA


sample_data_v2 <- sample_data

sample_data_v2$biopsySite = gsub("BiopsySite:","", sample_data_v2$biopsySite)
sample_data_v2$gender = gsub("Gender:","", sample_data_v2$gender)  ## new nov 2023

sample_data_v2$biopsySite[which(sample_data_v2$biopsySite=="plural effusion")] = "pleural effusion"

sample_data_v2$sampleType[which(sample_data_v2$sampleType=="SampleType:additional new primary tumor")]="Primary"
sample_data_v2$sampleType[which(sample_data_v2$sampleType=="SampleType:Primary")]="Primary"
sample_data_v2$sampleType[which(sample_data_v2$sampleType=="SampleType:Relapse primary site")]="Primary"
sample_data_v2$sampleType[which(sample_data_v2$sampleType=="SampleType:Primary metastasis")]="Metastatic"
sample_data_v2$sampleType[which(sample_data_v2$sampleType=="SampleType:Relapse metastasis")]="Metastatic"
sample_data_v2$sampleType[which(sample_data_v2$sampleType=="SampleType:normal")]="Normal"
sample_data_v2$sampleType[which(sample_data_v2$sampleType=="SampleType:not reported")]="NotReported"

## write.csv(sample_data_v2,"sample_data_v2.csv")

cell_data <- ccleinfo %>%
    mutate(age_group = case_when(
    Age >= 0 & Age <= 25 ~ "Young",
    Age > 25 & Age <= 65 ~ "Middle-aged",
    Age > 65 ~ "Elderly",
    TRUE ~ "Unknown" )) %>%
    as.data.frame() %>% select(-Age)

cell_data <- cell_data[,c("Name","PrimaryOrMetastasis","Sex","SampleCollectionSite","OncoTree1","OncoTree2","age_group")]
colnames(cell_data) <- c("Name","sampleType","gender","biopsySite","OncoTree1","OncoTree2","age_group")

## write.csv(cell_data,"cell_data.csv")    

patient_to_patient <- reactive ({ 
    shinyjs::toggle("balloonPlot", condition = FALSE)
    shinyjs::toggle("dynamicPatientHeading", condition = FALSE)
    shinyjs::toggle("patient_filter_label", condition = FALSE)
    shinyjs::toggle("selected_sample_info", condition = FALSE)
    
    # wexp = 0.8
    # wmut = 0.1
    # wmeth = 0
    # wprofile = 0.1
    # 
    # shiny::validate(need(wexp+wmut+wmeth+wprofile==1,"Total weight should be equal 1")) ## FE new
    
    start_time <- Sys.time() 
    selected_sample_value <- selected_sample()
    
    if (is.null(selected_sample_value)) {
      # Handle the case when no sample is selected
      return(NULL)
    }
    else{
    non_zero_mutation <- mutation_a %>%
      filter(!!sym(selected_sample_value) != 0)
    # test
    cat(selected_sample_value,"\t",dim(mutation_a), "\t", dim(non_zero_mutation), "\n")
    #rownames(sample_data) <- NULL
    patientDisease <- PatientSearchData[selected_row(), "disease"] 
    
    # Perform Correlation for Expression
    patient_filter <- input$patient_filter
    #print(patient_filter)
    if("Biomarkers" %in% patient_filter) {
      if (patientDisease %in% names(biomarkers)) {
        expression_a_filtered <- expression_a[rownames(expression_a) %in% biomarkers[[patientDisease]][["gene_expression"]], ]
        if(nrow(expression_a_filtered) < 3 && nrow(expression_a) != 0) {
          shinyalert::shinyalert(
            title = paste("Not enough biomarkers present"), 
            text = paste("We have", nrow(expression_a_filtered), " gene expression biomarkers present in this sample, thus the gene expression column will not be included"), 
            type = "info"
          )
        }
        #print(nrow(mutation_a))
        mutation_a_filtered <- non_zero_mutation[rownames(non_zero_mutation) %in% biomarkers[[patientDisease]][["mutation"]], ]
        #print(nrow(mutation_a_filtered))
        #print(head(mutation_a_filtered))
       ## old  if(nrow(mutation_a_filtered) < 3 && nrow(non_zero_mutation) != 0) {
          if(nrow(mutation_a_filtered) < 1) {
          shinyalert::shinyalert(
          title = paste("Not enough biomarkers present"), 
          text = paste("We have", nrow(mutation_a_filtered), " mutation biomarkers present in this sample, thus the mutation column will not be included"), 
          type = "info"
          )
        }
        #print(nrow(mutation_a_filtered))
        if (!is.null(hmethyl)) 
        {
        methylation_a_filtered <- methylation_a[rownames(methylation_a) %in% biomarkers[[patientDisease]][["methylation"]], ]
        if(nrow(methylation_a_filtered) < 3 && nrow(methylation_a) != 0) {
          shinyalert::shinyalert(
            title = paste("Not enough biomarkers present"), 
            text = paste("We have", nrow(methylation_a_filtered), " methylation biomarkers present in this sample, thus the methylation column will not be included"), 
            type = "info"
          )
        }
        }
        ## copy
        # copy_a_filtered <- copy_a[rownames(copy_a) %in% biomarkers[[patientDisease]][["copy_number"]], ]
        # if(nrow(copy_a_filtered) < 3 && nrow(copy_a) != 0) {
        #   shinyalert::shinyalert(
        #     title = paste("Not enough biomarkers present"), 
        #     text = paste("We have", nrow(copy_a_filtered), " copy number biomarkers present in this sample, thus the copy number column will not be included"), 
        #     type = "info"
        #   )
        # }
        ## end copy
      }
      else{
        shinyalert::shinyalert(
          title = paste("Biomarkers for ", patientDisease, " are not available"), 
          text = paste("We do not have biomarkers for the disease ", patientDisease, ". Please return to the Top 1000 and Bottom 1000 Genes Filter."), 
          type = "info"
        )
        shiny::validate(
          need(FALSE, paste("We do not have biomarkers for the disease ", patientDisease, ". Please return to the Top 1000 and Bottom 1000 Genes Filter."))
        )
      } 
      
    }   ## end biomarkers
        
    else {
    ## top 2000. genes except mutation
    ## expression selection ----------------------------------------
    expression_a_filtered <- expression_a
    mutation_a_filtered <- non_zero_mutation
    ## methylation_a_filtered <- methylation_a
    ## copy_a_filtered <- copy_a
    
    top_genes_expr <- expression_a_filtered %>%
      arrange(desc(!!sym(selected_sample_value))) %>%
      head(1000)
    bottom_genes_expr <- expression_a_filtered %>%      ## FE pat
      arrange(!!sym(selected_sample_value)) %>%
      slice(1:1000)
  
    expression_a_filtered <- rbind(top_genes_expr, bottom_genes_expr)
    
    ## methylation selection .......................
    if (!is.null(hmethyl)) 
    {
    methylation_a_filtered <- methylation_a
    top_genes_met <- methylation_a_filtered %>%
      #select({{ selected_sample_value }}) %>%
      arrange(desc(!!sym(selected_sample_value))) %>%
      head(1000)
    #print(head(top_genes_expr))
    bottom_genes_met <- methylation_a_filtered %>%
      #select({{ selected_sample_value }}) %>%
      arrange(!!sym(selected_sample_value)) %>%
      slice(1:1000)
    
    methylation_a_filtered <- rbind(top_genes_met, bottom_genes_met)
    }
    ## copy selection .......................Not used
    
    ## no mutation further selection

    }  ## end filtering based on variable genes or bio markers ..................................
    ## end selection ........................................................................
    
    ## compute now correlation with the selected sample .........................................
    ## MIN 3 observations .. could work with 2 but not p-value
    
    ## sample_vector <- expression_a_filtered[, selected_sample_value]  ## vector
    
 ## Final gene selection for our patient sample -------------------------------------------------------------   
    sample_vector <- expression_a_filtered %>%
      select({{ selected_sample_value }})  ## matrix + we have gene names
    
    sample_vector_mutation <- mutation_a_filtered %>%
      select({{ selected_sample_value }})
    
    if (!is.null(hmethyl)) 
    {
    sample_vector_methylation <- methylation_a_filtered %>%
      select({{ selected_sample_value }})
    }
 ## --------------------------------------------------------------------------------------------------------    
 
 ## Patient or cell line similarity
    
  if (input$sim_filter=="Patients") ## patient similarity
    {
    
    #if (nrow(expression_a_filtered) > 0 && nrow(expression_a_filtered) < 3) {
    if (nrow(expression_a_filtered) < 3) { 
      expression_cor <- NA
    }
    else{
    expression_cor <- cor(expression_a_filtered, sample_vector, use = "pairwise.complete.obs")
    rownames(expression_cor) <- NULL
    }
  
    
    # sample_vector_mutation <- mutation_a_filtered %>%
    #   select({{ selected_sample_value }})
    cat(nrow(mutation_a_filtered) ,"\t", nrow(sample_vector_mutation),"****mutation****\n" )
    write.csv(mutation_a_filtered,"mut_a_filtered.csv")
    
    ## new mutation calculation 
    if (nrow(mutation_a_filtered) < 1) {
      mutation_cor <- NA
    }
    else{
      mutation_cor <- apply(mutation_a_filtered, 2 ,function(x) {length(which(x>0)) / nrow(sample_vector_mutation) } )
      # mutation_cor[which(mutation_cor==0)] = NA
      rownames(mutation_cor) <- NULL
    }
   
     ## old computation
    # if (nrow(mutation_a_filtered) < 3) {
    #   mutation_cor <- NA
    # }
    # else{
    #   mutation_cor <- cor(mutation_a_filtered, sample_vector_mutation, use = "pairwise.complete.obs")
    #   rownames(mutation_cor) <- NULL
    # }
    ## end old
    
    # Check if there are still enough rows to calculate the correlation
    # cat(nrow(mutation_a_filtered) ,"\t", nrow(sample_vector_mutation),"****\n" )
    # if (nrow(mutation_a_filtered) > 1 && nrow(sample_vector_mutation) > 1) {
    #   mutation_cor <- cor(mutation_a_filtered, sample_vector_mutation, use = "pairwise.complete.obs")
    #   rownames(mutation_cor) <- NULL
    #   if(nrow(mutation_a_filtered) > 0 && nrow(mutation_a_filtered) < 3) {
    #     mutation_cor <- NA
    #   }
    # } else {
    #   # Handle the case when there are not enough data points for correlation
    #   mutation_cor <- NA
    # }
    
    # Perform Correlation for Expression
    # top_genes_met <- methylation_a_filtered %>%    ## not this place ## FE
    #   #select({{ selected_sample_value }}) %>%
    #   arrange(desc(!!sym(selected_sample_value))) %>%
    #   head(1000)
    # #print(head(top_genes_expr))
    # 
    # 
    # bottom_genes_met <- methylation_a_filtered %>%
    #   #select({{ selected_sample_value }}) %>%
    #   arrange(!!sym(selected_sample_value)) %>%
    #   slice(1:1000)
    # 
    # #print(head(bottom_genes_expr))
    # 
    # 
    # methylation_a_filtered <- rbind(top_genes_met, bottom_genes_met)
    
    # sample_vector_methylation <- methylation_a_filtered %>%
    #   select({{ selected_sample_value }})
    
    ## if (nrow(methylation_a_filtered) > 0 && nrow(methylation_a_filtered) < 3) {
    methylation_cor <- NA
    if (!is.null(hmethyl))
    {
    if (nrow(methylation_a_filtered) < 3) {
      methylation_cor <- NA
    }
    else{
      methylation_cor <- cor(methylation_a_filtered, sample_vector_methylation, use = "pairwise.complete.obs")
      rownames(methylation_cor) <- NULL
    }
    }
   
    data <- data.frame(Sample_Name = mutation_a_ids, Expression_Score = expression_cor, Mutation_Score = mutation_cor, Methylation_Score = methylation_cor)
    colnames(data) <- c("Sample_Name", "Expression_Score", "Mutation_Score", "Methylation_Score")
    
    #Compute Jaccard and profile score ................................................................
    
    # Compute Jaccard similarity with other samples
    selected_sample_info <- sample_data[sample_data$Name == selected_sample_value, ] %>%
      select(-Name)
    
    # selected_sample_info <- as.vector(selected_sample_info)
    selected_sample_info <- as.character(selected_sample_info)
    
    jaccard_scores_similarity <- numeric(length(sample_data$Name))
    #jaccard_scores_similarity[1] = 1
    jaccard_scores_profile <- numeric(length(sample_data$Name))
    #jaccard_scores_profile[1] = 1
  
    ## FE
    for (i in 1:length(sample_data$Name)) {
      sample <- sample_data$Name[i]
      sample_data_subset <- sample_data[i, -1] ## remove the name
      
     if(sample == selected_sample_value) {
        jaccard_scores_similarity[i] <- 1
        jaccard_scores_profile[i] <- 1
      }
      else{
        jaccard_scores_similarity[i] <- jaccard_similarity(selected_sample_info, as.character(sample_data_subset))
        jaccard_scores_profile[i]    <- jaccard_profile(selected_sample_info, as.character(sample_data_subset))
      }
    }
    
    
    # for (i in 1:length(sample_data$Name)) {
    #   sample <- sample_data$Name[i]
    #   sample_data_subset <- sample_data[sample_data$Name == sample, ]
    #   sample_data_subset <- sample_data_subset %>%
    #     select(-Name)
    #   # sample_data_subset <- sample_data[i, 2:10]
    #   
    #   # sample_data_subset <- sample_data_subset[, !names(sample_data_subset) %in% "Name"]
    #   if(sample == selected_sample_value) {
    #     jaccard_scores_similarity[i] <- 1
    #   }
    #   else{
    #   jaccard_scores_similarity[i] <- jaccard_similarity(selected_sample_info, as.vector(sample_data_subset))
    #   }
    # }
    
    # for (i in 1:length(sample_data$Name)) {
    #   sample <- sample_data$Name[i]
    #   sample_data_subset <- sample_data[sample_data$Name == sample, ]
    #   sample_data_subset <- sample_data_subset %>%
    #     select(-Name)
    #   if(sample == selected_sample_value) {
    #     jaccard_scores_profile[i] <- 1
    #   }
    #   else{
    #   jaccard_scores_profile[i] <- jaccard_profile(selected_sample_info, as.vector(sample_data_subset))
    #   }
    # }
    
    
    
    # Assign the Jaccard similarity scores to the data frame
    data$Jaccard_Similarity <- jaccard_scores_similarity
    data$Common_Feature_Score <- jaccard_scores_profile
    data <- subset(data, !duplicated(Sample_Name))
    
    data$Mean_Of_Genomic_Features  <- rowMeans(data[, c("Expression_Score", "Mutation_Score", "Methylation_Score")], na.rm = TRUE)
    
    data$SumScore  <- rowSums(data[, c("Expression_Score", "Mutation_Score", "Methylation_Score","Common_Feature_Score")], na.rm = TRUE)
    
    # old
    # data <- data %>%
    #   select("Sample_Name", "Jaccard_Similarity", "Common_Feature_Score", "Expression_Score", "Mutation_Score", "Methylation_Score" ,"Mean_Of_Genomic_Features") %>%
    #   arrange(desc(Mean_Of_Genomic_Features))
    
    # new
    data <- data %>%
      select("Sample_Name", "Common_Feature_Score", "Expression_Score", "Mutation_Score", "Methylation_Score" ,"SumScore","Mean_Of_Genomic_Features") %>%
      arrange(desc(SumScore))
    
    
    ## data <- left_join(data, sample_data[, c("Name", "OncoTree1", "OncoTree2")], by = c("Sample_Name" = "Name")) ## FE new

    ##    data <- left_join(data, PatientSampleData[, c("Sample_Name", "OncoTree1", "OncoTree2","Treatment_Name","Drug_Response")], by = c("Sample_Name" = "Sample_Name"))
        data <- left_join(data, PatientSampleData[, c("Sample_Name", "OncoTree1", "OncoTree2","Treatment","Response_to_treatment")], by = c("Sample_Name" = "Sample_Name"))
        
    
    # data <- data[, c("Sample_Name", "OncoTree1", "OncoTree2", setdiff(names(data), c("Sample_Name", "OncoTree1", "OncoTree2")))] %>%
    #   rename(Trt_Response = Drug_Response,
    #          Trt_Name = Treatment_Name)
    
    # todisplay = c("Sample_Name", "OncoTree1", "OncoTree2","Mean_Of_Genomic_Features","Jaccard_Similarity","Treatment_Name","Drug_Response")
    # data <- data[, c(todisplay, setdiff(names(data), todisplay))] %>%
    #   rename(Trt_Response = Drug_Response,
    #          Trt_Name = Treatment_Name)
    
    ## old
    # todisplay = c("Sample_Name", "OncoTree1", "OncoTree2","Mean_Of_Genomic_Features","Jaccard_Similarity","Treatment","Response_to_treatment")
    ## new
    todisplay = c("Sample_Name", "OncoTree1", "OncoTree2","SumScore","Treatment","Response_to_treatment")
    
    data <- data[, c(todisplay, setdiff(names(data), todisplay))]
      
    
    # numeric_columns <- sapply(data, is.numeric)
    # data[numeric_columns] <- round(data[numeric_columns], 3)
  }
    else ## cell line similariry +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    { 
      xsqgenes = intersect(rownames(expression_ccle), rownames(sample_vector))
      mutgenes = intersect(rownames(mutation_ccle), rownames(sample_vector_mutation))
      cat("xsq genes: ", length(xsqgenes), " mut genes: ", length(mutgenes), "\n")
       if (!is.null(hmethyl)) 
        metgenes = intersect(rownames(methylation_ccle), rownames(sample_vector_methylation))
       else  metgenes = c()
      
      if (length(xsqgenes) < 3) { 
        expression_cor <- NA
      }
      else{
        expression_cor <- cor(expression_ccle[xsqgenes,], sample_vector[xsqgenes,], use = "pairwise.complete.obs")
        rownames(expression_cor) <- NULL
      }
      ## old computation
      # if (length(mutgenes) < 3) {
      #   mutation_cor <- NA
      # }
      # else{
      #   mutation_cor <- cor(mutation_ccle[mutgenes,], sample_vector_mutation[mutgenes,], use = "pairwise.complete.obs")
      #   rownames(mutation_cor) <- NULL
      # }
      ## end old
      ## new mutation computation
      if (length(mutgenes) < 1) {
        mutation_cor <- NA
      }
      else{
        ## mutation_cor <- cor(mutation_ccle[mutgenes,], sample_vector_mutation[mutgenes,], use = "pairwise.complete.obs")
        mutation_cor <- apply(mutation_ccle[mutgenes,], 2 ,function(x) {length(which(x>0)) / length(mutgenes) } )
        rownames(mutation_cor) <- NULL
      }
      
      ##
      if (length(metgenes) < 3) {
        methylation_cor <- NA
      }
      else{
        methylation_cor <- cor(methylation_ccle[metgenes,], sample_vector_methylation[metgenes,], use = "pairwise.complete.obs")
        rownames(methylation_cor) <- NULL
      }
      
      data <- data.frame(Sample_Name = colnames(expression_ccle), Expression_Score = expression_cor, Mutation_Score = mutation_cor, Methylation_Score = methylation_cor)
      colnames(data) <- c("Sample_Name", "Expression_Score", "Mutation_Score", "Methylation_Score")
      
      #Compute Jaccard and profile score ................................................................
      
      # Compute Jaccard similarity with other samples
      selected_sample_info <- sample_data_v2[sample_data_v2$Name == selected_sample_value, ] %>%
        select(-Name, -disease, -race, -tumorStage)
      
      selected_sample_info <- as.character(selected_sample_info)
      
      jaccard_scores_similarity <- numeric(length(cell_data$Name))
      #jaccard_scores_similarity[1] = 1
      jaccard_scores_profile <- numeric(length(cell_data$Name))
      #jaccard_scores_profile[1] = 1
      
      ## FE
      for (i in 1:length(cell_data$Name)) {
        if (cell_data[i,1] == "COR-L95" )
        {
          c_info = as.character(cell_data[i, -1])
          p_info = selected_sample_info
        }
        sample_data_subset <- cell_data[i, -1] ## remove the name
       
        jaccard_scores_similarity[i] <- jaccard_similarity(selected_sample_info, as.character(sample_data_subset))
        jaccard_scores_profile[i]    <- jaccard_profile(selected_sample_info, as.character(sample_data_subset))
       
      }
      
      ##
      # Assign the Jaccard similarity scores to the data frame
      data$Jaccard_Similarity <- jaccard_scores_similarity
      data$Common_Feature_Score <- jaccard_scores_profile
      data <- subset(data, !duplicated(Sample_Name))
      
      data$Mean_Of_Genomic_Features  <- rowMeans(data[, c("Expression_Score", "Mutation_Score", "Methylation_Score")], na.rm = TRUE)
      
      # data <- data %>%
      #   select("Sample_Name", "Jaccard_Similarity", "Common_Feature_Score", "Expression_Score", "Mutation_Score", "Methylation_Score" ,"Mean_Of_Genomic_Features") %>%
      #   arrange(desc(Mean_Of_Genomic_Features))
      # 
      ## new
      
      data$SumScore  <- rowSums(data[, c("Expression_Score", "Mutation_Score", "Methylation_Score","Common_Feature_Score")], na.rm = TRUE)
   
      data <- data %>%
        select("Sample_Name", "Common_Feature_Score", "Expression_Score", "Mutation_Score", "Methylation_Score" ,"SumScore","Mean_Of_Genomic_Features") %>%
        arrange(desc(SumScore))
      
      ## end new
      
      
      ## data <- left_join(data, cell_data[, c("Name", "OncoTree1", "OncoTree2")], by = c("Sample_Name" = "Name")) ## FE new
      data <- left_join(data, ccleinfo[, c("Name", "OncoTree1", "OncoTree2","Top_Sensitive_Drugs","Bottom_Resistant_Drugs","nb_drugs")], by = c("Sample_Name" = "Name")) 
      ## data <- left_join(data, PatientSampleData[, c("Sample_Name", "OncoTree1", "OncoTree2","Treatment_Name","Drug_Response")], by = c("Sample_Name" = "Sample_Name"))
      ## old
      ## todisplay = c("Sample_Name", "OncoTree1", "OncoTree2", "Mean_Of_Genomic_Features","Jaccard_Similarity","Top_Sensitive_Drugs")
      
      todisplay = c("Sample_Name", "OncoTree1", "OncoTree2", "SumScore","Top_Sensitive_Drugs")
      
      ## data <- data[, c("Sample_Name", "OncoTree1", "OncoTree2", setdiff(names(data), c("Sample_Name", "OncoTree1", "OncoTree2")))] 
      data <- data[, c(todisplay, setdiff(names(data), todisplay))] 
      
      
      
      # numeric_columns <- sapply(data, is.numeric)
      # data[numeric_columns] <- round(data[numeric_columns], 3)
      ##
            
      } ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ end cell line 
    
    numeric_columns <- sapply(data, is.numeric)
    data[numeric_columns] <- round(data[numeric_columns], 3)
    
    end_time <- Sys.time()
    
    execution_time <- end_time - start_time
    cat("Execution Time:", execution_time, "\n")
    
    data <- data %>%
      select_if(~ !all(is.na(.))) %>%
      distinct() 
    
   
    shinyjs::toggle("balloonPlot", condition = TRUE)
    shinyjs::toggle("dynamicPatientHeading", condition = TRUE)
    shinyjs::toggle("patient_filter_label", condition = TRUE)
    shinyjs::toggle("selected_sample_info", condition = TRUE)
    
    
    write.csv(data,"test_balloon.csv")
    return(data)
    }
  })
  

output$balloonPlot <- renderPlot({
  res <- patient_to_patient()
  available_columns <- colnames(res)
  ## columns_to_select <- c("Sample_Name", "Jaccard_Similarity", "Common_Feature_Score")
  columns_to_select <- c("Sample_Name", "Common_Feature_Score")
  
  if (exists("Expression_Score", where = res)) {
    columns_to_select <- c(columns_to_select, "Expression_Score")
  }
  
  if (exists("Mutation_Score", where = res)) {
    columns_to_select <- c(columns_to_select, "Mutation_Score")
  }
  
  if (exists("Methylation_Score", where = res)) {
    columns_to_select <- c(columns_to_select, "Methylation_Score")
  }
  
  columns_to_select <- c(columns_to_select, "Mean_Of_Genomic_Features")
  columns_to_select <- c(columns_to_select, "SumScore")
  
  res2 <- res %>% select(all_of(columns_to_select))
  colnames(res2)[1] <- "Cell_Line_Name"
  
  ## top 20
  res2 <- res2[1:20, ]
  
  melt_res2 <- melt(res2, id.vars = "Cell_Line_Name")  
  #print(sum(available_columns %in% c("Expression_Score", "Mutation_Score", "Methylation_Score")))
  #old
  # melt_res2 <- transform(melt_res2, Cell_Line_Name = reorder(melt_res2$Cell_Line_Name, melt_res2$value, FUN = function(x) {mean(x[3:(2+sum(available_columns %in% c("Expression_Score", "Mutation_Score", "Methylation_Score")))], na.rm = TRUE)}))
  # new
  melt_res2 <- transform(melt_res2, Cell_Line_Name = reorder(melt_res2$Cell_Line_Name, melt_res2$value, FUN = function(x) {sum(x[1:(0+sum(available_columns %in% c("Expression_Score", "Mutation_Score", "Methylation_Score","Common_Feature_Score")))], na.rm = TRUE)}))
  
  p <- plot_balloon_plot(melt_res2, "Tumor To Tumor Similarity")
  p <- p + ylab("Top 20 Tumor Samples") + xlab("Unweighted Similarity Scores By Data Type") + theme(axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 17), axis.text = element_text(face="bold"), axis.title = element_text(face="bold"), title = element_text(face="bold")) 
  #p <- p + geom_label(aes_string(label = 'value', x = 'numeric_variable'), alpha = 1.0, size = 6, fill = "white",  nudge_x = 0.1)

 
  p
  
})

  
  output$find_similar_patient_data <- DT::renderDataTable({
    # datatable(
    #   patient_to_patient(), rownames = F,
    #   filter='top', style='bootstrap', selection = "none",extensions = "FixedColumns",
    #   options=list(scrollX = TRUE, fixedColumns = list(leftColumns = 1), lengthMenu = c(10, 25, 50, 100), pageLength = 10, language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')))
    # )
    
    datatable(
      patient_to_patient(), rownames = F,
      filter='top', style='bootstrap', selection = "none",
      options=list(scrollX = TRUE,lengthMenu = c(10, 25, 50, 100), pageLength = 10, language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')))
    )
    
    # datatable(
    #   filter = "top", 
    #   patient_to_patient(),
    #   style='bootstrap', selection = "none",
    #   options = list(
    #     columnDefs = list(
    #       list(visible = FALSE, targets = c(0, 5, 6, 7, 8))  
    #     ),
    #     dom = 'Bfrtip',  # Displays the column visibility button
    #     buttons = list('colvis'),  # Adds the column visibility button
    #     pageLength = 10, 
    #     lengthMenu = c(10, 25, 50, 100)
    #   ))
    
  })


  
  #Display Gene Expression Data
  
  # Get selected patientID
  selected_row <- reactive(input$myPatientSearchData_rows_selected)
  selected_pat <- reactive({
    res = input$patID
    if (is.null(res))
        {
          sel_row  = input$myPatientSearchData_rows_selected
          if (!is.null(sel_row)) {
             res = PatientSearchData[sel_row, "patientID"]
          }
         }
     return(res)
    }) ## new
  
  # Display patientID
  output$gene_expression_info <- renderUI({
    # req(selected_row())  # Ensure that a row is selected
    req(selected_pat())  # Ensure that a row is selected
    tags$div(
      style = "font-size: 16px; color: #333333;",
      tags$p(
        style = "font-size: 20px; color: blue;",
      #  paste("Patient ID: ",  PatientSearchData[selected_row(), "patientID"])
        paste("Patient ID: ",  PatientSearchData[which(PatientSearchData$patientID==selected_pat()), "patientID"])
      )
    )
  })
  
  
  # Reset input values when a new patient is selected
  observeEvent(selected_row(), {
    updateRadioButtons(session, "gene_filter", selected = "Protein Coding")  # Reset the gene filter radio buttons to default
    updateRadioButtons(session, "slider_type", selected = "Gene Expression Range")  # Reset the slider type radio buttons to default
    #updateSliderInput(session, "gene_expr_range", value = c(-1, 15))  # Reset the gene expression range slider to default values
    updateSliderInput(session, "num_genes_displayed", value = 10)  # Reset the number of genes displayed slider to default value
  })

  
  output$gene_filter_label <- renderUI({
    # patientDisease <- PatientSearchData[selected_row(), "disease"] 
    patientDisease <- PatientSearchData[which(PatientSearchData$patientID == selected_pat() ), "disease"] 
    disease_label <- paste("Disease: ", patientDisease)
    tagList(
      tags$label(disease_label, style = "font-weight: bold;"),
      br()
    )
  })
  
  my_gene_expression <- reactive({
    # if (!is.null(selected_row())) {  # Check if a row is selected
    #   patientID <- PatientSearchData[selected_row(), "patientID"]  # Get the patient ID of the selected row
    if (!is.null(selected_pat())) {  # Check if a row is selected
      patientID <- selected_pat()  # Get the patient ID of the selected row
      selected_samples <- PatientSampleData$Sample_Name[PatientSampleData$patientID == patientID]  # Get the selected patient's sample names
      gene_expression_cols <- c("ID", "Symbol", "chr", "start", "end", "strand", "EnsID", "biotype", selected_samples)  # Define the columns to select from gene expression data
      selected_gene_expression <- combined_gene_expression[, gene_expression_cols, drop = FALSE]  # Select the specified columns from the combined gene expression data
      
      sample_columns_expression <- selected_gene_expression[, selected_samples, drop = FALSE] %>%
        select_if(~ !all(is.na(.)))
      shiny::validate(need(ncol(sample_columns_expression)>0,"No expression data for this patient")) ## FE new
      
      selected_samples <- colnames(sample_columns_expression)
      
      gene_expression_cols <- c("ID", "Symbol", "chr", "start", "end", "strand", "EnsID", "biotype", selected_samples)
      selected_gene_expression <- selected_gene_expression[, gene_expression_cols, drop = FALSE]
      
      
      ## patientDisease <- PatientSearchData[selected_row(), "disease"] 
      patientDisease <- PatientSearchData[which(PatientSearchData$patientID==selected_pat()), "disease"] 
      
      # Apply gene filter
      gene_filter <- input$gene_filter  # Get the selected gene filter option
      
      # Read the biomarkers JSON file
      biomarkers <- fromJSON("biomarkers.json")
      
      # Access the gene_expression lists for ACC and SCLC
      if ("Protein Coding" %in% gene_filter) {  # If "Protein Coding" option is selected
        filtered_gene_expression <- selected_gene_expression[selected_gene_expression$biotype == "protein_coding", ]  # Filter the gene expression data to include only protein coding genes
      } else if ("Biomarkers" %in% gene_filter) {  # If "Biomarkers" option is selected
        if (patientDisease %in% names(biomarkers)) {
          filtered_gene_expression <- selected_gene_expression[selected_gene_expression$ID %in% biomarkers[[patientDisease]][["gene_expression"]], ]
        } else {
          #Display shinyalert with "No biomarkers found" message
          shinyalert::shinyalert(
            title = "No Biomarkers Found",
            text = paste("We currently do not have gene expression biomarkers data for the disease:", patientDisease, ". Please return to either the 'All Genes' or 'Protein Coding' filter"),
            type = "info"
          ) # Use the original selected gene expression data without any filter
          shiny::validate(
            need(FALSE, "No biomarkers found. Please return to either the protein coding or all genes filter.")
          )
        }
      } else {
        filtered_gene_expression <- selected_gene_expression  # Use the original selected gene expression data without any filter
      }
      
      # Compute mean and standard deviation columns
      num_samples <- length(selected_samples)  # Get the number of selected samples
      
      if (num_samples >= 2) {  # If there are at least 2 samples
        filtered_gene_expression$mean <- round(rowMeans(filtered_gene_expression[, selected_samples], na.rm = TRUE), 3)  # Compute the mean expression across samples for each gene
      }
      
      if (num_samples >= 3) {  # If there are at least 3 samples
        filtered_gene_expression$std_dev <- round(apply(filtered_gene_expression[, selected_samples], 1, sd, na.rm = TRUE), 3)  # Compute the standard deviation across samples for each gene
      }

        selected_gene_expression <- filtered_gene_expression
        
        if (num_samples == 1) {  # If there is only 1 selected sample
          # Handle the case when num_samples == 1
          # Only default columns and selected sample
          columns_to_display <- c("ID", "Symbol", "chr", "start", "end", "strand", "EnsID", "biotype", selected_samples)
          selected_gene_expression <- selected_gene_expression[, columns_to_display, drop = FALSE]  # Select the specified columns
          selected_gene_expression[, selected_samples] <- round(selected_gene_expression[, selected_samples], 3)  # Round the expression values
          
          # Order the selected_gene_expression based on the selected_samples column
          selected_gene_expression <- selected_gene_expression[order(selected_gene_expression[, selected_samples], decreasing = TRUE), ]  # Order the gene expression data based on the selected sample column
        } else {
          if (num_samples == 2) {  # If there are exactly 2 selected samples
            # Compute mean for each gene across two samples
            selected_gene_expression$mean <- round(rowMeans(selected_gene_expression[, selected_samples], na.rm = TRUE), 3)  # Compute the mean expression across samples for each gene
            selected_gene_expression$"FC (Sample 1 vs Sample 2)" <- selected_gene_expression[, selected_samples[1]] - selected_gene_expression[, selected_samples[2]]
            selected_gene_expression <- selected_gene_expression[order(selected_gene_expression$mean, decreasing = TRUE), ]  # Order the gene expression data based on the mean expression column
          } else {
            # Compute standard deviation for each gene across three or more samples
            selected_gene_expression$std_dev <- round(apply(selected_gene_expression[, selected_samples], 1, sd, na.rm = TRUE), 3)  # Compute the standard deviation across samples for each gene
            selected_gene_expression <- selected_gene_expression[order(selected_gene_expression$std_dev, decreasing = TRUE), ]  # Order the gene expression data based on the standard deviation column
          }
        }

      
      
      
      # Select the top and bottom expression values based on the total number of rows (genes)
      ##num_genes <- if (input$slider_type == "Gene Expression Range For At Least 1 Sample") {
      num_genes <- if (input$slider_type == "Gene Expression Range") {
        req(selected_gene_expression)  # Ensure selected_gene_expression is available
        nrow(selected_gene_expression)
      } else {
        req(input$num_genes_displayed)  # Ensure input$num_genes_displayed is available
        if(nrow(selected_gene_expression) < input$num_genes_displayed) {
          nrow(selected_gene_expression) 
        }
        else{
          input$num_genes_displayed 
        }
      }
      
      
      # Determine the number of top and bottom genes based on the even value
      num_top_genes <- num_genes %/% 2  # Use integer division to get the floor value
      num_bottom_genes <- num_genes - num_top_genes
      
      # If num_genes is odd, adjust the values to make both num_top_genes and num_bottom_genes cover the entire range of genes
      if (num_genes %% 2 == 1) {
        num_top_genes <- num_genes %/% 2 + 1
        num_bottom_genes <- num_genes %/% 2
      }
      
            
      if(num_top_genes <= nrow(selected_gene_expression) | gene_filter == "Biomarkers") {
        if (num_top_genes == 1 & num_bottom_genes == 0) {
          top_bottom_genes <- selected_gene_expression[1:num_top_genes, ]
        }
        else {
          top_bottom_genes <- rbind(
            selected_gene_expression[1:num_top_genes, ],
            selected_gene_expression[(nrow(selected_gene_expression) - num_bottom_genes + 1):nrow(selected_gene_expression), ] )
        }
        
      }
      else{
        top_bottom_genes <- selected_gene_expression
        
      } 

    
      top_bottom_genes <- top_bottom_genes %>%
        mutate_if(is.numeric, round, digits = 3)  # Round numeric values to 3 decimal places
      
      rownames(top_bottom_genes) <- NULL
      #print(min(top_bottom_genes[, selected_samples]))
      #print(max(top_bottom_genes[, selected_samples]))
      
      
      #return(top_bottom_genes)
      return(list(my_data = top_bottom_genes, sample_list = selected_samples))
    }
    
  }
)
  
  filter_gene_expression <- reactive({
    my_gene_expression <- my_gene_expression()
    gene_data <- my_gene_expression$my_data
    gene_samples <- my_gene_expression$sample_list
  
    
    ## if (input$slider_type == "Gene Expression Range For At Least 1 Sample") {
    if (input$slider_type == "Gene Expression Range") {
      req(input$gene_expr_range)
      #print("Lower Bound")
      #print(input$gene_expr_range[1])
      #print("Upper Bound")
      #print(input$gene_expr_range[2])
      if (length(gene_samples) == 1) {
        sample_col <- gene_samples[1]
        gene_data <- gene_data[gene_data[, sample_col] >= (input$gene_expr_range[1] - 0.01) &
                                 gene_data[, sample_col] <= (input$gene_expr_range[2] + 0.01), ]
        
      } else {
        samples_in_range <- rowSums(gene_data[, gene_samples] >= (input$gene_expr_range[1] - 0.01) &
                                      gene_data[, gene_samples] <= (input$gene_expr_range[2] + 0.01), na.rm = TRUE) > 0
        gene_data <- gene_data[samples_in_range, ]
      }
    }
    
  
  
    return(gene_data)
    
    
  })

  output$gene_expression_data <- DT::renderDataTable({
    datatable(
      filter = "top",  style='bootstrap' ,
      filter_gene_expression(),
      extensions = 'Buttons',
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = c(0, 2, 3, 4, 5, 6, 7, 8))  # Hides the first three columns
        ),
        dom = 'Bfrtip',  # Displays the column visibility button
        buttons = list('colvis'),  # Adds the column visibility button
        ## pageLength = input$pageLengthInputGene
        pageLength = 10, 
        lengthMenu = c(10, 25, 50, 100)
      )
    )
  })
  
  # MUTATION DATA
  # Display patientID
  output$mutation_info <- renderUI({
    # req(selected_row())  # Ensure that a row is selected
    req(selected_pat())
    tags$div(
      style = "font-size: 16px; color: #333333;",
      tags$p(
        style = "font-size: 20px; color: blue;",
        ## paste("Patient ID: ",  PatientSearchData[selected_row(), "patientID"])
        paste("Patient ID: ",  PatientSearchData[which(PatientSearchData$patientID==selected_pat()), "patientID"])
        
      )
    )
  })
  
  # Reset input values when a new patient is selected
  observeEvent(selected_row(), {
    updateRadioButtons(session, "mutation_filter", selected = "All Mutations")  # Reset the mutation filter radio buttons to default
    updateRadioButtons(session, "slider_type_mutation", selected = "Top and Bottom Mutations")  # Reset the slider type radio buttons to default
    #updateSliderInput(session, "mut_range", value = c(0, 100))  # Reset the mutation slider to default values
    updateSliderInput(session, "num_mut_displayed", value = 10)  # Reset the number of mutations displayed slider to default value
  })
  
  output$mutation_filter_label <- renderUI({
    ## patientDisease <- PatientSearchData[selected_row(), "disease"]
    patientDisease <- PatientSearchData[which(PatientSearchData$patientID==selected_pat()), "disease"]
    disease_label <- paste("Disease: ", patientDisease)
    tagList(
      tags$label(disease_label, style = "font-weight: bold;"),
      br()
    )
  })
  
#Mutation Data
  my_mutation <- reactive ({
    # Get the selected patient's sample column names
    # if (!is.null(selected_row())) {  # Check if a row is selected
    #   patientID <- PatientSearchData[selected_row(), "patientID"]  # Get the patient ID of the selected row
    if (!is.null(selected_pat())) {  # Check if a row is selected
        patientID <- selected_pat()  # Get the patient ID of the selected row
        
      
      selected_samples <- PatientSampleData$Sample_Name[PatientSampleData$patientID == patientID]  # Get the selected patient's sample names
      mutation_cols <- c("Symbol", "Chromosome", "Start", "End", "Strand", selected_samples)  # Define the columns to select from gene expression data
      selected_mutation <- combined_mutation[, mutation_cols, drop = FALSE]  # Select the specified columns from the combined gene expression data
      sample_columns <- selected_mutation[, selected_samples, drop = FALSE] %>%
        select_if(~ !all(is.na(.)))
      shiny::validate(need(ncol(sample_columns)>0,"No mutation data for this patient")) ## FE new
      selected_samples <- colnames(sample_columns)
      
      mutation_cols <- c("Symbol", "Chromosome", "Start", "End", "Strand", selected_samples)
      selected_mutation <- selected_mutation[, mutation_cols, drop = FALSE]
      
      
      if (length(selected_samples) == 1) {
        selected_mutation <- selected_mutation[selected_mutation[, selected_samples] != 0, ]
      } else {
        selected_mutation <- selected_mutation[rowSums(selected_mutation[, selected_samples, drop = FALSE] != 0, na.rm = TRUE) > 0, ]
      }
      
      
      # Get the patient's disease
      ## patientDisease <- PatientSearchData[selected_row(), "disease"]
      patientDisease <- PatientSearchData[which(PatientSearchData$patientID==selected_pat()), "disease"]
      
      # Apply gene filter
      mutation_filter <- input$mutation_filter  # Get the selected gene filter option
      
      # Read the biomarkers JSON file
      biomarkers <- fromJSON("biomarkers.json")
      
      if ("Biomarkers" %in% mutation_filter) {  # If "Biomarkers" option is selected
        if (patientDisease %in% names(biomarkers)) {
          filtered_mutation <- selected_mutation[selected_mutation$symbol %in% biomarkers[[patientDisease]][["mutation"]], ]
          if(nrow(filtered_mutation) == 0){
            shiny::validate(
              need(FALSE, "There are no non-zero biomarkers available for this sample. Please return to the all mutations filter")
            )
          }
        } else {
          #Display shinyalert with "No biomarkers found" message
          shinyalert::shinyalert(
            title = "No Biomarkers Found",
            text = paste("We currently do not have mutation biomarkers data for the disease:", patientDisease, ". Please return to the 'All Mutations' filter"),
            type = "info"
          ) # Use the original selected gene expression data without any filter
          shiny::validate(
            need(FALSE, "No biomarkers found. Please return to the all mutations filter.")
          )
        }
      } else {
        filtered_mutation <- selected_mutation  # Use the original selected gene expression data without any filter
      }
      
      
      # Compute mean and standard deviation columns
      num_samples <- length(selected_samples)  # Get the number of selected samples
      
      if (num_samples >= 2) {  # If there are at least 2 samples
        filtered_mutation$mean <- round(rowMeans(filtered_mutation[, selected_samples], na.rm = TRUE), 3)  # Compute the mean expression across samples for each gene
      }
      
      if (num_samples >= 3) {  # If there are at least 3 samples
        filtered_mutation$std_dev <- round(apply(filtered_mutation[, selected_samples], 1, sd, na.rm = TRUE), 3)  # Compute the standard deviation across samples for each gene
      }
      

        selected_mutation <- filtered_mutation
        
        if (num_samples == 1) {  # If there is only 1 selected sample
          # Handle the case when num_samples == 1
          # Only default columns and selected sample
          columns_to_display <- c("Symbol", "Chromosome", "Start", "End", "Strand", selected_samples)
          selected_mutation <- selected_mutation[, columns_to_display, drop = FALSE]  # Select the specified columns
          selected_mutation[, selected_samples] <- round(selected_mutation[, selected_samples], 3)  # Round the expression values
          
          # Order the selected_gene_expression based on the selected_samples column
          selected_mutation <- selected_mutation[order(selected_mutation[, selected_samples], decreasing = TRUE), ]  # Order the gene expression data based on the selected sample column
        } else {
          if (num_samples == 2) {  # If there are exactly 2 selected samples
            # Compute mean for each gene across two samples
            selected_mutation$mean <- round(rowMeans(selected_mutation[, selected_samples], na.rm = TRUE), 3)  # Compute the mean expression across samples for each gene
            #selected_mutation$"FC (Sample 1 vs Sample 2)" <- selected_mutation[, selected_samples[1]] - selected_mutation[, selected_samples[2]]
            selected_mutation <- selected_mutation[order(selected_mutation$mean, decreasing = TRUE), ]  # Order the gene expression data based on the mean expression column
          } else {
            # Compute standard deviation for each gene across three or more samples
            if (ncol(selected_mutation[, selected_samples]) > 1) {
              selected_mutation$std_dev <- round(apply(selected_mutation[, selected_samples], 1, sd, na.rm = TRUE), 3)  # Compute the standard deviation across samples for each gene
              selected_mutation <- selected_mutation[order(selected_mutation$std_dev, decreasing = TRUE), ]  # Order the gene expression data based on the standard deviation column
            }
            
          }
        }
      #}
      
      # Select the top and bottom expression values based on the total number of rows (genes)
      
      num_mutation <- if (input$slider_type_mutation == "Gene Mutation Range For At Least 1 Sample") {
        req(selected_mutation)  # Ensure selected_gene_expression is available
        nrow(selected_mutation)
      } else {
        req(input$num_mut_displayed)  # Ensure input$num_genes_displayed is available
        if(nrow(selected_mutation) < input$num_mut_displayed) {
          nrow(selected_mutation) 
        }
         else{
           input$num_mut_displayed 
         }
      }
      
      
      # if ((input$num_mut_displayed > nrow(selected_mutation) & (input$mutation_filter == "Biomarkers"))) {
      #   num_mutation <- nrow(selected_mutation)  # Set num_genes to the total number of genes if it exceeds the available biomarkers
      #   shinyalert::shinyalert(
      #     title = "Limited Number of Biomarkers",
      #     text = paste("There are only", num_mutation, "biomarkers available."),
      #     type = "info"
      #   )
      # }
      
      # Determine the number of top and bottom genes based on the even value
      num_top_mutation <- num_mutation %/% 2  # Use integer division to get the floor value
      num_bottom_mutation <- num_mutation - num_top_mutation
      
      # If num_genes is odd, adjust the values to make both num_top_genes and num_bottom_genes cover the entire range of genes
      if (num_mutation %% 2 == 1) {
        num_top_mutation <- num_mutation %/% 2 + 1
        num_bottom_mutation <- num_mutation %/% 2
      }
      
      
      
      if(num_top_mutation <= nrow(selected_mutation) | mutation_filter == "Biomarkers") {
        if ((num_top_mutation == 1 & num_bottom_mutation == 0)) {
          top_bottom_mutation <- selected_mutation[1:num_top_mutation, ]
        }
        else {
          top_bottom_mutation <- rbind(
            selected_mutation[1:num_top_mutation, ],
            selected_mutation[(nrow(selected_mutation) - num_bottom_mutation + 1):nrow(selected_mutation), ] )
        }
        
      }
      else{
        top_bottom_mutation <- selected_mutation
        
      }
      # Remove NA entries and round numeric values
      #top_bottom_mutation <-  top_bottom_mutation[complete.cases(top_bottom_mutation), ] %>%  # Remove rows with NA values
      top_bottom_mutation <- top_bottom_mutation %>%
        mutate_if(is.numeric, round, digits = 3)  
        
      
      if (all(is.na(top_bottom_mutation[, selected_samples]))) {
        shiny::validate(
          need(FALSE, "We do not have mutation data for this sample")
        )
      }
      
      
      rownames(top_bottom_mutation) <- NULL
      
      return(list(my_data_mutation = top_bottom_mutation, sample_list_mutation = selected_samples))
      
    }
    
  })

  
  filter_mutation <- reactive({
    my_mutation_data <- my_mutation()
    mutation_data_display <- my_mutation_data$my_data_mutation
    mutation_samples <- my_mutation_data$sample_list_mutation
    
    
    if (input$slider_type_mutation == "Gene Mutation Range For At Least 1 Sample") {
      req(input$mut_range)
      if (length(mutation_samples) == 1) {
        sample_col_mutation <- mutation_samples[1]
        mutation_data_display <- mutation_data_display[mutation_data_display[, sample_col_mutation] >= (input$mut_range[1] - 0.01) &
                                 mutation_data_display[, sample_col_mutation] <= (input$mut_range[2] + 0.01), ]
        
      } else {
        samples_in_range_mut <- rowSums(mutation_data_display[, mutation_samples] >= (input$mut_range[1] - 0.01) &
                                      mutation_data_display[, mutation_samples] <= (input$mut_range[2] + 0.01), na.rm = TRUE) > 0
        mutation_data_display <- mutation_data_display[samples_in_range_mut, ]
      }
    }
    
    
    # cat(dim(mutation_data_display)," mutation table dimension \n")
    # cat(colnames(mutation_data_display)," mutation table names \n")
    return(mutation_data_display)
    
    
  })
  
  
  output$mutation_data <- DT::renderDataTable({
    datatable(
      filter_mutation(),
      filter = "top", style='bootstrap' ,
      extensions = 'Buttons',
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = c(0, 2, 3, 4, 5)) # Hides the first three columns
        ),
        dom = 'Bfrtip',  # Displays the column visibility button
        buttons = list('colvis'),  # Adds the column visibility button
        #pageLength = input$pageLengthInputMutation, 
        pageLength = 10, 
        lengthMenu = c(10, 25, 50, 100)
      )
    )
  })
  
  #DISPLAY METHYLATION DATA --------------------------------------
  # Display patientID
  output$methylation_info <- renderUI({
    req(selected_row())  # Ensure that a row is selected
    tags$div(
      style = "font-size: 16px; color: #333333;",
      tags$p(
        style = "font-size: 20px; color: blue;",
        paste("Patient ID: ",  PatientSearchData[selected_row(), "patientID"])
      )
    )
  })
  
  # Reset input values when a new patient is selected
  observeEvent(selected_row(), {
    updateRadioButtons(session, "methylation_filter", selected = "All Genes")  # Reset the gene filter radio buttons to default
    updateRadioButtons(session, "slider_type_methylation", selected = "Top and Bottom Methylation Scores")  # Reset the slider type radio buttons to default
    updateSliderInput(session, "methylation_range", value = c(0, 1))  # Reset the gene expression range slider to default values
    updateSliderInput(session, "num_methylation_displayed", value = 10)  # Reset the number of genes displayed slider to default value
  })
  
  output$methylation_filter_label <- renderUI({
    patientDisease <- PatientSearchData[selected_row(), "disease"]
    disease_label <- paste("Disease: ", patientDisease)
    tagList(
      tags$label(disease_label, style = "font-weight: bold;"),
      br()
    )
  })
  
  my_methylation <- reactive ({
    # Get the selected patient's sample column names
    if (!is.null(selected_row())) {  # Check if a row is selected
      shiny::validate(need(!is.null(hmethyl),"No methylation data for this patient")) ## FE new 02/25
      
      patientID <- PatientSearchData[selected_row(), "patientID"]  # Get the patient ID of the selected row
      selected_samples <- PatientSampleData$Sample_Name[PatientSampleData$patientID == patientID]  # Get the selected patient's sample names
      ## methylation_cols <- c("gene", "accession", "priority", "probe.count", selected_samples)  # Define the columns to select from gene expression data
      
      methylation_cols <- c("symbol", "chr", "start", "end", selected_samples)  # Define the columns to select from gene expression data
      selected_methylation <- combined_methylation[, methylation_cols, drop = FALSE]  # Select the specified columns from the combined gene expression data
      
      sample_columns <- selected_methylation[, selected_samples, drop = FALSE] %>%
        select_if(~ !all(is.na(.)))
      shiny::validate(need(ncol(sample_columns)>0,"No methylation data for this patient")) ## FE new
      selected_samples <- colnames(sample_columns)
      ## methylation_cols <- c("gene", "accession", "priority", "probe.count", selected_samples)  ## FE to add
      
      methylation_cols <- c("symbol", "chr", "start", "end", selected_samples)  ## FE to add
      selected_methylation <- selected_methylation[, methylation_cols, drop = FALSE]
      
      # Get the patient's disease
      patientDisease <- PatientSearchData[selected_row(), "disease"]
      
      # Apply gene filter
      methylation_filter <- input$methylation_filter  # Get the selected gene filter option
      
      # Access the gene_expression lists for ACC and SCLC
      
      if ("Biomarkers" %in% methylation_filter) {  # If "Biomarkers" option is selected
        if (patientDisease %in% names(biomarkers)) {
          filtered_methylation <- selected_methylation[selected_methylation$gene %in% biomarkers[[patientDisease]][["methylation"]], ]
          if(nrow(filtered_methylation) == 0){
            shiny::validate(
              need(FALSE, "There are no biomarkers available for this sample. Please return to the all genes filter")
            ) }
        } else {
          #Display shinyalert with "No biomarkers found" message
          shinyalert::shinyalert(
            title = "No Biomarkers Found",
            text = paste("We currently do not have methylation biomarkers data for the disease:", patientDisease, ". Please return to the 'All Genes' filter"),
            type = "info"
          ) # Use the original selected gene expression data without any filter
          shiny::validate(
            need(FALSE, "No biomarkers found. Please return to the all genes filter.")
          )
        }
      } else {
        filtered_methylation<- selected_methylation  # Use the original selected gene expression data without any filter
      }
      
      
      # Compute mean and standard deviation columns
      num_samples <- length(selected_samples)  # Get the number of selected samples
      
      if (num_samples >= 2) {  # If there are at least 2 samples
        filtered_methylation$mean <- round(rowMeans(filtered_methylation[, selected_samples], na.rm = TRUE), 3)  # Compute the mean expression across samples for each gene
      }
      
      if (num_samples >= 3) {  # If there are at least 3 samples
        filtered_methylation$std_dev <- round(apply(filtered_methylation[, selected_samples], 1, sd, na.rm = TRUE), 3)  # Compute the standard deviation across samples for each gene
      }
      
        selected_methylation <- filtered_methylation
        
        if (num_samples == 1) {  # If there is only 1 selected sample
          # Handle the case when num_samples == 1
          # Only default columns and selected sample
          ## columns_to_display <- c("gene", "accession", "priority", "probe.count", selected_samples)
          
          columns_to_display <- c("symbol", "chr", "start", "end", selected_samples)
          selected_methylation <- selected_methylation[, columns_to_display, drop = FALSE]  # Select the specified columns
          selected_methylation[, selected_samples] <- round(selected_methylation[, selected_samples], 3)  # Round the expression values
          
          # Order the selected_gene_expression based on the selected_samples column
          selected_methylation <- selected_methylation[order(selected_methylation[, selected_samples], decreasing = TRUE), ]  # Order the gene expression data based on the selected sample column
        } else {
          if (num_samples == 2) {  # If there are exactly 2 selected samples
            # Compute mean for each gene across two samples
            selected_methylation$mean <- round(rowMeans(selected_methylation[, selected_samples], na.rm = TRUE), 3)  # Compute the mean expression across samples for each gene
            #selected_methylation$"FC (Sample 1 vs Sample 2)" <- selected_methylation[, selected_samples[1]] - selected_methylation[, selected_samples[2]]
            selected_methylation <- selected_methylation[order(selected_methylation$mean, decreasing = TRUE), ]  # Order the gene expression data based on the mean expression column
          } else {
            # Compute standard deviation for each gene across three or more samples
            selected_methylation$std_dev <- round(apply(selected_methylation[, selected_samples], 1, sd, na.rm = TRUE), 3)  # Compute the standard deviation across samples for each gene
            selected_methylation <- selected_methylation[order(selected_methylation$std_dev, decreasing = TRUE), ]  # Order the gene expression data based on the standard deviation column
          }
        }
      #}
      
     selected_methylation <- selected_methylation[complete.cases(selected_methylation), ] ## why? FE
     
    # Select the top and bottom expression values based on the total number of rows (genes)
        num_methylation <- if (input$slider_type_methylation == "Methylation Range For At Least 1 Sample") {
          req(selected_methylation)  # Ensure selected_methylation is available
          nrow(selected_methylation)
        } else {
          req(input$num_methylation_displayed)  # Ensure input$num_num_methylation_displayed is available
          if(nrow(selected_methylation) < input$num_methylation_displayed) {
            nrow(selected_methylation) 
          }
          else{
            input$num_methylation_displayed
          }
        }
        
  
        # Determine the number of top and bottom genes based on the even value
        num_top_methylation <- num_methylation %/% 2  # Use integer division to get the floor value
        num_bottom_methylation <- num_methylation - num_top_methylation
        
        # If num_genes is odd, adjust the values to make both num_top_genes and num_bottom_genes cover the entire range of genes
        if (num_methylation %% 2 == 1) {
          num_top_methylation <- num_methylation %/% 2 + 1
          num_bottom_methylation <- num_methylation %/% 2
        }
        
        
        # Print the values for debugging
        
        # Combine the top and bottom genes into a single dataframe
        if (num_top_methylation <= nrow(selected_methylation) | methylation_filter == "Biomarkers") {
          if (num_top_methylation == 1 & num_bottom_methylation == 0) {
            top_bottom_methylation <- selected_methylation[1:num_top_methylation, ]
          } else {
            top_bottom_methylation <- rbind(
              selected_methylation[1:num_top_methylation, ],
              selected_methylation[(nrow(selected_methylation) - num_bottom_methylation + 1):nrow(selected_methylation), ]
            )
          }
        } else {
          top_bottom_methylation <- selected_methylation
        }
         
          
        # Remove NA entries and round numeric values
        top_bottom_methylation <- top_bottom_methylation %>%
          mutate_if(is.numeric, round, digits = 3)  

        
        if (all(is.na(top_bottom_methylation[, selected_samples]))) {
          shiny::validate(
            need(FALSE, "We do not have methylation data for this sample")
          )
        }

        # top_bottom_methylation <- top_bottom_methylation %>%
        #   select_if(~ !all(is.na(.))) %>%
        #   distinct()  # Remove duplicate rows
    
      
      
      rownames(top_bottom_methylation) <- NULL
      
      return(list(my_data_methylation = top_bottom_methylation, sample_list_methylation = selected_samples))
      
    }
   
    
  
})
  
  filter_methylation <- reactive({
    my_methylation_data <- my_methylation()
    methylation_data_display <- my_methylation_data$my_data_methylation
    methylation_samples <- my_methylation_data$sample_list_methylation
    
   # cat(" dim1 m: ",dim(methylation_data_display),"\n")
    
    if (input$slider_type_methylation == "Methylation Range For At Least 1 Sample") {  # If the slider type is "Gene Mutation Range For At Least 1 Sample"
      if (length(methylation_samples) == 1) {
        methylation_data_display <- methylation_data_display[methylation_data_display[, methylation_samples[1]] >= input$methylation_range[1] &
                                                               methylation_data_display[, methylation_samples[1]] <= input$methylation_range[2], ]
      } else {
        methylation_data_display <- methylation_data_display[rowSums(methylation_data_display[, methylation_samples] >= input$methylation_range[1] &
                                                                       methylation_data_display[, methylation_samples] <= input$methylation_range[2], na.rm = TRUE) > 0, ]
      }
      
    }
    
   # cat(" dim2 m: ",dim(methylation_data_display),"\n")
    return(methylation_data_display)
    
    
  })
 
  
  output$methylation_data <- DT::renderDataTable({
    datatable(
      filter_methylation(),
      filter = "top" , style='bootstrap' ,
      extensions = 'Buttons',
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = c(0, 2, 3, 4)) # Hides the first three columns
        ),
        dom = 'Bfrtip',  # Displays the column visibility button
        buttons = list('colvis'),  # Adds the column visibility button
        # pageLength = input$pageLengthInputMethylation
        pageLength = 10, 
        lengthMenu = c(10, 25, 50, 100)
      )
    )
  })
  
  #DISPLAY COPY NUMBER DATA --------------------------------------
 
  output$copy_info <- renderUI({
    req(selected_row())  # Ensure that a row is selected
    tags$div(
      style = "font-size: 16px; color: #333333;",
      tags$p(
        style = "font-size: 20px; color: blue;",
        paste("Patient ID: ",  PatientSearchData[selected_row(), "patientID"])
      )
    )
  })
 
  output$copy_filter_label <- renderUI({
    patientDisease <- PatientSearchData[selected_row(), "disease"]
    disease_label <- paste("Disease: ", patientDisease)
    tagList(
      tags$label(disease_label, style = "font-weight: bold;"),
      br()
    )
  })
  
  my_copy <- reactive ({
    # Get the selected patient's sample column names
    if (!is.null(selected_row())) {  # Check if a row is selected
      shiny::validate(need(!is.null(hcopy),"No copy number data for this patient")) ## FE new 02/25
      
      patientID <- PatientSearchData[selected_row(), "patientID"]  # Get the patient ID of the selected row
      selected_samples <- PatientSampleData$Sample_Name[PatientSampleData$patientID == patientID]  # Get the selected patient's sample names
      ## methylation_cols <- c("gene", "accession", "priority", "probe.count", selected_samples)  # Define the columns to select from gene expression data
      
      copy_cols <- c("gene", "chr", "start", "end", selected_samples)  # Define the columns to select from gene expression data
      selected_copy <- combined_copy[, copy_cols, drop = FALSE]  # Select the specified columns from the combined gene expression data
      
      sample_columns <- selected_copy[, selected_samples, drop = FALSE] %>%
        select_if(~ !all(is.na(.)))
      shiny::validate(need(ncol(sample_columns)>0,"No copy data for this patient")) ## FE new
      selected_samples <- colnames(sample_columns)
      ## copy_cols <- c("gene", "accession", "priority", "probe.count", selected_samples)  ## FE to add
      
      copy_cols <- c("gene", "chr", "start", "end", selected_samples)  ## FE to add
      selected_copy <- selected_copy[, copy_cols, drop = FALSE]
      
      # Get the patient's disease
      patientDisease <- PatientSearchData[selected_row(), "disease"]
      
      # Apply gene filter
      copy_filter <- input$copy_filter  # Get the selected gene filter option
      
      # Access the gene_expression lists for ACC and SCLC
      
      if ("Biomarkers" %in% copy_filter) {  # If "Biomarkers" option is selected
        if (patientDisease %in% names(biomarkers)) {
          filtered_copy <- selected_copy[selected_copy$gene %in% biomarkers[[patientDisease]][["copy_number"]], ]
          if(nrow(filtered_copy) == 0){
            shiny::validate(
              need(FALSE, "There are no biomarkers available for this sample. Please return to the all genes filter")
            ) }
        } else {
          #Display shinyalert with "No biomarkers found" message
          shinyalert::shinyalert(
            title = "No Biomarkers Found",
            text = paste("We currently do not have copy biomarkers data for the disease:", patientDisease, ". Please return to the 'All Genes' filter"),
            type = "info"
          ) # Use the original selected gene expression data without any filter
          shiny::validate(
            need(FALSE, "No biomarkers found. Please return to the all genes filter.")
          )
        }
      } else {
        filtered_copy<- selected_copy  # Use the original selected gene expression data without any filter
      }
      
      
      # Compute mean and standard deviation columns
      num_samples <- length(selected_samples)  # Get the number of selected samples
      
      if (num_samples >= 2) {  # If there are at least 2 samples
        filtered_copy$mean <- round(rowMeans(filtered_copy[, selected_samples], na.rm = TRUE), 3)  # Compute the mean expression across samples for each gene
      }
      
      if (num_samples >= 3) {  # If there are at least 3 samples
        filtered_copy$std_dev <- round(apply(filtered_copy[, selected_samples], 1, sd, na.rm = TRUE), 3)  # Compute the standard deviation across samples for each gene
      }
      
      selected_copy <- filtered_copy
      
      if (num_samples == 1) {  # If there is only 1 selected sample
        # Handle the case when num_samples == 1
        # Only default columns and selected sample
        ## columns_to_display <- c("gene", "accession", "priority", "probe.count", selected_samples)
        
        columns_to_display <- c("gene", "chr", "start", "end", selected_samples)
        selected_copy <- selected_copy[, columns_to_display, drop = FALSE]  # Select the specified columns
        selected_copy[, selected_samples] <- round(selected_copy[, selected_samples], 3)  # Round the expression values
        
        # Order the selected_gene_expression based on the selected_samples column
        selected_copy <- selected_copy[order(selected_copy[, selected_samples], decreasing = TRUE), ]  # Order the gene expression data based on the selected sample column
      } else {
        if (num_samples == 2) {  # If there are exactly 2 selected samples
          # Compute mean for each gene across two samples
          selected_copy$mean <- round(rowMeans(selected_copy[, selected_samples], na.rm = TRUE), 3)  # Compute the mean expression across samples for each gene
          #selected_copy$"FC (Sample 1 vs Sample 2)" <- selected_copy[, selected_samples[1]] - selected_copy[, selected_samples[2]]
          selected_copy <- selected_copy[order(selected_copy$mean, decreasing = TRUE), ]  # Order the gene expression data based on the mean expression column
        } else {
          # Compute standard deviation for each gene across three or more samples
          selected_copy$std_dev <- round(apply(selected_copy[, selected_samples], 1, sd, na.rm = TRUE), 3)  # Compute the standard deviation across samples for each gene
          selected_copy <- selected_copy[order(selected_copy$std_dev, decreasing = TRUE), ]  # Order the gene expression data based on the standard deviation column
        }
      }
      #}
      
      selected_copy <- selected_copy[complete.cases(selected_copy), ] ## why ? FE
      
      # Select the top and bottom expression values based on the total number of rows (genes)
      num_copy <- if (input$slider_type_copy == "Copy Range For At Least 1 Sample") {
        req(selected_copy)  # Ensure selected_copy is available
        nrow(selected_copy)
      } else {
        req(input$num_copy_displayed)  # Ensure input$num_num_copy_displayed is available
        if(nrow(selected_copy) < input$num_copy_displayed) {
          nrow(selected_copy) 
        }
        else{
          input$num_copy_displayed
        }
      }
      
      
      # Determine the number of top and bottom genes based on the even value
      num_top_copy <- num_copy %/% 2  # Use integer division to get the floor value
      num_bottom_copy <- num_copy - num_top_copy
      
      # If num_genes is odd, adjust the values to make both num_top_genes and num_bottom_genes cover the entire range of genes
      if (num_copy %% 2 == 1) {
        num_top_copy <- num_copy %/% 2 + 1
        num_bottom_copy <- num_copy %/% 2
      }
      
      
      # Print the values for debugging
      
      # Combine the top and bottom genes into a single dataframe
      if (num_top_copy <= nrow(selected_copy) | copy_filter == "Biomarkers") {
        if (num_top_copy == 1 & num_bottom_copy == 0) {
          top_bottom_copy <- selected_copy[1:num_top_copy, ]
        } else {
          top_bottom_copy <- rbind(
            selected_copy[1:num_top_copy, ],
            selected_copy[(nrow(selected_copy) - num_bottom_copy + 1):nrow(selected_copy), ]
          )
        }
      } else {
        top_bottom_copy <- selected_copy
      }
      
      
      # Remove NA entries and round numeric values
      top_bottom_copy <- top_bottom_copy %>%
        mutate_if(is.numeric, round, digits = 3)  
      
      
      if (all(is.na(top_bottom_copy[, selected_samples]))) {
        shiny::validate(
          need(FALSE, "We do not have copy data for this sample")
        )
      }
      
      # top_bottom_copy <- top_bottom_copy %>%
      #   select_if(~ !all(is.na(.))) %>%
      #   distinct()  # Remove duplicate rows
      
     # cat(range(top_bottom_copy[,selected_samples])," range copy \n")
     # write.csv(top_bottom_copy,"copy_number_selected.csv")
      
      rownames(top_bottom_copy) <- NULL
      
      return(list(my_data_copy = top_bottom_copy, sample_list_copy = selected_samples))
      
    }
    
    
    
  })
  
  filter_copy <- reactive({
    my_copy_data <- my_copy()
    copy_data_display <- my_copy_data$my_data_copy
    copy_samples <- my_copy_data$sample_list_copy
    
   # cat(" dim1: ",dim(copy_data_display),"\n")
   #  cat(" samples: ",copy_samples,"\n")
    
    if (input$slider_type_copy == "Copy Range For At Least 1 Sample") {  # If the slider type is "Gene Mutation Range For At Least 1 Sample"
      if (length(copy_samples) == 1) {
       # cat(copy_samples[1], ":copy sample 1 \n")
        copy_s <- copy_samples[1]
        copy_data_display <- copy_data_display[copy_data_display[, copy_s] >= input$copy_range[1] &
                                                               copy_data_display[, copy_s] <= input$copy_range[2], ]
      } else {
        copy_data_display <- copy_data_display[rowSums(copy_data_display[, copy_samples] >= input$copy_range[1] &
                                                                       copy_data_display[, copy_samples] <= input$copy_range[2], na.rm = TRUE) > 0, ]
      }
      
    }
   # cat(" dim2: ",dim(copy_data_display),"\n")
    
    
    return(copy_data_display)
    
    
  })
  
  output$copy_data <- DT::renderDataTable({
    datatable(
      filter_copy(),
      filter = "top" , style='bootstrap' ,
      extensions = 'Buttons',
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = c(0, 2, 3, 4)) # Hides the first three columns
        ),
        dom = 'Bfrtip',  # Displays the column visibility button
        buttons = list('colvis'),  # Adds the column visibility button
        # pageLength = input$pageLengthInputMethylation
        pageLength = 10, 
        lengthMenu = c(10, 25, 50, 100)
      )
    )
  })
  
  ## specific genes ***

  output$spg_info <- renderUI({
    req(selected_row())  # Ensure that a row is selected
    tags$div(
      style = "font-size: 16px; color: #333333;",
      tags$p(
        style = "font-size: 20px; color: blue;",
        paste("Patient ID: ",  PatientSearchData[selected_row(), "patientID"])
      )
    )
  })
  
  output$spg_filter_label <- renderUI({
    ## patientDisease <- PatientSearchData[selected_row(), "disease"]
    patientDisease <- PatientSearchData[which(PatientSearchData$patientID == selected_pat()), "disease"]
    disease_label <- paste("Disease: ", patientDisease)
    tagList(
      tags$label(disease_label, style = "font-weight: bold;"),
      br()
    )
  })
  
  my_spg <- reactive ({
    # Get the selected patient's sample column names
    # get exp, mut, met and cop for selected gene(s)
    
    # if (!is.null(selected_row())) {  # Check if a row is selected
    #   patientID <- PatientSearchData[selected_row(), "patientID"]  # Get the patient ID of the selected row
    #   patientDisease <- PatientSearchData[selected_row(), "disease"]

    if (!is.null(selected_pat())) {  # Check if a row is selected
      patientID <- selected_pat()
      patientDisease <- PatientSearchData[which(PatientSearchData$patientID == patientID), "disease"]
       
      # patientID <- PatientSearchData[selected_row(), "patientID"]  # Get the patient ID of the selected row
      # patientDisease <- PatientSearchData[selected_row(), "disease"]
      
      ## selected disease samples !!!!
      disease_samples <- PatientSampleData$Sample_Name[PatientSampleData$disease == patientDisease]
      ## end disease samples
      
      selected_samples <- PatientSampleData$Sample_Name[PatientSampleData$patientID == patientID]  # Get the selected patient's sample names
      
      # selgenes = trimws(input$spg_genes) ## for now 1 gene
      
      selgenes = unlist(strsplit(trimws(input$spg_genes)," ")) ## many genes
      
      
      df = data.frame(matrix( 
        vector(), 0, length(selected_samples)+5, dimnames=list(c(),  c(selected_samples,"median","min","max","gene","type"))), 
        stringsAsFactors=F) 
      
     selxsq = intersect(selgenes, rownames(xsq_a))
      selmut = intersect(selgenes, rownames(mutation_a)) 
      selmet = intersect(selgenes, rownames(methylation_a)) 
      selcop = intersect(selgenes, rownames(copy_a)) 
      
      gexp <- xsq_a[selxsq, selected_samples, drop = FALSE] 
      if (nrow(gexp)!=0) {
        gdisease <- xsq_a[selxsq, disease_samples, drop = FALSE] 
        vmin = apply(gdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(min(x, na.rm =T))  )
        vmax = apply(gdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(max(x, na.rm =T))  )
        vmedian = apply(gdisease, 1,function(x) median(x, na.rm =T))
        gexp$median.disease = vmedian; gexp$min.disease = as.numeric(vmin) ; gexp$max.disease = as.numeric(vmax)  
        gexp = round(gexp,3)
        gexp$gene = rownames(gexp)
        gexp$type = "xsq"
        df = rbind(df, gexp)
      }
      gmut <- mutation_a[selmut, selected_samples, drop = FALSE]  
      if (nrow(gmut)!=0) {
        mudisease <- mutation_a[selmut, disease_samples, drop = FALSE] 
        
        remain_columns <- gmut %>% select_if(~ !all(is.na(.)))
        if (ncol(remain_columns)>0)
        {
        vmin3 = apply(mudisease, 1, function(x) if (all(is.na(x))) return(NA) else return(min(x, na.rm =T))  )
        vmax3 = apply(mudisease, 1, function(x) if (all(is.na(x))) return(NA) else return(max(x, na.rm =T))  )
        vmedian3 = apply(mudisease, 1,function(x) median(x, na.rm =T))
        gmut$median.disease = vmedian3; gmut$min.disease = as.numeric(vmin3) ; gmut$max.disease = as.numeric(vmax3) 
        gmut = round(gmut,3)
        ## gmut$min.disease = NA; gmut$max.disease = NA; gmut$median.disease = NA
        gmut$gene = rownames(gmut)
        gmut$type = "mut"
        df = rbind(df, gmut)
        }
      }
      
      if (!is.null(hmethyl) )
      {
      gmet <- methylation_a[selmet, selected_samples, drop = FALSE] 
      if (nrow(gmet)!=0) {
        metdisease <- methylation_a[selmet, disease_samples, drop = FALSE] 
        vmin4 = apply(metdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(min(x, na.rm =T))  )
        vmax4 = apply(metdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(max(x, na.rm =T))  )
        vmedian4 = apply(metdisease, 1,function(x) median(x, na.rm =T))
        gmet$median.disease = vmedian4; gmet$min.disease = as.numeric(vmin4) ; gmet$max.disease = as.numeric(vmax4)  
        gmet = round(gmet,3)
        # gmet$min.disease = NA; gmet$max.disease = NA; gmet$median.disease = NA
        gmet$gene = rownames(gmet)
        gmet$type = "met"
        df = rbind(df, gmet)
      }
      }
      
      if (!is.null(hcopy) )
      {
      gcop <- copy_a[selcop, selected_samples, drop = FALSE] ; 
      if (nrow(gcop)!=0) {
        cpdisease <- copy_a[selcop, disease_samples, drop = FALSE] 
        vmin2 = apply(cpdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(min(x, na.rm =T)) )
        vmax2 = apply(cpdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(max(x, na.rm =T))  )
        vmedian2 = apply(cpdisease, 1,function(x) median(x, na.rm =T))
        gcop$median.disease = vmedian2; gcop$min.disease = as.numeric(vmin2); gcop$max.disease = as.numeric(vmax2)
        gcop = round(gcop,3)
        gcop$gene = rownames(gcop)
        gcop$type = "cop"
        df = rbind(df, gcop)
      }
      }
      
      shiny::validate(need(nrow(df)>0,"No data found for this gene"))
      nb = ncol(df)
      ## if (nb >= 3 ) 
      df = df[,c((nb-1), nb, 1:(nb-2))]
      df = df[order(df[,1]),]
      return(df)
      
    }
    
    
    
  })
  
  output$spg_data <- DT::renderDataTable({
    datatable(
      my_spg(), rownames = F,
      filter = "top" , style='bootstrap' ,
  #    extensions = 'Buttons',
      options = list(
      #   columnDefs = list(
      #     list(visible = FALSE, targets = c(0, 2, 3, 4)) # Hides the first three columns
      #   ),
        dom = 'frtip',  # Displays the column visibility button dom = 'Bfrtip'
 #       buttons = list('colvis'),  # Adds the column visibility button
        # pageLength = input$pageLengthInputMethylation
        pageLength = 10, 
        lengthMenu = c(10, 25, 50, 100)
      )
    )
  })
  
 ## end specific genes ***
  
## Immune Therapy Treatment *** --------------------------------
  output$trt_info <- renderUI({
    # req(selected_row())  # Ensure that a row is selected
    req(selected_pat())  # Ensure that a patient is selected
    
    tags$div(
      style = "font-size: 16px; color: #333333;",
      tags$p(
        style = "font-size: 20px; color: blue;",
        # paste("Patient ID: ",  PatientSearchData[selected_row(), "patientID"])
        paste("Patient ID: ",  selected_pat())
      )
    )
  })
  
  output$trt_filter_label <- renderUI({
    # patientDisease <- PatientSearchData[selected_row(), "disease"]
    patientDisease <- PatientSearchData[which(PatientSearchData$patientID ==selected_pat()), "disease"]
    disease_label <- paste("Disease: ", patientDisease)
    tagList(
      tags$label(disease_label, style = "font-weight: bold;"),
      br()
    )
  })
  
  my_trt <- reactive ({
    # Get the selected patient's sample column names
    # get exp, mut, met and cop for selected gene(s)
    
    # if (!is.null(selected_row())) {  # Check if a row is selected
    #   patientID <- PatientSearchData[selected_row(), "patientID"]  # Get the patient ID of the selected row
    #  patientDisease <- PatientSearchData[selected_row(), "disease"]
    cat(selected_pat(), " >> patient for trt entry...\n")
    if (!is.null(selected_pat())) {  # Check if a row is selected
      patientID <- selected_pat()  # Get the patient ID of the selected row
      patientDisease <- PatientSearchData[which(PatientSearchData$patientID == selected_pat()), "disease"]
        
      ## selected disease samples !!!!
      disease_samples <- PatientSampleData$Sample_Name[PatientSampleData$disease == patientDisease]
      ## end disease samples
      
      selected_samples <- PatientSampleData$Sample_Name[PatientSampleData$patientID == patientID]  # Get the selected patient's sample names
      
      # selgenes = trimws(input$spg_genes) ## for now 1 gene
      
      selgenes = rownames(adcgenes) ## 44 genes
      
      
      df = data.frame(matrix( 
        vector(), 0, length(selected_samples)+4, dimnames=list(c(),  c(selected_samples,"median", "min","max","gene"))), 
        stringsAsFactors=F) 
      
      selxsq = intersect(selgenes, rownames(xsq_a))
      
      ## gann <- adcgenes[selxsq,]
      
      gexp <- xsq_a[selxsq, selected_samples, drop = FALSE] 
      
      if (nrow(gexp)!=0) {
        ## sort by max samples values
        
        mmax = apply(gexp,1, function(x) if (all(is.na(x))) return(NA) else return(max(x, na.rm =T)))
        gexp$sample.max = mmax
        ## write.csv(gexp,"gexp_test.csv")
        ## cat(colnames(gexp), "\n")
        gexp = gexp[with(gexp,order(sample.max, decreasing = T)),]
        gexp$sample.max = NULL
        ##
        ## gdisease <- xsq_a[selxsq, disease_samples, drop = FALSE] 
        gdisease <- xsq_a[rownames(gexp), disease_samples, drop = FALSE] 
        ##
        ## write.csv(gdisease,"gdisease_test.csv")
        vmin = apply(gdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(min(x, na.rm =T))  )
        vmax = apply(gdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(max(x, na.rm =T))  )
        vmedian = apply(gdisease, 1,function(x) if (all(is.na(x))) return(NA) else return(median(x, na.rm =T)) )
        gexp$median.disease = vmedian; gexp$min.disease = as.numeric(vmin) ; gexp$max.disease = as.numeric(vmax)
        gexp = round(gexp,3)
        gexp$gene = rownames(gexp)
        # gexp$type = "xsq"
        df = rbind(df, gexp)
      }

      
      
      shiny::validate(need(nrow(df)>0,"No data found for this patient"))
      nb = ncol(df)
      ## if (nb >= 3 ) 
      df = df[,c(nb, 1:(nb-1))]
      
      ####   df = df[order(df[,1]),] ## sort by gene name
      
      ## filter for gene value above median
     
      if (length(selected_samples) == 1 )  { 
        # df = df[which(df[,selected_samples] >= df[,"median.disease"]),] 
        df = df[which(df[,selected_samples] >= df[,"median.disease"] & df[,selected_samples] >= 1),] 
       }
      else  {
        # genes_median <- rowSums(df[, selected_samples] >= df[,"median.disease"], na.rm=T) > 0
        genes_median <- rowSums(df[, selected_samples] >= df[,"median.disease"] & df[,selected_samples] >= 1 , na.rm=T) > 0
        df <- df[genes_median, ]
           }
      ## df = df[which(df[,selected_samples] >= df[,"median.disease"]),] ## if only one sample
      
      
      
      shiny::validate(need(nrow(df)>0,"No gene has value greater or equal to 1 and the median disease"))
      ## annotate with drug information
      gann <- adcgenes[df$gene,]
      df = cbind(gann,df)
      ## df = df[, c(5,1:4,6:ncol(df))]
      return(df)
      
    }
    
    
    
  })
  
  my_trt_v2 <- reactive ({
    # Get the selected patient's sample column names
    # get exp, mut, met and cop for selected gene(s)
    
    # if (!is.null(selected_row())) {  # Check if a row is selected
    #   patientID <- PatientSearchData[selected_row(), "patientID"]  # Get the patient ID of the selected row
    #   patientDisease <- PatientSearchData[selected_row(), "disease"]
    if (!is.null(selected_pat())) {  # Check if a row is selected
      patientID <- selected_pat()  # Get the patient ID of the selected row
      patientDisease <- PatientSearchData[which(PatientSearchData$patientID == patientID), "disease"]
        
      ## selected disease samples !!!!
      disease_samples <- PatientSampleData$Sample_Name[PatientSampleData$disease == patientDisease]
      ## end disease samples
      
      selected_samples <- PatientSampleData$Sample_Name[PatientSampleData$patientID == patientID]  # Get the selected patient's sample names
      
      selected_types <- PatientSampleData$sampleType[PatientSampleData$patientID == patientID]  # Get the selected patient's sample types
      
      inormal = which(selected_types == "normal")
      if (length(inormal) > 0) {
        selected_samples = selected_samples[-inormal]
        selected_types = selected_types[-inormal]
      }
      
      # selgenes = trimws(input$spg_genes) ## for now 1 gene
      
      ### selgenes = rownames(adcgenes) ## 43 genes
      selgenes = unique(adcgenes$gene) ## 43 genes
      
      
      
      selxsq = intersect(selgenes, rownames(xsq_a))
      
      ## gann <- adcgenes[selxsq,]
      
      gexp <- xsq_a[selxsq, selected_samples, drop = FALSE]  %>% select_if(~ !all(is.na(.)))
      cat("gexp dim 1: ", dim(gexp), " \n")
      
      shiny::validate(need(ncol(gexp)>0,"No expression data for this patient")) ## FE new
      
      if (ncol(gexp) != length(selected_samples)) {
         selected_types = selected_types[match(colnames(gexp), selected_samples)]
         selected_samples = colnames(gexp)
      }
      
      newnames = paste(selected_samples, selected_types, sep= " ")
      df = data.frame(matrix( 
        vector(), 0, length(selected_samples)+2, dimnames=list(c(),  c(newnames,"median","gene"))), 
        stringsAsFactors=F) 
      
      if (nrow(gexp)!=0) {
        ## sort by max samples values
        colnames(gexp) = newnames
        
        mmax = apply(gexp,1, function(x) if (all(is.na(x))) return(NA) else return(max(x, na.rm =T)))
        gexp$sample.max = mmax
        ## write.csv(gexp,"gexp_test.csv")
        ## cat(colnames(gexp), "\n")
        gexp = gexp[with(gexp,order(sample.max, decreasing = T)),]
        gexp$sample.max = NULL
        ##
        ## gdisease <- xsq_a[selxsq, disease_samples, drop = FALSE] 
        gdisease <- xsq_a[rownames(gexp), disease_samples, drop = FALSE] 
        ##
        ## write.csv(gdisease,"gdisease_test.csv")
  #      vmin = apply(gdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(min(x, na.rm =T))  )
  #      vmax = apply(gdisease, 1, function(x) if (all(is.na(x))) return(NA) else return(max(x, na.rm =T))  )
        vmedian = apply(gdisease, 1,function(x) if (all(is.na(x))) return(NA) else return(median(x, na.rm =T)) )
        gexp$median.disease = vmedian ##; gexp$min.disease = as.numeric(vmin) ; gexp$max.disease = as.numeric(vmax)
        gexp = round(gexp,3)
        gexp$gene = rownames(gexp)
        # gexp$type = "xsq"
        df = rbind(df, gexp)
      }
      
      cat("df dim 1: ", dim(df), " \n")
      
      
      shiny::validate(need(nrow(df)>0,"No data found for this patient"))
      nb = ncol(df)
      ## if (nb >= 3 ) 
      df = df[,c(nb, 1:(nb-1))]
      
      ####   df = df[order(df[,1]),] ## sort by gene name
      
      ## filter for gene value above median
      
      if (length(selected_samples) == 1 )  { 
        df = df[which(df[,newnames] >=  6),] 
        ## df = df[which(df[,selected_samples] >= df[,"median.disease"] & df[,selected_samples] >= 1),] 
      }
      else  {
         
        # genes_median <- rowSums(df[, selected_samples] >= df[,"median.disease"] & df[,selected_samples] >= 1 , na.rm=T) > 0
        # df <- df[genes_median, ]
        
        genes_sel <- rowSums(df[,newnames] >= 6 , na.rm=T) > 0
        df <- df[genes_sel, ]
        
      }
      ## df = df[which(df[,selected_samples] >= df[,"median.disease"]),] ## if only one sample
      
      cat("df dim 2 after filter: ", dim(df), " \n")
      # View(df)
      
      shiny::validate(need(nrow(df)>0,"No gene has value greater or equal to 6"))
      
      ## annotate with drug information
          resu = merge(adcgenes,df, by.x=1, by.y=1)
      # gann <- adcgenes[df$gene,]
      # df = cbind(gann,df)
      cat("Final dim df: ", dim(df), " dim resu: ", dim(resu),"\n")
      resu = resu[order(resu[,4], resu[,1]),]
      return(resu)
      
    }
    
    
    
  })
  
  output$trt_data <- DT::renderDataTable({
    datatable(
      my_trt(), rownames = F,
      filter = "top" , style='bootstrap' ,
      extensions = 'Buttons',
      options = list(
          columnDefs = list(
            list(visible = FALSE, targets = c(0, 2, 3)) # 
          ),
        ## dom = '<"top" pifB>',  # Displays the column visibility button dom = 'Bfrtip'
        dom='lipBt',
        scrollX = TRUE,
              buttons = list('colvis'),  # Adds the column visibility button
        # pageLength = input$pageLengthInputMethylation
        pageLength = 10, 
        lengthMenu = c(10, 25, 50, 100)
      )
    )
  })
  
  output$trt_data_v2 <- DT::renderDataTable({
    
    #### new custom head **** -------------------
    myframe = my_trt_v2()
    
    sketch = htmltools::withTags(table(
      class = 'display',
      thead(
        # Define the grouping of your df
        tr(
        #  th(rowspan = 1, ''),
          th(colspan = 2, 'OverExpressed Genes'),
          th(colspan = 4 , 'Candidate drugs'),
          th(colspan = ncol(myframe)-6, 'Patient Sample Expression')
        ),
        # Repeat column names 8 times
        tr(
          lapply(colnames(myframe), th)
        )
      )
    ))
    # Using JS for adding CSS, i.e., coloring your heading
    # Get the corresponding table header (th) from a table cell (td) and apply color to it
    headjs <- "function(thead) {
	  $(thead).closest('thead').find('th').eq(0).css('background-color', '#FFFFFF');
  $(thead).closest('thead').find('th').eq(1).css('background-color', '#8EA9DB');
   $(thead).closest('thead').find('th').eq(1).css('text-align', 'center');
   $(thead).closest('thead').find('th').eq(2).css('background-color', '#D9E1F2');
    $(thead).closest('thead').find('th').eq(2).css('text-align', 'center');
     }"
    #### end custom head  ****
    
    datatable(
     #  my_trt_v2(), 
      myframe,
      container = sketch, ## new head ****
      rownames = F,
      filter = "top" , style='bootstrap' ,
      extensions = 'Buttons',
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = c( 1,4)) # 
        ),
        ## dom = '<"top" pifB>',  # Displays the column visibility button dom = 'Bfrtip'
        dom='lipBt',
        scrollX = TRUE,
        headerCallback = DT::JS(headjs), ## new head ****
        buttons = list('colvis'),  # Adds the column visibility button
        # pageLength = input$pageLengthInputMethylation
        pageLength = 10, 
        lengthMenu = c(10, 25, 50, 100)
      )
    )
  })
## end ADC treatment ***  -----------------------------------
  
  
  
  #SLIDERS
  
 
 output$gene_slider <- renderUI({
   if (input$slider_type == "Top and Bottom Genes") {
     sliderInput(
       session$ns("num_genes_displayed"),
       "Top and Bottom Genes (based on mean or stdev in case of 2 or more than 2 samples)",
       min = 10,
       max = 2000,
       value = 20,  ## was 10
       step = 10
     )
   } else {
     geneRange <- my_gene_expression()
     gene_data <- geneRange$my_data
     gene_samples <- geneRange$sample_list
     
     min_expr <- round(min(gene_data[, gene_samples]), 3)
     max_expr <- round(max(gene_data[, gene_samples]), 3)
     
     sliderInput(
       session$ns("gene_expr_range"),
       "Gene Expression Range (For At Least 1 Sample)",
       min = min_expr, 
       max = max_expr,
       value = c(min_expr,  max_expr),
       step = 0.5
     )
   }
 })
 
 output$mutation_slider <- renderUI({
   if (input$slider_type_mutation == "Top and Bottom Mutations") {
     sliderInput(
       session$ns("num_mut_displayed"),
       "Top and Bottom Mutations",
       min = 10,
       max = 2000,
       value = 20,
       step = 10
     )
   } else {
     mutRange <- my_mutation()
     mut_data <- mutRange$my_data_mutation
     mut_samples <- mutRange$sample_list_mutation
     
     min_mut <- round(min(mut_data[, mut_samples]), 3)
     max_mut <- round(max(mut_data[, mut_samples]), 3)
     
     sliderInput(
       session$ns("mut_range"),
       "Mutation Range For At Least 1 Sample",
       min = min_mut, 
       max = max_mut,
       value = c(min_mut,  max_mut),
       step = 0.1
     )
   }
 })

 output$methylation_slider <- renderUI({
   if (input$slider_type_methylation == "Top and Bottom Methylation Scores") {
     sliderInput(
       session$ns("num_methylation_displayed"),
       "Top and Bottom Methylation Scores",
       min = 10,
       max = 2000,
       value = 20,
       step = 10
     )
   } else {
     metRange <- my_methylation()
     met_data <- metRange$my_data_methylation
     met_samples <- metRange$sample_list_methylation

     min_met <- round(min(met_data[, met_samples]), 3)
     max_met <- round(max(met_data[, met_samples]), 3)

     sliderInput(
       session$ns("methylation_range"),
       "Methylation Range For At Least 1 Sample",
       min = min_met,
       max = max_met,
       value = c(min_met,  max_met),
       step = 0.05
     )
   }
 })
 
 output$copy_slider <- renderUI({
   if (input$slider_type_copy == "Top and Bottom copy Scores") {
     sliderInput(
       session$ns("num_copy_displayed"),
       "Top and Bottom Copy Scores",
       min = 10,
       max = 2000,
       value = 20,
       step = 10
     )
   } else {
     copRange <- my_copy()
     cop_data <- copRange$my_data_copy
     cop_samples <- copRange$sample_list_copy
     
     min_cop <- round(min(cop_data[, cop_samples]), 3)
     max_cop <- round(max(cop_data[, cop_samples]), 3)
     
     sliderInput(
       session$ns("copy_range"),
       "Copy Range For At Least 1 Sample",
       min = min_cop,
       max = max_cop,
       value = c(min_cop,  max_cop),
       step = 0.05
     )
   }
 })
 

}

  
  
  

  
  
