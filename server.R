library(shiny)
#library(d3heatmap)
library(rcellminer)
library(rcellminerElasticNet)
library(geneSetPathwayAnalysis)
library(jsonlite)
library(stringr)
library(glmnet)
library(ggplot2)
library(plotly)
##library(xlsx)
library(shinycssloaders)
library(gplots)
library(heatmaply)
library(memoise)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(maftools)

## library(survminer) ## for forest plot


## library(shinyHeatmaply)
library(dplyr)
library(shinyalert)
library(shinyjs)
#
library(DT)
library(reshape2)
library(tumorcomparer)

#library(tooltipsterR)

if (!require(rcellminerUtilsCDB)){
	warning("rcellminerUtilsCDB package must be installed for full cross-database functionality.")
}


#--------------------------------------------------------------------------------------------------
# LOAD CONFIGURATION AND REQUIRED DATA SOURCE PACKAGES.
#--------------------------------------------------------------------------------------------------
config <- jsonlite::fromJSON("config.json")
appConfig <- jsonlite::fromJSON("appConfig.json")
metaConfig <- jsonlite::fromJSON("configMeta.json")

oncolor <- read.delim("oncotree1_colors.txt",row.names = 1,stringsAsFactors = F)
rownames(oncolor)=toupper(rownames(oncolor))

mygeneset = geneSetPathwayAnalysis::geneSets
list.genesets = names(mygeneset)
ll = unlist(lapply(mygeneset, function(x) length(mygeneset[x])))
list.genesets = paste0(list.genesets," (",ll,")")
names(mygeneset) = list.genesets
  
 
acc_list = read.delim("acc_genes_4997.txt",stringsAsFactors = F)[,1]
accgeneset = list("ACC variable genes (4997)" = acc_list)
mygeneset = append(accgeneset, mygeneset)

## source("modal1.R")
source("appUtils.R")
source("dataLoadingFunctions.R")

#if (!is.null(appConfig$appName)){
#	appTitle <- appConfig$appName
#} else{
#	appTitle <- "CellMiner"
#}
cacheDir <- appConfig$cacheDir
downloadDir <- appConfig$downloadDir

# Construct named character vector mapping displayed data source names to
# internally used source identifiers.
dataSourceChoices <- setNames(names(config),
															vapply(config, function(x) { x[["displayName"]] }, 
																		 character(1)))

metaChoices <- setNames(names(metaConfig),
												vapply(metaConfig, function(x) { x[["displayName"]] }, 
															 character(1)))
if(!file.exists("srcContent.rds")) {
	
 for (configSrcId in names(config)){
	srcName <- config[[configSrcId]][["displayName"]]
	srcPackages <- names(config[[configSrcId]][["packages"]])
	for (pkgName in srcPackages){
		 if (!require(pkgName, character.only = TRUE)){
			  dataSourceChoices[srcName] <- NA
		  	break
	  	}
	 }
  }

 if (any(is.na(dataSourceChoices))){
	stop("Check configuration file: one or more required data source packages must be installed.")
 } 
  ## new staff -------------------------------------------------------
  cat("creating RDS content file\n")
  ## srcContent <- lapply(config, loadSourceContent)
  srcContent <- lapply(config, loadSourceContentFiltered,onco1="Lung",onco2="Small Cell Lung Cancer (SCLC)")
  isLoadedSrc <- vapply(srcContent, function(x) { !is.null(x) }, logical(1))
  if (any(!isLoadedSrc)){
    srcContent <- srcContent[isLoadedSrc]
  }
  
  
  saveRDS(srcContent, "srcContent.rds", compress = FALSE)
  cat("RDS content file created ! \n")
  ## end new staff ---------------------------------------------------
  
	
} else {
	srcContent <- readRDS("srcContent.rds")
}
#--------------------------------------------------------------------------------------------------

#if("rCharts" %in% installed.packages()) {
#	options(RCHART_LIB='highcharts')	
#	library(rCharts)
#	hasRCharts <- TRUE
#} else {
#	hasRCharts <- FALSE
#}

colorSet <- loadNciColorSet(returnDf=TRUE)

###--------

options("DT.TOJSON_ARGS" = list(na = "string")) ## try dev version of DT

#--------------------------------------------------------------------------------------------------
sysinfo <- Sys.info()
if (sysinfo["nodename"]=="discovery.nci.nih.gov" | sysinfo["nodename"]=="ncias-d2059-v.nci.nih.gov") {
##db <- cache_filesystem("/srv/shiny-server/cellminercdb_internal/.rcache")
db <- cache_filesystem("/data/patientMiner/.rcache")
downloadpath <- "/data/patientMiner-downloads"
 } else {
  if (sysinfo["nodename"]=="discover.nci.nih.gov" | sysinfo["nodename"]=="ncias-p2122-v.nci.nih.gov")  {
    db <- cache_filesystem("/data/patientMiner/.rcache") 
    downloadpath <- "/data/patientMiner-downloads"
    }
   else {
     # db <- cache_filesystem("/Users/elloumif/patientCache/.rcache")
     # downloadpath <- "/Users/elloumif/patientMiner-downloads"
     db <- cache_filesystem(cacheDir)
     downloadpath <- downloadDir
   }
}
patternComparison <- memoise(rcellminer::patternComparison, cache = db)
# patternComparison <- memoise(rcellminer::patternComparison) # cache = cache_memory()
getMolDataType <- memoise(rcellminer::getMolDataType, cache = db)
removeMolDataType <- memoise(rcellminer::removeMolDataType, cache = db)

# sink("sessioninfo.txt")
# print(sessionInfo())
# sink()

# library(bigrquery)
# library(httpuv)
aproject <-"isb-cgc-fathi"

## enableBookmarking("url") ### bookmarking >> global.R

shinyServer(function(input, output, session) {
  ######### bookmarking
  # Automatically bookmark every time an input changes
  observe({
    # reactiveValuesToList(input)
    
    toExclude = setdiff(names(input), c("xId","xPrefix","yId","yPrefix"))
    setBookmarkExclude(toExclude)
    session$doBookmark()
  })
  # Update the query string
  onBookmarked(updateQueryString)
  ####### end bookmarking
  ##########-------------------#################
  distPlot <-eventReactive(input$subtcga,{
    prefixChoices <- srcContent[[input$cmpSource]][["featurePrefixes"]]
    pf="xsq"
    if (is.na(match("xsq",prefixChoices))) 
    {
      if (is.na(match("exp",prefixChoices)))  pf=NA else pf="exp"
    }
    shiny::validate(need(!is.na(pf), 
                         "There is no gene expression for selected cell line set."))
    wdata=srcContent[[input$cmpSource]][["molPharmData"]][[pf]]
    rownames(wdata)=substr(rownames(wdata),4,nchar(rownames(wdata)))
    ## search for gene in cell line set - for now no synonym
    shiny::validate(need(toupper(input$vgene)  %in% rownames(wdata), 
                         "gene not found in selected cell line set"))
    resu = cbind(genexp=wdata[toupper(input$vgene),],cohort=srcContent[[input$cmpSource]][["sampleData"]][["OncoTree1"]])
    resu = data.frame(resu,stringsAsFactors = F)
    resu[,1] = as.numeric(resu[,1])
    resu[,2] = as.factor(resu[,2])
    # new
    if (input$zsid) {
    resu = resu[order(resu[,2]),]
    resu[,1]=unlist(tapply(resu[,1],resu[,2],scale))
    }
    #
    print(dim(resu))
    ## add TCGA
    dat <- geneExpTcga(toupper(input$vgene), aproject)
    shiny::validate(need(!is.null(dat), 
                         "gene not found in TCGA"))
    ## multiplatform >> batch effect ??
    if (length(levels(factor(dat$platform))>1)) showNotification("TCGA platform is not unique, possible batch effect.", duration=10, type="warning")
    ## new
    print(dim(dat))
    vind = which(as.numeric(dat$ttype)<10)
    tcga = dat[vind,c("normalized_count","project_short_name")]
    colnames(tcga)=c("genexp","cohort")
    #
    if (input$zsid) {
    tcga  = tcga[order(tcga$cohort),]
    tcga$genexp = unlist(tapply(tcga$genexp,tcga$cohort,scale))
    }
    ## select only tumor samples
    
    resu= rbind(resu,tcga)
    print(dim(resu))
    ## 
    p<-ggplot(resu, aes(x=cohort, y=genexp)) +
     
     ##geom_boxplot(aes(fill = cohort)) + theme(plot.title = element_text(size=14, face="bold"), legend.text = element_text(size=14),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank() ) + ggtitle(paste0("Distribution of ",toupper(input$vgene)," gene expression (z-score)"))
  
      geom_boxplot(aes(fill = cohort)) + theme(plot.title = element_text(size=14, face="bold"), legend.text = element_text(size=14),axis.title.x = element_blank(), axis.text.x = element_text(angle = 45) ) + ggtitle(paste0("Distribution of ",toupper(input$vgene)," gene expression"))
    
     g <- ggplotly(p,plotWidth=1000,plotHeight=1400)

        ### g <- ggplotly(p,plotWidth=1400,plotHeight=2000)
        ### g <- layout(g, margin=list(t = 75, b = 0))
    g <- layout(g, margin=list(t = 75), legend = list(font = list(size = 18)))
    g1 <- config(p = g, cloud=FALSE, displaylogo=FALSE, displayModeBar=TRUE,
                 modeBarButtonsToRemove=c("select2d", "sendDataToCloud", "pan2d", "resetScale2d",
                                          "hoverClosestCartesian", "hoverCompareCartesian",
                                          "lasso2d", "zoomIn2d", "zoomOut2d"))
    g1
      })
  
  ## output$tcgaPlot <- renderPlot({distPlot()}, width=1100, height=1100)
  getDistData <- reactive({
    wdata=srcContent[[input$dataDist]][["molPharmData"]][["xsq"]]
    rownames(wdata)=substr(rownames(wdata),4,nchar(rownames(wdata)))
    ## search for gene in cell line set - for now no synonym
    shiny::validate(need(toupper(input$geneDist)  %in% rownames(wdata), 
                         "gene not found !!!"))
    # resu = cbind(genexp=wdata[toupper(input$geneDist),],disease_source=srcContent[[input$dataDist]][["sampleData"]][["dataSource"]],sampletype = , priortreatment = )
    resu = cbind(genexp=wdata[toupper(input$geneDist),],disease_source=srcContent[[input$dataDist]][["sampleData"]][["dataSource"]])
    ## write.csv(resu, "resu0.csv")
    
    sinfo = cbind(sampleType=srcContent[[input$dataDist]][["sampleData"]][["sampleType"]], priorTreatment=srcContent[[input$dataDist]][["sampleData"]][["priorTreatment"]])
    resu = data.frame(resu,stringsAsFactors = F)
    resu[,1] = as.numeric(resu[,1])
    resu[,2] = as.factor(resu[,2])
    
    if (input$sampletype2) {
      ii= which(sinfo[,1]=="SampleType:Primary")
      shiny::validate(need(length(ii)>0, 
                           "No primary sample found !"))
      resu = resu[ii,]
      sinfo = sinfo[ii,]
    }
    # cat("1", dim(geneexp), " ", dim(clin), "\n")
    if (input$therapy2) {
      jj= which(sinfo[,2]=="PriorTreatment:none")
      shiny::validate(need(length(jj)>0, 
                           "No sample without prior treatment found !"))
      resu = resu[jj,]
      sinfo = sinfo[jj,]
    }
    kk = which(is.na(resu[,1]))
    if (length(kk)>0) {
      resu = resu[-kk,]
      sinfo = sinfo[-kk,]
    }
    return(list(resu=resu,sinfo=sinfo))
  })
  
  output$genePlot <- renderPlotly({
    
    # wdata=srcContent[[input$dataDist]][["molPharmData"]][["xsq"]]
    # rownames(wdata)=substr(rownames(wdata),4,nchar(rownames(wdata)))
    # ## search for gene in cell line set - for now no synonym
    # shiny::validate(need(toupper(input$geneDist)  %in% rownames(wdata), 
    #                      "gene not found !!!"))
    # # resu = cbind(genexp=wdata[toupper(input$geneDist),],disease_source=srcContent[[input$dataDist]][["sampleData"]][["dataSource"]],sampletype = , priortreatment = )
    # resu = cbind(genexp=wdata[toupper(input$geneDist),],disease_source=srcContent[[input$dataDist]][["sampleData"]][["dataSource"]])
    # ## write.csv(resu, "resu0.csv")
    # 
    # sinfo = cbind(sampleType=srcContent[[input$dataDist]][["sampleData"]][["sampleType"]], priorTreatment=srcContent[[input$dataDist]][["sampleData"]][["priorTreatment"]])
    # resu = data.frame(resu,stringsAsFactors = F)
    # resu[,1] = as.numeric(resu[,1])
    # resu[,2] = as.factor(resu[,2])
    # 
    # if (input$sampletype2) {
    #   ii= which(sinfo[,1]=="SampleType:Primary")
    #   shiny::validate(need(length(ii)>0, 
    #                        "No primary sample found !"))
    #   resu = resu[ii,]
    #   sinfo = sinfo[ii,]
    # }
    # # cat("1", dim(geneexp), " ", dim(clin), "\n")
    # if (input$therapy2) {
    #   jj= which(sinfo[,2]=="PriorTreatment:none")
    #   shiny::validate(need(length(jj)>0, 
    #                        "No sample without prior treatment found !"))
    #   resu = resu[jj,]
    #   sinfo = sinfo[jj,]
    #  }
    # kk = which(is.na(resu[,1]))
    # if (length(kk)>0) {
    #   resu = resu[-kk,]
    #   sinfo = sinfo[-kk,]
    # }
    ### resu = na.omit(resu)
    resu = getDistData()$resu
    sinfo = getDistData()$sinfo
    print(dim(resu))
    origin_group = resu$disease_source
    ## write.csv(cbind(resu,sinfo), "resu1.csv")
    group_ordered <- with(resu,                       # Order boxes by median
                          reorder(disease_source,
                                  genexp,
                                  function(x) median(x,na.rm=T)))
    resu$disease_source = factor(resu$disease_source,levels = levels(group_ordered))
    ##       ## 
    ## mycols = c("red","green","blue","brown")
    ## mycols = c("green","blue","brown","red")
    mycols = c("green","brown","red","blue")
    
    p<-ggplot(resu, aes(x=disease_source, y=genexp)) +
      
      ##geom_boxplot(aes(fill = source)) + theme(plot.title = element_text(size=14, face="bold"), legend.text = element_text(size=14),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank() ) + ggtitle(paste0("Distribution of ",toupper(input$vgene)," gene expression (z-score)"))
    
      geom_boxplot(aes(fill = origin_group)) +   
    # geom_boxplot(aes(fill = disease_source)) + 
      theme_bw() +  theme(plot.title = element_text(size=18, face="bold"), legend.text = element_text(size=14),legend.title=element_text(size = 14),axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, size =14) , axis.title.y = element_text(angle = 0, size =14), axis.text.y = element_text(angle = 0, size =14)) + ggtitle(paste0("Distribution of ",toupper(input$geneDist)," gene expression for ",nrow(resu)," patient samples")) +
    stat_summary(fun.data = function(x){
      return(data.frame(y = max(x)*1.05, label = length(x)))
    }, geom="text", size=5.5) + scale_fill_manual(values = mycols)
    
    # g <- ggplotly(p,plotWidth=1000,plotHeight=1400)
    g <- ggplotly(p)
    # g <- ggplotly(p,plotWidth=1400,plotHeight=2000)
    ### g <- layout(g, margin=list(t = 75, b = 0))
    g <- layout(g, margin=list(t = 75), legend = list(font = list(size = 18)) , hoverlabel = list(font=list(size=20))  )
    g1 <- config(p = g, cloud=FALSE, displaylogo=FALSE, displayModeBar=TRUE,
                 modeBarButtonsToRemove=c("select2d", "sendDataToCloud", "pan2d", "resetScale2d",
                                          "hoverClosestCartesian", "hoverCompareCartesian",
                                          "lasso2d", "zoomIn2d", "zoomOut2d"))
    g1
    
    })
  ##--------------
  output$downloadPlotData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste0("Expression_Distribution_",input$geneDist,"_primary_",input$sampletype2,"_NOpriortreatment_",input$therapy2,".csv")
    },
    
    content = function(file) {
     
      resu = getDistData()$resu
      sinfo = getDistData()$sinfo
      write.csv(cbind(resu,sinfo), file)
      # write.table(dataMatrix, file, sep = "\t", col.names = NA,quote=F)  
      
    }
  )
  
  ## survival plot  
  output$SurvplotNew <- renderPlot({
    cohort   <- as.character(input$onco3ch)
    ## cat(cohort ,"...\n")
    geneexp  <- srcContent[[input$survSource]][["molPharmData"]][[input$optsurv]] 
    clin     <- srcContent[[input$survSource]][["sampleData"]]
    ## cat(input$maxSurvPredictorsLasso,"  maxsurv ..\n")
    # drawPlotExp3(cohort, input$varname, geneexp, clin, input$sampletype,input$therapy,input$optsurv)
    ## drawPlotExp4(cohort, input$varname, geneexp, clin, input$sampletype,input$therapy,input$optsurv, input$Lasso,input$inputGSsurv ,input$maxSurvPredictorsLasso)
    drawPlotExpNew(cohort, input$varname, geneexp, clin, input$sampletype,input$therapy,input$optsurv, input$Lasso,input$inputGSsurv ,input$maxSurvPredictorsLasso, input$optgender)
  },  width=1400, height=800)
  
  ### new survical reactive functions -----------------------------------------
  
  findInfoExpReact <- reactive({
    #function(varName, cohort, geneexp, clin,sampletype,priortrt,optsurv,chkLasso,gsets,nbpred,gender) {
    # find gene and survival information
    # sampleinfo is the oncotree3 >> dataSource
    varName = input$varname
    cohort = as.character(input$onco3ch)
    geneexp = srcContent[[input$survSource]][["molPharmData"]][[input$optsurv]] 
    clin = srcContent[[input$survSource]][["sampleData"]]
    sampletype = input$sampletype
    priortrt = input$therapy
    optsurv = input$optsurv
    chkLasso = input$Lasso
    gsets = input$inputGSsurv
    nbpred = input$maxSurvPredictorsLasso
    gender = input$optgender
    
    if (sampletype) { 
      geneexp = geneexp[, which(clin$sampleType=="SampleType:Primary")] 
      clin = clin[which(clin$sampleType=="SampleType:Primary"),]
      # sampleinfo = sampleinfo[which(clin$sampleType=="SampleType: Primary")]
    }
    # cat("1", dim(geneexp), " ", dim(clin), "\n")
    if (priortrt) {
      geneexp = geneexp[, which(clin$priorTreatment=="PriorTreatment:none")]
      clin    = clin[which(clin$priorTreatment=="PriorTreatment:none"),]
      # sampleinfo = sampleinfo[which(clin$priorTreatemnt=="PriorTreatment: none")]
    }
    # cat("2", dim(geneexp), " ", dim(clin), "\n")
    
    if (gender == "Female") {
      geneexp = geneexp[, which(clin$gender=="Gender:Female")]
      clin    = clin[which(clin$gender=="Gender:Female"),]
    } else
    {
      if (gender == "Male") {
        geneexp = geneexp[, which(clin$gender=="Gender:Male")]
        clin    = clin[which(clin$gender=="Gender:Male"),]
      }
    }
    # cat("3", dim(geneexp), " ", dim(clin), "\n")
    # Focus on non redundant samples
    geneexp = geneexp[, which(clin$survivalStatus==TRUE)]
    clin    = clin[which(clin$survivalStatus==TRUE),]
    # cat("3", dim(geneexp), " ", dim(clin), "\n") 
    #
    dat =  data.frame()
    ## new stuff
    rownames(geneexp) = toupper(rownames(geneexp))
    gl = trimws(toupper(varName))
    if (optsurv=="xsq") gl = unlist(strsplit(gl,"\\s+"))
    
    if (chkLasso & optsurv=="xsq" ) {
      shiny::validate(need(nbpred>=1 & nbpred <=10, "Error: The number of genes should be between 1 and 10"))
      shiny::validate(need(length(gsets) > 0, "Please select one or more gene sets."))
      ## cat(gsets, "...gsets..\n")
      ## gl = union(gl, geneSetPathwayAnalysis::geneSets[[gsets]])
      # gl = union(gl, unique(c(geneSetPathwayAnalysis::geneSets[gsets], recursive = TRUE))) ## without ACC
      gl = union(gl, unique(c(mygeneset[gsets], recursive = TRUE))) ## new ACC
      comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp)) ## adding prefix
      ## cat(length(gl), ".. \n")
      ## ge <- geneexp[comgenes, which(clin$OncoTree3==cohort), drop=F]
      # ge <- geneexp[comgenes, which(clin$OncoTree3 %in% cohort), drop=F]
      ge <- geneexp[comgenes, which(clin$dataSource %in% cohort), drop=F]
      rownames(ge) =  paste0(tolower(substr(rownames(ge),1,3)), substr(rownames(ge),4,nchar(rownames(ge))))
      ## cat(dim(ge), "..\n")
      
      # check if samples have NA for all genes so remove them for Lasso --------------
      vna = apply(ge,2, function(x) length(which(is.na(x))))
      kna = which(vna==nrow(ge))
      if (length(kna)!=0) ge = ge[,-kna]
      # ------------------------------------------------------------------------------
      shiny::validate(need(ncol(ge)>=10, "Error: Min number of samples is 10"))
      
      months_to_last_known_alive <- as.numeric(as.character(clin[colnames(ge),"osMonths"]))
      vital_status <- clin[colnames(ge),"vitalStatus"]
      vital_status <- ifelse(vital_status == "alive", 0, 1)
      ge = data.frame(ge)
      
      # length(vital_status)
      ## cat(dim(ge), " dim ge..................\n") 
      # ww = cbind(t(ge), time = months_to_last_known_alive, status = vital_status)
      # write.csv(ww,"tempcox.csv")
      # set.seed(1)
      ## cat(months_to_last_known_alive, " days alive .....\n")
      
      # checking NA in status or time
      jna = which(is.na(months_to_last_known_alive) | is.na(vital_status))
      if (length(jna)!=0) {
        ge = ge[,-jna]
        months_to_last_known_alive = months_to_last_known_alive[-jna]
        vital_status = vital_status[-jna]
        shiny::validate(need(ncol(ge)>=10, "Error: Min number of samples is 10"))
      }
      
      shiny::validate(need(sum(vital_status)!=0, 
                           "Can not compute HR, no events found."))
      
      y = cbind(time = months_to_last_known_alive, status = vital_status)
      #test
      # write.csv(y,"y.csv")
      #
      set.seed(1) ## right one 
      fit1_cv = cv.glmnet(t(ge), y, family = "cox", alpha = 1) # default alpha=1 >> Lasso
      
      # summary(fit1_cv)
      # plot(fit1_cv)
      # title("TCGA ACC Lasso, Cox Family", line = 2.5)
      lassoPredictorWts0 <- coef(fit1_cv, s = "lambda.min")[,1] # deviance
      lassoPredictorWts0 <- lassoPredictorWts0[lassoPredictorWts0 != 0]
      lassoPredictorWts0 <- lassoPredictorWts0[order(abs(lassoPredictorWts0), decreasing = TRUE)]
      
      ## nb max of genes
      if (length(lassoPredictorWts0) > nbpred) lassoPredictorWts0 = lassoPredictorWts0[1:nbpred]
      shiny::validate(need(length(lassoPredictorWts0)>0, "No Lasso predictors found"))
      ## cat(lassoPredictorWts0, "selected beta ......")
      # based on lambda.min
      res = data.frame(beta = lassoPredictorWts0, HR = exp(lassoPredictorWts0), stringsAsFactors = F)
      res$p.value = NA
      rownames(res) = names(lassoPredictorWts0)
      beta = lassoPredictorWts0
      selgenes = names(lassoPredictorWts0)
      ## cat(selgenes, "selected genes .....\n")
      ## dat = data.frame(t(ge[selgenes,]))
      dat = data.frame(t(ge[selgenes,]))
      ## cat(dim(dat), "  dat  ..\n")
      rscore = t(ge[selgenes,]) %*% beta ## risk score
      
      dat$variable = rowMeans(t(ge[selgenes,]))
      ## dat$variable = rscore
      dat$rscore = rscore
      
      ## cat(rscore, " ..rscore  ..\n")
      ## dat$variable = rscore
      dat$months_to_last_known_alive = months_to_last_known_alive
      dat$vital_status = vital_status
     ## write.csv(dat, "dat_test_v3_LASSO.csv")
      ## cat(dim(dat), "  ..dat\n")
      ## cat(dim(res), "  ..cox multi\n")
      cesp0 = coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)
      beta0 = unlist(summary(cesp0)[7])[1]
      HR0 = unlist(summary(cesp0)[7])[2]
      p.value0 = unlist(summary(cesp0)[7])[5]
      # resavg = data.frame(beta = beta0, HR = HR0,  p.value = p.value0, stringsAsFactors = F)
      # rownames(resavg) = "Average_gene(s)"
      
      cesp1 = coxph(Surv(months_to_last_known_alive, vital_status)~rscore,data=dat)
      beta1 = unlist(summary(cesp1)[7])[1]
      HR1 = unlist(summary(cesp1)[7])[2]
      p.value1 = unlist(summary(cesp1)[7])[5]
      
      resavg = data.frame(beta = c(beta1,beta0), HR = c(HR1,HR0),  p.value = c(p.value1,p.value0), stringsAsFactors = F)
      rownames(resavg) = c("Risk_Score (weighted average genes using beta)", "Mean_Score (arithmetic average genes)")
      
      return(list(data=dat,comgenes=selgenes,coxmulti = res, coxaverage = resavg, coxmodel = NULL))
    }
    
    else
    {  
      
      # mutation1
      if (optsurv=="mut") geneexp[geneexp>0] = 1
      
      nb=length(gl)
      shiny::validate(need(nb<11, 
                           "Error: Max number of genes is 10"))
      comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp))
      
      shiny::validate(need(length(comgenes)> 0, 
                           "Error: Gene(s) not found !"))
      
      if (length(comgenes)> 0) 
        ## 
        ## if(!is.na(match(toupper(varName), rownames(geneexp))))
      { 
        ##   ge <- geneexp[toupper(varName), which(clin$OncoTree3==cohort)]
        ## new
        ## ge <- geneexp[comgenes, which(clin$OncoTree3==cohort), drop=F]
        # ge <- geneexp[comgenes, which(clin$OncoTree3 %in% cohort), drop=F]
        ge <- geneexp[comgenes, which(clin$dataSource %in% cohort), drop=F]
        # dat = data.frame(patient = clin[colnames(ge),"Name"], stringsAsFactors = F)
        
        dat = data.frame(t(ge))
        colnames(dat) = paste0(tolower(substr(colnames(dat),1,3)), substr(colnames(dat),4,nchar(colnames(dat))))
        comgenes = colnames(dat)
        ## dat$normalized_count <- colMeans(ge)
        # then we need to do a little data-cleaning to get ready for our survival model
        dat$months_to_last_known_alive <- as.numeric(as.character(clin[colnames(ge),"osMonths"]))
        dat$vital_status <- clin[colnames(ge),"vitalStatus"]
        dat$vital_status <- ifelse(dat$vital_status == "alive", 0, 1)
        
        ## cat("4 dat", dim(dat), class(dat), "\n") 
      }
      # testing
      # print(nrow(dat))
      # write.csv(dat,"dat_test.csv")
      # 
      shiny::validate(need(nrow(na.omit(dat))>1, 
                           "Error: Many missing values, found 1 sample"))
      
      shiny::validate(need(length(na.omit(dat[,-match(c("months_to_last_known_alive","vital_status"),colnames(dat))]))>0, 
                           "all values are NA"))
      #
      #print(dat$normalized_count)
      # first compute cox
      
      nbc = length(comgenes)
      
      if (nbc==1) {
        dat$variable = dat[,comgenes]
        cesp = coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)
        beta = unlist(summary(cesp)[7])[1]
        HR = unlist(summary(cesp)[7])[2]
        p.value = unlist(summary(cesp)[7])[5]
        res = data.frame(beta = beta, HR = HR,  p.value = p.value, stringsAsFactors = F)
        rownames(res) = comgenes
        resavg = res
        rownames(resavg) = "Average_genes"
        ## to use formatC(c(1.543,6,8900), format = "e", digits = 2)
        
      } else {
        cesp=coxph(Surv(months_to_last_known_alive, vital_status)~.,data=dat)
        # print(summary(cesp))
        beta = unlist(summary(cesp)[7])[1:nbc]
        
        shiny::validate(need(!all(is.na(beta)), 
                             "Can not compute HR, check number of events."))
        
        HR = unlist(summary(cesp)[7])[(1+nbc):(nbc+nbc)]
        p.value = unlist(summary(cesp)[7])[(1+4*nbc):(nbc+4*nbc)]
        res = data.frame(beta = beta, HR = HR,  p.value = p.value, stringsAsFactors = F)
        rownames(res) = comgenes
        rscore = as.matrix(dat[,comgenes]) %*% beta ## risk score
        ## print(rscore)
        dat$variable = rowMeans(as.matrix(dat[,comgenes]))
        ## dat$variable = rscore
        dat$rscore = rscore
        
        write.csv(dat, "dat_test_v3.csv")
        
        cesp0 = coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)
        beta0 = unlist(summary(cesp0)[7])[1]
        HR0 = unlist(summary(cesp0)[7])[2]
        p.value0 = unlist(summary(cesp0)[7])[5]
        
        cesp1 = coxph(Surv(months_to_last_known_alive, vital_status)~rscore,data=dat)
        beta1 = unlist(summary(cesp1)[7])[1]
        HR1 = unlist(summary(cesp1)[7])[2]
        p.value1 = unlist(summary(cesp1)[7])[5]
         
        ## print(summary(cesp1))
         
        resavg = data.frame(beta = c(beta1,beta0), HR = c(HR1,HR0),  p.value = c(p.value1,p.value0), stringsAsFactors = F)
        rownames(resavg) = c("Risk_Score (weighted average genes using beta)", "Mean_Score (arithmetic average genes)")
      }
      return(list(data=dat,comgenes=comgenes,coxmulti = res,coxaverage = resavg, coxmodel = cesp))
      # new for Lasso
    }
    
    #
    
  }
  )
  
  
  ### end survival reactive ------------------------------------------
  
  output$SurvplotReact <- renderPlot({
    
    resu <- findInfoExpReact()
    rcox <- resu$coxmulti ## new to add table in the same plot using gridExtra !!!! to do
    dat = resu$data
    cohort = as.character(input$onco3ch)
    optsurv = input$optsurv
    chkLasso = input$Lasso
    gender = input$optgender
    optgroup = input$optgroup
    # covariates = resu$comgenes
    comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
    
     tlabchange=F
    response.M25 = dat$variable
    # mutation2
    if (!((comgenes %in% c("CGA.STAIN","HORMONE.PROD","SYP.POSITIVE", "INSM1","MALEGENDER"))[1]) & (optsurv!="mut") ) {
      if (optgroup == "25") {
      response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[2],0,NA)
      response.M25[which(dat$variable>quantile(dat$variable,na.rm=T)[4])]=1
      } else
      {
        if (optgroup == "50") 
        {
          response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[3],0,NA)
          response.M25[which(dat$variable>quantile(dat$variable,na.rm=T)[3])]=1
        } else
        {
          response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T, probs = seq(0, 1, 0.33) )[2],0,NA)
          response.M25[which(dat$variable>quantile(dat$variable,na.rm=T, probs = seq(0, 1, 0.33))[3])]=1
        }
        
      }
      tlabchange = T
    } 
    
    
    
    # nb1 = length(which(response.M25==0))
    # nb2 = length(which(response.M25==1))
    
    # new counts
    inb1 = which(response.M25==0)
    inb2 = which(response.M25==1)
    nb1 = length(which(!is.na(dat$months_to_last_known_alive[inb1])))
    nb2 = length(which(!is.na(dat$months_to_last_known_alive[inb2])))
    # 
    ## now using cox model
    c0=coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)
    nbobs= unlist(summary(c0)[4])

    issue = F ## new
    ## c=coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)
    c = tryCatch({
      coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)
    }, warning = function(w) {
      
      ## print(w)
      shinyalert::shinyalert(
        # title = "Warning in Cox 2 groups",
        # text = paste("Message: ",w),
        title = paste("Warning in Cox 2 groups:\n",w),
        text = "",
        type = "warning",
        imageHeight = 150,
        size ="m"
      ) 
      ##
      # mess <<- w
      # issue <<- T
    }, error = function(e) {
      shinyalert::shinyalert(
        # title = "Error in Cox 2 groups",
        # text = paste("Message: ",e),
        title = paste("Error in Cox 2 groups:\n",e),
        text = "",
        type = "error",
        imageHeight = 150,
        size ="m"
      ) 
      mess <<- e
      issue <<- T
    }
    )
    ## error case 
    ##
    if (class(c) != "coxph") c= coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)
    mfit.subtype1=survfit(c,newdata=data.frame(response.M25=c(0,1)))
    
    if (length(comgenes)<=10) genlab= paste(comgenes,collapse=",") else genlab = paste(paste(comgenes[1:10],collapse=","),paste(comgenes[11:length(comgenes)],collapse=","), sep="\n")
    labcohort = paste(cohort,collapse=",")
    ## improve cohort display
    ## integer div  % / % ,  mod %% >> rest 
    if (length(cohort)>10) {
      labcohort = ""
      nbl = length(cohort) %/% 10
      reste = length(cohort) %% 10
      for (k in seq(1,(nbl*10), by=10)) {
        labcohort = paste0(labcohort,paste(cohort[k:(k+9)],collapse=","),"\n")
      }
      if (reste!=0) labcohort = paste0(labcohort,paste(cohort[(k+10):(k+9+reste)],collapse=","))
    }
    
    ## end improvments
    
    if (optsurv=="xsq") {
      if (chkLasso) { stitle = " Average Expression of Lasso selected genes: " } else {
        if (length(comgenes)==1) stitle = " gene expression of "  else stitle = " Average Expression of entered genes: " }
    }
    else
    {
      if (optsurv=="mda") 
          stitle = " metadata using "  else 
        { if (optsurv=="mut")
            stitle = " mutation of " 
          else
           { if (optsurv=="sig")
                stitle = " signature using "  else 
                stitle = " NMF cluster using "
           }
        }
    }
    
    if (gender=="Both") mesgender = "All" else mesgender = gender
    
    ##plot(mfit.subtype1,main=paste("Survival prediction in ",cohort," based on \n",stitle, paste(comgenes,collapse=",")," with ",nbobs, " samples.",sep=""),lwd=2,xlab="Time (days)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
    par(mar = c(5.1,4.5,7.1,2.1))
    plot(mfit.subtype1,main=paste("Survival Prediction including ",labcohort," based on \n",stitle, genlab," with ",mesgender," ",nbobs, " patients.",sep=""),lwd=3,xlab="Time (months)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=2,col.main="blue", cex.axis=1.7, cex.lab=1.7)
    
    ## mtext(paste("Log-rank P = ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1)
    
    # mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
    # mtext(paste("HR 2groups: ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep=""),side=1,line=-2,adj=0,cex=1.2)
    # mtext(paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep=""),side=1,line=-3,adj=0,cex=1.2)
    
    #  mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
    mtext((paste("HR continuous: ",round(summary(c0)$coefficients[2],2)," p-val: ",formatC(summary(c0)$coefficients[5], format = "e", digits = 2),sep="")),side=1,line=-2,adj=0.05,cex=1.7,col="red")
    # if (issue)  { 
    #   mtext((mess),side=1,line=-2,adj=0.05,cex=1.7,col="black")
    #   }
    mtext((paste("HR 2 groups: ",round(summary(c)$coefficients[2],2)," p-val: ",formatC(summary(c)$coefficients[5],format = "e", digits = 2 ),sep="")),side=1,line=-4,adj=0.05,cex=1.7,col="red") 

    mtext((paste("Hazard Ratio (HR): ")),side=1,line=-6,adj=0.05,cex=1.7,col="red")
     
    
    # legend('topright', c("Q1","Q4"),lty=1, col=c("blue","red"), bty='n', cex=.75,title=paste(comgenes,collapse=","))
    if (tlabchange==F) 
    {  if(optsurv!="mut") 
       {
      stitle21 = "No   : "
      stitle22 = "Yes  : " } 
      else 
      {
        stitle21 = "Wt   : "
        stitle22 = "Mut  : " 
        }
    } else 
      { 
      if (optgroup=="25") 
      {
        stitle21 = "Q1 : "
        stitle22 = "Q4 : "
      } 
      else {
          if (optgroup=="50")
          {
            stitle21 = "Lower than Median: "
            stitle22 = "Higher than Median: "
          } else
          {
            stitle21 = "T1 : "
            stitle22 = "T3 : "
          }
        }
    
      # stitle21 = "Exp <= 25% percentile: "
      # stitle22 = "Exp  > 75%  percentile: "
    }
    ## legend('topright', c(paste0(stitle21,nb1," patients"),paste0(stitle22,nb2," patients")),lty=1, col=c("blue","red"), bty='n', cex=1.2,title=paste(comgenes,collapse=","))
    legend('topright', c(paste0(stitle21,nb1," patients"),paste0(stitle22,nb2," patients")),lty=1, col=c("blue","red"), bty='n', cex=1.7,title=genlab , inset=0.03)
    
    #######
    # if (issue) {
    # shinyalert::shinyalert(
    #   title = "Warning in Cox 2 groups",
    #   text = paste("Message: ",mess),
    #   type = "warning",
    #   imageHeight = 150,
    #   size ="m"
    # ) 
    # }
    ####################################################
    
  },  width=1400, height=800)
  
  ## new survplotReact with KM 06212024 ...................................................
  output$SurvplotReact_km <- renderPlot({
    
    resu <- findInfoExpReact()
    rcox <- resu$coxmulti ## new to add table in the same plot using gridExtra !!!! to do
    dat = resu$data
    cohort = as.character(input$onco3ch)
    optsurv = input$optsurv
    chkLasso = input$Lasso
    gender = input$optgender
    optgroup = input$optgroup
    # covariates = resu$comgenes
    comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
    
    tlabchange=F
    response.M25 = dat$variable
    ## 10/26
    shiny::validate(need(sd(response.M25, na.rm=T) > 0, 
                         "The gene values are the same"))
    # mutation2
    if (!((comgenes %in% c("CGA.STAIN","HORMONE.PROD","SYP.POSITIVE", "INSM1","MALEGENDER"))[1]) & (optsurv!="mut") ) {
      if (optgroup == "25") {
        response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[2],0,NA)
        response.M25[which(dat$variable>quantile(dat$variable,na.rm=T)[4])]=1
      } else
      {
        if (optgroup == "50") 
        {
          response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[3],0,NA)
          response.M25[which(dat$variable>quantile(dat$variable,na.rm=T)[3])]=1
        } else
        {
          response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T, probs = seq(0, 1, 0.33) )[2],0,NA)
          response.M25[which(dat$variable>quantile(dat$variable,na.rm=T, probs = seq(0, 1, 0.33))[3])]=1
        }
        
      }
      tlabchange = T
    } else
    {
      shiny::validate(need(nlevels(factor(response.M25)) > 1, 
                           "There are no two groups (eg. Wt vs Mut)"))
    }
    
    
    
    # nb1 = length(which(response.M25==0))
    # nb2 = length(which(response.M25==1))
    
    # new counts
    inb1 = which(response.M25==0)
    inb2 = which(response.M25==1)
    nb1 = length(which(!is.na(dat$months_to_last_known_alive[inb1])))
    nb2 = length(which(!is.na(dat$months_to_last_known_alive[inb2])))
    # 
    ## now using cox model
    c0=coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)
    nbobs= unlist(summary(c0)[4])
    
    issue = F ## new
    ## c=coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)
    c = tryCatch({
      coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)
    }, warning = function(w) {
      
      ## print(w)
      shinyalert::shinyalert(
        # title = "Warning in Cox 2 groups",
        # text = paste("Message: ",w),
        title = paste("Warning in Cox 2 groups:\n",w),
        text = "",
        type = "warning",
        imageHeight = 150,
        size ="m"
      ) 
      ##
      # mess <<- w
      # issue <<- T
    }, error = function(e) {
      shinyalert::shinyalert(
        # title = "Error in Cox 2 groups",
        # text = paste("Message: ",e),
        title = paste("Error in Cox 2 groups:\n",e),
        text = "",
        type = "error",
        imageHeight = 150,
        size ="m"
      ) 
      mess <<- e
      issue <<- T
    }
    )
    ## error case 
    ##
    if (class(c) != "coxph") c= coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)

    ##  mfit.subtype1=survfit(c,newdata=data.frame(response.M25=c(0,1)))
    if (length(comgenes) == 1) genaff = comgenes else genaff ="Average genes" 
    if (length(comgenes)<=10) genlab= paste(comgenes,collapse=",") else genlab = paste(paste(comgenes[1:10],collapse=","),paste(comgenes[11:length(comgenes)],collapse=","), sep="\n")
    
    labcohort = paste(cohort,collapse=",")
    ## improve cohort display
    ## integer div  % / % ,  mod %% >> rest 
    if (length(cohort)>10) {
      labcohort = ""
      nbl = length(cohort) %/% 10
      reste = length(cohort) %% 10
      for (k in seq(1,(nbl*10), by=10)) {
        labcohort = paste0(labcohort,paste(cohort[k:(k+9)],collapse=","),"\n")
      }
      if (reste!=0) labcohort = paste0(labcohort,paste(cohort[(k+10):(k+9+reste)],collapse=","))
    }
    
    ## end improvments
    
    # if (optsurv=="xsq") {
    #   if (chkLasso) { stitle = " Average Expression of Lasso selected genes: " } else {
    #     if (length(comgenes)==1) stitle = " gene expression of "  else stitle = " Average Expression of entered genes: " }
    # }
    # else
    # {
    #   if (optsurv=="mda") 
    #   {stitle = " metadata using " } else 
    #   {stitle = " mutation of " }
    # }
    
    if (optsurv=="xsq") {
      if (chkLasso) { stitle = " Average Expression of Lasso selected genes: " } else {
        if (length(comgenes)==1) stitle = " gene expression of "  else stitle = " Average Expression of entered genes: " }
    }
    else
    {
      if (optsurv=="mda") 
        stitle = " metadata "  else 
        { if (optsurv=="mut")
          stitle = " mutation of " 
        else
        { if (optsurv=="sig")
          stitle = " signature "  else 
            stitle = " NMF cluster "
        }
        }
    }

    if (gender=="Both") mesgender = "All" else mesgender = gender
    ###
    if (tlabchange==F) 
    {  if(optsurv!="mut") 
    {
      vleg = c("No","Yes")
      # stitle21 = "No   : "
      # stitle22 = "Yes  : " 
      } 
      else 
      {
        vleg=c("Wt","Mut")
        # stitle21 = "Wt   : "
        # stitle22 = "Mut  : " 
      }
    } else 
    { 
      if (optgroup=="25") 
      {
        vleg = c("Q1","Q4")
        # stitle21 = "Q1 : "
        # stitle22 = "Q4 : "
      } 
      else {
        if (optgroup=="50")
        {
          vleg = c("Lower than Median","Higher than Median")
         # stitle21 = "Lower than Median: "
        #  stitle22 = "Higher than Median: "
        } else
        {
          vleg = c("T1","T3")
          # stitle21 = "T1 : "
          # stitle22 = "T3 : "
        }
      }
      
    }
    
    ## new months_to_last_known_alive
    mydata2 = data.frame(vm= dat$months_to_last_known_alive, vs= dat$vital_status, response.M25)
    fit2 <- survfit(Surv(vm, vs) ~ response.M25, data=mydata2 )
    # finally visualize the survival model using ggsurvplot.
    ggsurv <- survminer::ggsurvplot(fit=fit2, data=mydata2, pval=T,pval.method=T,
                                   pval.size = 7,
                                   title = paste("Survival prediction in ",labcohort," based on",stitle,"\n" ,genlab," with ",mesgender," ",nbobs, " patients.",sep=""),
                                   subtitle = paste0("Hazard Ratio (HR) 2 groups:    ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4)," and HR continuous:    ",round(summary(c0)$coefficients[2],4),' p-val: ',round(summary(c0)$coefficients[5],4)),
                                   risk.table.col = "strata",
                                   ## xlim = c(0,500),
                                   ## legend.labs = c("low","high"),
                                   legend.labs = vleg,
                                   xlab = "Time in months",
                                   risk.table=T,
                                   legend.title = genaff, 
                                   surv.median.line = "hv",
                                   # conf.int=T,
                                   palette = c("blue","red"),
                                   size = 1, ## size of line
                                   font.x = c(20, face = "bold"),
                                   font.y = c(20, face = "bold"),
                                   font.tickslab = c(20),
                                   font.legend = c(20),
                                   fontsize = 8, ## inside table
                                   # tables.theme = theme_survminer() + theme(plot.title = element_text(size = 40)),
                                   tables.theme = theme_survminer(font.main = 20, font.x = 20, font.y = 20, font.tickslab = 20),
                                   ## font.tickslab = c(12),
                                   ggtheme = theme_bw() + theme(plot.title = element_text(size = 28, hjust = 0.5, face = "bold" , family = "Helvetica"), 
                                                                plot.subtitle = element_text(size = 20, hjust = 0.5,  family = "Helvetica")) 
                                    ) # 
    
 ## ggsurv$plot <- ggsurv$plot+ ggplot2::annotate("text",x=30, y = 0.75, label = paste0("Hazard Ratio (HR): \n HR 2 groups:    ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),"\n HR continuous:    ",round(summary(c0)$coefficients[2],4),' p-val: ',round(summary(c0)$coefficients[5],4)), size = 7)
 ggsurv
    
    ##
    
    # par(mar = c(5.1,4.5,7.1,2.1))
    # plot(mfit.subtype1,main=paste("Survival Prediction including ",labcohort," based on \n",stitle, genlab," with ",mesgender," ",nbobs, " patients.",sep=""),lwd=3,xlab="Time (months)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=2,col.main="blue", cex.axis=1.7, cex.lab=1.7)
    # mtext((paste("HR continuous: ",round(summary(c0)$coefficients[2],2)," p-val: ",formatC(summary(c0)$coefficients[5], format = "e", digits = 2),sep="")),side=1,line=-2,adj=0.05,cex=1.7,col="red")
    # mtext((paste("HR 2 groups: ",round(summary(c)$coefficients[2],2)," p-val: ",formatC(summary(c)$coefficients[5],format = "e", digits = 2 ),sep="")),side=1,line=-4,adj=0.05,cex=1.7,col="red") 
    # mtext((paste("Hazard Ratio (HR): ")),side=1,line=-6,adj=0.05,cex=1.7,col="red")
    
    
  },  width=1400, height=800)
  
  ## end KM 06212024
  
  output$CoxDetailsReact <- DT::renderDataTable({
    resu <- findInfoExpReact()$coxmulti
    comgenes <- findInfoExpReact()$comgenes
    resu = as.matrix(resu)
    resu[,1:2] = round(resu[,1:2],2)
    resu[,3] = as.numeric(formatC(resu[,3], format = "e", digits = 2))
    ## 
    rownames(resu) = substr(rownames(resu),4,nchar(rownames(resu)))
    if (input$optsurv == "xsq") {
    avgres = findInfoExpReact()$coxaverage
    ## avgres = round(as.matrix(avgres),4)
    avgres = as.matrix(avgres)
    avgres[,1:2] = round(avgres[,1:2],2)
    avgres[,3] = as.numeric(formatC(avgres[,3], format = "e", digits = 2))
    if (length(comgenes) > 1 ) resu = rbind(resu,avgres) 
    }
    
    if (input$Lasso) selsource = "Lasso Cox" else selsource ="Cox Multivariate"
    colnames(resu)[2] = "Hazard Ratio (HR)"
    resu = resu[,c(2,3,1), drop =F]
    DT::datatable(resu, filter='top', selection = "none", style='bootstrap4', rownames = T,caption=htmltools::tags$caption(paste0("Genes Hazard ratio using ",selsource),style="color:dodgerblue; font-size: 18px") )
    
  })
  
  output$CoxHeatmapReact <- renderPlotly({
    resu <- findInfoExpReact()
    mydata <- resu$dat
    nbgenes = length(resu$comgenes)
    mydata = mydata[,c((nbgenes+1), 1:nbgenes)]
    colnames(mydata)[1] = "Overall.Survival"
    mydata <- mydata[order(mydata[,1], decreasing = TRUE), ]
    mydata <- na.omit(mydata)
    mydata = t(mydata)
    
    
    scaledDataMatrix <- scaleDataForHeatmap(mydata, TRUE)
    
    # #save(scaledDataMatrix, file = "~/Downloads/scaledDataMatrix.RData")
    # 
    # 
    dataMatrix2=matrix(paste("original_value:",mydata), nrow=nrow(mydata),ncol=ncol(mydata))
     g1 <-heatmaply::heatmaply(x = scaledDataMatrix,custom_hovertext = dataMatrix2, grid_gap = 0.5, cellnote_color="black",
                              colors= colorpanel(75,low="blue",mid="white",high="red"),fontsize_col = 10 ,fontsize_row = 14,
                              dendrogram="none",label_names=c("row","column","scaled_value"))   %>% plotly::layout(margin=list(t = 10))
     g2 <- config(p = g1,  cloud=FALSE, displaylogo=FALSE, displayModeBar=TRUE,
                 modeBarButtonsToRemove=c("select2d", "sendDataToCloud", "pan2d", "resetScale2d",
                                          "hoverClosestCartesian", "hoverCompareCartesian",
                                          "lasso2d", "zoomIn2d", "zoomOut2d","autoScale2d","toggleSpikelines","zoom2d"))
    g2
  })
## ----------------------------------------------------------------------
 output$Survplot <- renderPlot({
   cohort   <- as.character(input$onco3ch)
   ## cat(cohort ,"...\n")
   geneexp  <- srcContent[[input$survSource]][["molPharmData"]][[input$optsurv]] 
   clin     <- srcContent[[input$survSource]][["sampleData"]]
   ## cat(input$maxSurvPredictorsLasso,"  maxsurv ..\n")
   # drawPlotExp3(cohort, input$varname, geneexp, clin, input$sampletype,input$therapy,input$optsurv)
   ## drawPlotExp4(cohort, input$varname, geneexp, clin, input$sampletype,input$therapy,input$optsurv, input$Lasso,input$inputGSsurv ,input$maxSurvPredictorsLasso)
   ## drawPlotExp6(cohort, input$varname, geneexp, clin, input$sampletype,input$therapy,input$optsurv, input$Lasso,input$inputGSsurv ,input$maxSurvPredictorsLasso, input$optgender)
   
  ##  drawPlotExp5(cohort, input$varname, geneexp, clin, input$sampletype,input$therapy,input$optsurv, input$Lasso,input$inputGSsurv ,input$maxSurvPredictorsLasso, input$optgender)
   drawPlotExp_km(cohort, input$varname, geneexp, clin, input$sampletype,input$therapy,input$optsurv, input$Lasso,input$inputGSsurv ,input$maxSurvPredictorsLasso, input$optgender)
   
    },  width=1400, height=800)
 
 output$CoxDetails <- renderPrint({
   cohort   <- as.character(input$onco3ch)
   geneexp  <- srcContent[[input$survSource]][["molPharmData"]][[input$optsurv]] 
   clin     <- srcContent[[input$survSource]][["sampleData"]]
   resu <- findInfoExp2(toupper(input$varname), cohort, geneexp,clin,input$sampletype,input$therapy,input$optsurv)
   dat = resu$data
   comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
   # write.csv(dat,"temp.csv")
   # cat(dim(dat), "\n")
   shiny::validate(need(nrow(dat)>0, 
                        "data not found"))
   
   shiny::validate(need(length(na.omit(dat$normalized_count))>0, 
                        "all values are NA"))
   cat("Summary continuous Cox results for ",comgenes,"\n")
   summary(coxph(Surv(days_to_last_known_alive, vital_status)~normalized_count,data=dat))
 })
 
 output$CoxDetails3 <- renderPrint({
   cohort   <- as.character(input$onco3ch)
   geneexp  <- srcContent[[input$survSource]][["molPharmData"]][[input$optsurv]] 
   clin     <- srcContent[[input$survSource]][["sampleData"]]
   resu <- findInfoExp3(toupper(input$varname), cohort, geneexp,clin,input$sampletype,input$therapy,input$optsurv)
   dat = resu$data
   comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
   # write.csv(dat,"temp.csv")
   # cat(dim(dat), "\n")
   cat("Notes: Event means Death, the beta value is coef and the Hazard ratio is exp(coef) \n")
   cat("------------------\n")
  if (length(comgenes)==1) {
   cat("Summary continuous Cox results for ",comgenes,"\n")
     summary(coxph(Surv(days_to_last_known_alive, vital_status)~variable,data=dat)) 
   }
   else {
     cat("Summary Cox Multivariate based for",comgenes,"\n")
     ## print(summary(coxph(Surv(days_to_last_known_alive, vital_status)~.-variable,data=dat)))
     print(summary(resu$coxmulti))
     cat("------------------\n")
     cat("Summary continuous Cox results based on risk scores for ",comgenes,"\n")
     print(summary(coxph(Surv(days_to_last_known_alive, vital_status)~variable,data=dat)))
   }
     
 })

 output$CoxDetails4 <- renderPrint({
   cohort   <- as.character(input$onco3ch)
   geneexp  <- srcContent[[input$survSource]][["molPharmData"]][[input$optsurv]] 
   clin     <- srcContent[[input$survSource]][["sampleData"]]
   
   ## resu <- findInfoExp4(toupper(input$varname), cohort, geneexp,clin,input$sampletype,input$therapy,input$optsurv, input$Lasso,input$inputGSsurv ,input$maxSurvPredictorsLasso)
   
   resu <- findInfoExp5(toupper(input$varname), cohort, geneexp,clin,input$sampletype,input$therapy,input$optsurv, input$Lasso,input$inputGSsurv ,input$maxSurvPredictorsLasso, input$optgender)
   dat = resu$data
   comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
   # write.csv(dat,"temp.csv")
   
   percentile.25.75 = dat$variable
   if (!((comgenes %in% c("CGA.STAIN","HORMONE.PROD","SYP.POSITIVE", "INSM1","MALEGENDER"))[1]) & (input$optsurv!="mut") ) {
     percentile.25.75=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[2],0,NA)
     percentile.25.75[which(dat$variable>quantile(dat$variable,na.rm=T)[4])]=1
     
   } 
   
   # cat(dim(dat), "\n")
   cat("Notes: Event means Death, the beta value is coef and the Hazard ratio is exp(coef) \n")
   cat("----------------------------------------------------------------------------------\n")
   if (input$Lasso) {
     cat("Lasso selected genes: ",comgenes,"\n")
     ## print(summary(coxph(Surv(days_to_last_known_alive, vital_status)~.-variable,data=dat)))
     print(resu$coxmulti)
     cat("----------------------------------------------------------------------------------\n")
     cat("Summary continuous Cox results based on Lasso average gene expression for ",comgenes,"\n")
     print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)))
     cat("----------------------------------------------------------------------------------\n")
     cat("Summary 2 groups Cox results based on Lasso average gene expression for ",comgenes,"\n")
     print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~percentile.25.75,data=dat)))
   }
   else {
   if (length(comgenes)==1) {
     cat("Summary continuous Cox results for ",comgenes,"\n")
     print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)))
     cat("----------------------------------------------------------------------------------\n")
     cat("Summary 2 groups Cox results for ",comgenes,"\n")
     print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~percentile.25.75,data=dat)))
   }
   else {
     cat("Summary Cox Multivariate based for",comgenes,"\n")
     ## print(summary(coxph(Surv(days_to_last_known_alive, vital_status)~.-variable,data=dat)))
     print(summary(resu$coxmulti))
     cat("----------------------------------------------------------------------------------\n")
     cat("Summary continuous Cox results based on average gene expression for ",comgenes,"\n")
     print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)))
     cat("----------------------------------------------------------------------------------\n")
     cat("Summary 2 groups Cox results based on average gene expression for ",comgenes,"\n")
     print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~percentile.25.75,data=dat)))
   }
   }
   
 })
 
 output$CoxPlusDetailsReact <- renderPrint({
   
   ## resu <- findInfoExp4(toupper(input$varname), cohort, geneexp,clin,input$sampletype,input$therapy,input$optsurv, input$Lasso,input$inputGSsurv ,input$maxSurvPredictorsLasso)
   
   resu <- findInfoExpReact()
   ## dat = resu$data
   comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
   
   
   # cat(dim(dat), "\n")
   cat("Notes: Event means Death, the beta value is coef and the Hazard ratio is exp(coef) \n")
   cat("----------------------------------------------------------------------------------\n")
   if (input$Lasso & input$optsurv == "xsq") {
     cat("Lasso selected genes: ",comgenes,"\n")
     ## print(summary(coxph(Surv(days_to_last_known_alive, vital_status)~.-variable,data=dat)))
     print(resu$coxmulti)
     # cat("----------------------------------------------------------------------------------\n")
     # cat("Summary continuous Cox results based on Lasso average gene expression for ",comgenes,"\n")
     # print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)))
     # cat("----------------------------------------------------------------------------------\n")
     # cat("Summary 2 groups Cox results based on Lasso average gene expression for ",comgenes,"\n")
     # print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~percentile.25.75,data=dat)))
   }
   else {
     if (length(comgenes)==1) {
       cat("Summary continuous Cox results for ",comgenes,"\n")
       print(summary(resu$coxmodel))
       # print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)))
       cat("\n----------------------------------------------------------------------------------\n")
       cat("Checking for proportional assumption: violation if pvalue (p) is significant\n")
       cat("----------------------------------------------------------------------------------\n")
      
       ##### catch error
       issue = F ## new
       czph = tryCatch({
         cox.zph(resu$coxmodel)
       }, warning = function(w) {
           mess <<- w
           issue <<- T
       }, error = function(e) {
           mess <<- e
           issue <<- T
       }
       )
       if (issue) print(mess) else print(czph)
       
      #####
       ## print(cox.zph(resu$coxmodel)) ## old code
     }
     else {
       cat("Summary Cox Multivariate based on",comgenes,"\n")
       ## print(summary(coxph(Surv(days_to_last_known_alive, vital_status)~.-variable,data=dat)))
       print(summary(resu$coxmodel))
       cat("\n----------------------------------------------------------------------------------\n")
       cat("Checking for proportional assumption: violation if pvalue (p) is significant\n")
       cat("----------------------------------------------------------------------------------\n")
       ##### catch error
       issue = F ## new
       czph = tryCatch({
         cox.zph(resu$coxmodel)
       }, warning = function(w) {
         mess <<- w
         issue <<- T
       }, error = function(e) {
         mess <<- e
         issue <<- T
       }
       )
       if (issue) print(mess) else print(czph)
       
       #####
       ### print(cox.zph(resu$coxmodel))
       
       # cat("----------------------------------------------------------------------------------\n")
       # cat("Summary continuous Cox results based on average gene expression for ",comgenes,"\n")
       # print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)))
       # cat("----------------------------------------------------------------------------------\n")
       # cat("Summary 2 groups Cox results based on average gene expression for ",comgenes,"\n")
       # print(summary(coxph(Surv(months_to_last_known_alive, vital_status)~percentile.25.75,data=dat)))
     }
   }
   
 })
 
	#----[Reactive Variables]---------------------------------------------------------------
	# Record current input validity status, data type prefix values.
 
	## globalReactiveValues <- reactiveValues(xPrefix = NULL, yPrefix = NULL)
  globalReactiveValues <- 	 reactiveValues(xPrefix = NULL, yPrefix = NULL, xPrefixpb = NULL, yPrefixpb = NULL)
 
	# Provides a srcContent (list-based) data structure containing all molecular profiling
	# drug response, and feature/sample annotation data required for the application 
	# (for data sources specified in the config.json file).
	srcContentReactive <- reactive({
		if(!exists("srcContent")) {
			srcContent <- lapply(config, loadSourceContent)
			isLoadedSrc <- vapply(srcContent, function(x) { !is.null(x) }, logical(1))
			if (any(!isLoadedSrc)){
				srcContent <- srcContent[isLoadedSrc]
			}

			# For NCI-60, replace default color map to use CellMiner tissue type colors.
			nci60ColorTab <- loadNciColorSet(returnDf=TRUE)
			nci60ColorTab$OncoTree1 <- srcContent$nci60$sampleData$OncoTree1
			srcContent$nci60$tissueColorMap <- c(by(nci60ColorTab, nci60ColorTab$OncoTree1,
																							FUN = function(x) unique(x$colors)))
		}
		
		return(srcContent)
	})
	
	isPackageLoadingComplete <- reactive({
		srcContentReactive()
		
		return(TRUE)
	})

	# Provides the set of all possible feature identifiers available for plotting along
	# the x-axis of the 2-D scatter plot (and for use in pattern comparisons).
	# xIdChoices <- reactive({
	# 	srcContent <- srcContentReactive()
	# 
	# 	t1 <- lapply(srcContent[[input$xDataset]][["molPharmData"]], function(x) {
	# 		return(unname(removeMolDataType(rownames(x))))
	# 	})
	# 	l1 <- unique(unname(unlist(t1)))
	# 
	# 	t2 <- srcContent[[input$xDataset]][["drugInfo"]]
	# 	l2 <- unname(removeMolDataType(rownames(t2)))
	# 
	# 	l3 <- c(l1, l2)
	# 
	# 	return(l3)
	# })

	# Provides the set of all possible feature identifiers available for plotting along
	# the y-axis of the 2-D scatter plot.
	# yIdChoices <- reactive({
	# 	srcContent <- srcContentReactive()
	# 
	# 	t1 <- lapply(srcContent[[input$yDataset]][["molPharmData"]], function(x) {
	# 		return(unname(removeMolDataType(rownames(x))))
	# 	})
	# 	l1 <- unique(unname(unlist(t1)))
	# 
	# 	t2 <- srcContent[[input$yDataset]][["drugInfo"]]
	# 	l2 <- unname(removeMolDataType(rownames(t2)))
	# 
	# 	l3 <- c(l1, l2)
	# 
	# 	return(l3)
	# })
	
	
	# Returns all valid OncoTree types with respect to
	# the xDataset cell line OncoTree types AND user-selected tissue
	# types for inclusion or exclusion. These types need not be
	# mutually exclusive, given the nested OncoTree structure.
	# NOTE: with added data packages, always verify that matched cell
	# line OncoTree tissue type annotations are consistent across packages.
	selectedTissues <- reactive({
	  tree <- input$tree
	  req(tree)
	  ## zz = get_selected(tree, format = "names")
	  selectedTissues = names(unlist(get_selected(tree, format = "slices")))
	  selectedTissues = gsub("\\.",":", selectedTissues)
	  selectedTissues = gsub("all:","", selectedTissues)
	  selectedTissues = gsub("none:","", selectedTissues)
	  cat(selectedTissues,"\n")
	  return(sort(selectedTissues))
	  
	})
	
	analysisTissueTypes <- reactive({
		# Note: We want this code to re-run whenever either 
		# input$tissueSelectionMode OR  input$selectedTissues change.
		# BUT, *reactivity can be based on the selectedTissues alone*,
		# because switching the tissueSelectionMode will always trigger
		# a change to a distinct (default) selectedTissues value (e.g.,
		# "all" in the case of "Include", or "none" in the case of 
		# "Exclude"). 
		# Without the isolate() around input$tisssueSelectionMode (below),
		# this code actually runs twice upon a change to the
		# tissueSelectionMode, with the first run having an invalid
		# selectedTissues value (because the reativity relative
		# to the tissueSelectionMode appears to trigger the code below
		# before the selectedValues changes, for the first run). Then
		# the update of the selectedValues to the default value triggers
		# a second run. This behavior causes a bug-like re-drawing of
		# the 2D plot, etc.)
		tissueSelectionMode <- isolate(input$tissueSelectionMode)
		## new stuff ----------------------------------------------------------
		# tree <- input$tree
		# req(tree)
		#           ## zz = get_selected(tree, format = "names")
		# selectedTissues = names(unlist(get_selected(tree, format = "slices")))
		# selectedTissues = gsub(".",":", selectedTissues)
		
		selectedTissues = selectedTissues()
		## end new stuff ------------------------------------------------------
		
		# selectedTissues <- input$selectedTissues
		
		# cat("--- Entering analysisTissueTypes()", sep = "\n")
		# cat(paste0("Selection Mode: ", tissueSelectionMode), sep = "\n")
		# cat(paste0("Selected Tissues: ", selectedTissues), sep = "\n")
		
		srcContent <- srcContentReactive()
		tissueToSamplesMap <- srcContent[[input$xDataset]][["tissueToSamplesMap"]]
		tissueTypes <- names(tissueToSamplesMap)
		
		if (tissueSelectionMode == "To include"){
			if (!("all" %in% selectedTissues)){
				selectedLines <- unique(c(tissueToSamplesMap[selectedTissues], recursive = TRUE))
				# For which tissue types are ALL lines in selectedLines?
				allInSelectedLines <- vapply(tissueToSamplesMap, function(x){
					all(x %in% selectedLines)
				}, logical(1))
				tissueTypes <- names(tissueToSamplesMap[allInSelectedLines])
			}
		} else{ # tissueSelectionMode == "Exclude"
			if (!("none" %in% selectedTissues)){
				selectedLines <- unique(c(tissueToSamplesMap[selectedTissues], recursive = TRUE))
				# For which tissue types are NO lines in selectedLines?
				notInSelectedLines <- vapply(tissueToSamplesMap, function(x){
					length(intersect(x, selectedLines)) == 0
				}, logical(1))
				tissueTypes <- names(tissueToSamplesMap[notInSelectedLines])
			}
		}
		
		#cat("--- LEAVING analysisTissueTypes()", sep = "\n")
		
		return(sort(unique(tissueTypes)))
	})
	
	# Provides a data frame with columns indicating the matched cell lines between
	# the input$xDataset (column 1) and the input$yDataset (column 2).
	# matchedCellLinesTab <- reactive({
	# 	shiny::validate(need(require(rcellminerUtils),
	# 											 "ERROR: x and y axis data sets must be the same."))
	# 	matchedCellLinesTab <- getMatchedCellLines(c(input$xDataset, input$yDataset))
	# 	shiny::validate(need(nrow(matchedCellLinesTab) > 0, 
	# 											 "There are no shared cell lines between the selected datasets."))
	# 	colnames(matchedCellLinesTab) <- c("xDataset", "yDataset")
	# 	return(matchedCellLinesTab)
	# })
	
	# Provides a data frame with columns indicating the matched cell lines between
	# the input$xDataset (column 1) and the input$yDataset (column 2).
	# The cell line pairing will be updated to reflect restrictions based on:
	# --- matched cell lines across databases (if input$xDataset != input$yDataset)
	# --- user tissue type selections.
	matchedCellLinesTab <- reactive({
		srcContent <- srcContentReactive()
		analysisTissueTypes <- analysisTissueTypes()
		
		if (input$xDataset == input$yDataset){
			matchedCellLinesTab <- data.frame(
				xDataset = srcContent[[input$xDataset]]$sampleData[, "Name"],
				stringsAsFactors = FALSE
			)
			matchedCellLinesTab$yDataset <- matchedCellLinesTab$xDataset
		} else{
			shiny::validate(need(require(rcellminerUtilsCDB),
													 "ERROR: x and y axis data sets must be the same."))
			matchedCellLinesTab <- getMatchedCellLines(c(input$xDataset, input$yDataset))
			shiny::validate(need(nrow(matchedCellLinesTab) > 0, 
													 "There are no shared cell lines between the selected datasets."))
			colnames(matchedCellLinesTab) <- c("xDataset", "yDataset")
		}
		stopifnot(all(!duplicated(matchedCellLinesTab$xDataset)))
		rownames(matchedCellLinesTab) <- matchedCellLinesTab$xDataset
		
		tissueMatchedLines <- getTissueTypeSamples(analysisTissueTypes, input$xDataset, srcContent)
		tissueMatchedLines <- intersect(tissueMatchedLines, rownames(matchedCellLinesTab))
		
		shiny::validate(need(length(tissueMatchedLines) > 0, 
												 "There are no cell lines of the selected tissue type(s)."))
		
		matchedCellLinesTab <- matchedCellLinesTab[tissueMatchedLines, ]
		
		return(matchedCellLinesTab)
	})
	
	##------------
	PatternCompTable <- reactive({
	  srcContent <- srcContentReactive()
	  ## new
	  if (input$crossdb == "Yes")
	  {
	    # if (input$patternComparisonSeed == "xPattern"){
	      dat <- xData()
	      pcDataset <- input$yDataset
	      # selectedLines <- names(yData()$data)
	      selectedLines <- as.character(matchedCellLinesTab()[, "yDataset"])
	    # } else{
	    #   dat <- yData()
	    #   pcDataset <- input$xDataset
	    #   # selectedLines <- names(xData()$data)
	    #   selectedLines <- as.character(matchedCellLinesTab()[, "xDataset"])
	    # }
	    names(dat$data) = selectedLines
	  }
	  
	  else
	  {
	  
	  # if (input$patternComparisonSeed == "xPattern"){
	    dat <- xData()
	    pcDataset <- input$xDataset
	  # } else{
	  #   dat <- yData()
	  #   pcDataset <- input$yDataset
	  # }
	    selectedLines <- names(dat$data)
	  } # end new
	  
	  
	  shiny::validate(need(length(selectedLines)>0, paste("ERROR:", " No common complete data found.")))
	  shiny::validate(need(length(selectedLines)>2, paste("ERROR:", " No display for less than 3 observations.")))
	  
	  if(input$patternComparisonType == "drug") {
	    shiny::validate(need(srcContent[[pcDataset]][["molPharmData"]][["act"]], "No drug available for this cell line set"))
	    #if (is.null(srcContent[[pcDataset]][["molPharmData"]][["act"]])) stop("No drug available for this cell line set")
	    results <- patternComparison(dat$data,
	                                 srcContent[[pcDataset]][["molPharmData"]][["act"]][, selectedLines, drop=FALSE])
	    results$ids <- rownames(results)
	    results$NAME <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "NAME"]
	    # OLD
	    # if ("MOA" %in% colnames(srcContent[[pcDataset]][["drugInfo"]])){
	    #   results$MOA <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "MOA"]
	    #   results <- results[, c("ids", "NAME", "MOA", "COR", "PVAL")]
	    #   colnames(results) <- c("ID", "Name", "MOA", "Correlation", "P-Value")
	    # } else{
	    #   results <- results[, c("ids", "NAME", "COR", "PVAL")]
	    #   colnames(results) <- c("ID", "Name", "Correlation", "P-Value")
	    # }
	    
	    #NEW
	    results$MOA <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "MOA"]
	    results$CLINICAL.STATUS <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "CLINICAL.STATUS"]
	    results <- results[, c("ids", "NAME", "MOA", "CLINICAL.STATUS", "COR", "PVAL")]
	    colnames(results) <- c("ID", "Name", "MOA", "CLINICAL.STATUS", "Correlation", "P-Value")
	    # end new
	    
	    DataType <- getMolDataType(results$ID)
	    DrugID <- removeMolDataType(results$ID)
	    results$ID <- NULL
	    results$FDR=p.adjust(results[,"P-Value"],method="BH",nrow(results))
	    results=cbind(DataType,DrugID,results)
	    colnames(results)[1:2]=c("Data Type", "Drug ID")
	    
	  } else {
	    molPharmData <- srcContent[[pcDataset]][["molPharmData"]]
	    molData <- molPharmData[setdiff(names(molPharmData), c("act","copA","mutA","metA","expA","xaiA","proA","mirA","mdaA","swaA","xsqA","mthA","hisA","criA","mtbA","rrbA","tmmA","tpmA","varA","var","mucA","hllA","sigA","nmfA"))]
	    shiny::validate(need(length(molData)>0, "No molecular data available for this cell line set"))
	    ##if (length(molData)==0) stop("No molecular data available for this cell line set")
	    ## old: molData <- lapply(molData, function(X) X[, selectedLines])
	    ## new selection in case a dataset has only one row
	    ## one solution: molData <- lapply(molData, function(X) subset(X, select=selectedLines))
	    molData <- lapply(molData, function(X) X[, selectedLines,drop=FALSE])
	    results <- patternComparison(dat$data, molData)
	    results$ids <- rownames(results)
	    
	    results$molDataType <- getMolDataType(results$ids)
	    results$gene <- removeMolDataType(results$ids)
	    
	    # Reorder columns
	    results <- results[, c("ids", "molDataType", "gene", "COR", "PVAL")]
	    colnames(results) <- c("ID", "Data Type", "Gene", "Correlation", "P-Value")
	    
	    if (require(rcellminerUtilsCDB)){
	      chromLocs <- character(nrow(results))
	      haveLoc <- results$Gene %in% names(geneToChromBand)
	      chromLocs[haveLoc] <- geneToChromBand[results$Gene[haveLoc]]
	      
	      results$Location <- chromLocs
	      results <- results[, c("ID", "Data Type", "Gene", "Location", "Correlation", "P-Value")]
	    }
	    results$FDR=p.adjust(results[,"P-Value"],method="BH",nrow(results))
	    
	    if (require(geneSetPathwayAnalysis)){
	      # old :results$Annotation <- geneSetPathwayAnalysis::geneAnnotTab[results$Gene, "SHORT_ANNOT"]
	      # issue related to prefix subsetting ex: "age" >> gives "ager"
	      #st=Sys.time()
	      results$Annotation <- geneSetPathwayAnalysis::geneAnnotTab[match(results$Gene,rownames(geneSetPathwayAnalysis::geneAnnotTab)), "SHORT_ANNOT"]
	      #et=Sys.time()
	      #cat(et-st,"\n")
	      results$Annotation[is.na(results$Annotation)] <- ""
	    }
	    
	    results$ID <- NULL
	    colnames(results)[2]="ID"
	  }
	  
	  results[, "Correlation"] <- round(results[, "Correlation"], 3)
	  results[, "P-Value"] <- signif(results[, "P-Value"], 3)
	  results[, "FDR"] <- signif(results[, "FDR"], 3)
	  ## sort by p-value
	  results <- results[order(results[, "P-Value"]),]
	  
	  return(results)
	})
	##--------------
	
	# Explanation of xData, yData reactive variables -------------------------------------------------
	# The xData and yData reactive variables provide list objects (accessed via xData() and yData()) 
	# that store the essential information about a data source feature that the application code
	# requires (e.g., to make plots, for pattern comparisons, etc.).
	# For example, if the x-axis feature is SLFN11 NCI-60 mRNA expression, the xData() 
	# list object would contain:
	# xData()$dataSource: "nci60"
	# xData()$name: "expSLFN11" (the feature identifier prepended with a data type specifier)
	# xData()$uniqName: "expSLFN11_nci60" (the above identifier appended with the data source)
	# xData()$plotLabel: "SLFN11 (exp, nci60)"
	# xData()$data: A vector of numeric feature data.
	#
	# Note: the reactive variable construction code ensures that:
	# (1) User-requested feature data is available for the specified data source,
	# (2) if the x-axis and y-axis features are derived from different data sources,
	#     the xData()$data and yData()$data vectors *have feature data from matched cell lines*.
	#
	# Through the reactivity, these variables are always 'current', reflecting the latest
	# user selections, with regard to feature identifiers, data source(s), etc.
	# ------------------------------------------------------------------------------------------------
	
	# Provides an a list object with x-axis feature-related data, including numeric data,
	# data type prefix, data source, and plot label.
	# Note that use of reactive matchedLinesTab() ensures that xData() and yData() are
	# current and coupled, reflecting matched cell lines (even when different x any y data 
	# sources are selected) of whatever tissue types are selected.
	xData <- reactive({
		## shiny::validate(need(length(input$selectedTissues) > 0, "Please select tissue types."))
	  shiny::validate(need(length(selectedTissues()) > 0, "Please select tissue types."))
    ## new
	  shiny::validate(need(!is.na(match(input$xPrefix,srcContentReactive()[[input$xDataset]][["featurePrefixes"]])), "Non valid data type"))
	  ##
		xPrefix <- input$xPrefix
		if (!is.character(xPrefix)){
			xPrefix <- srcContentReactive()[[input$xDataset]][["defaultFeatureX"]]
		}
		#
		originalId <- trimws(input$xId)
		# 
		xId <- getMatchedIds(xPrefix, trimws(input$xId), input$xDataset, srcContent = srcContentReactive())
		
		if (length(xId) == 0){
			shiny::validate(need(FALSE, paste("ERROR:", paste0(xPrefix, input$xId), "not found. Please use the Search IDs tab to find available IDs for each dataset.")))
		} else{
			globalReactiveValues$xPrefix <- xPrefix
			if (length(xId) > 1){
				warningMsg <- paste0("Other identifiers matching x-axis ID: ",
														 paste0(xId[-1], collapse = ", "), ".")
				showNotification(warningMsg, duration = 10, type = "message")
				xId <- xId[1]
			}
			# xData <- getFeatureData(xPrefix, xId, input$xDataset, srcContent = srcContentReactive())
			xData <- getFeatureData(xPrefix, xId, input$xDataset, srcContent = srcContentReactive(), originalId)
			
			matchedLinesTab <- matchedCellLinesTab()
			# xData$data <- xData$data[matchedLinesTab[, "xDataset"]]
			xData$data <- xData$data[as.character(matchedLinesTab[, "xDataset"])]
			
		}
		
		return(xData)
	})
	
	# Provides an a list object with y-axis feature-related data, including numeric data,
	# data type prefix, data source, and plot label.
	# Note that use of reactive matchedLinesTab() ensures that xData() and yData() are
	# current and coupled, reflecting matched cell lines (even when different x any y data 
	# sources are selected) of whatever tissue types are selected.
	yData <- reactive({
		## shiny::validate(need(length(input$selectedTissues) > 0, "Please select tissue types."))
	  shiny::validate(need(length(selectedTissues()) > 0, "Please select tissue types."))
	  ## new
	  shiny::validate(need(!is.na(match(input$yPrefix,srcContentReactive()[[input$yDataset]][["featurePrefixes"]])), "Non valid data type"))
	  ##
	  
		yPrefix <- input$yPrefix
		if (!is.character(yPrefix)){
			yPrefix <- srcContentReactive()[[input$yDataset]][["defaultFeatureY"]]
		}
		#
		originalId <- trimws(input$yId)
		# 
		
		yId <- getMatchedIds(yPrefix, trimws(input$yId), input$yDataset, srcContent = srcContentReactive())
		
		if (length(yId) == 0){
			shiny::validate(need(FALSE, paste("ERROR:", paste0(yPrefix, input$yId), "not found. Please use the Search IDs tab to find available IDs for each dataset.")))
		} else{
			globalReactiveValues$yPrefix <- yPrefix
			if (length(yId) > 1){
				warningMsg <- paste0("Other identifiers matching y-axis ID: ",
														 paste0(yId[-1], collapse = ", "), ".")
				showNotification(warningMsg, duration = 10, type = "message")
				yId <- yId[1]
			}
			# yData <- getFeatureData(yPrefix, yId, input$yDataset, srcContent = srcContentReactive())
			yData <- getFeatureData(yPrefix, yId, input$yDataset, srcContent = srcContentReactive(), originalId)
			
			matchedLinesTab <- matchedCellLinesTab()
			# yData$data <- yData$data[matchedLinesTab[, "yDataset"]]
			yData$data <- yData$data[as.character(matchedLinesTab[, "yDataset"])]
		}
		
		return(yData)
	})
	
	#----[outputs]--------------------------------------------------------------------------

  #----[Render Application Title]------------------------------------------------------
  output$public <- renderText("CellMiner")
  #------------------------------------------------------------------------------------

  #----[Render 2D Plot in 'Plot Data' Tab]---------------------------------------------
#   if(require(rCharts)) {
# 		output$rCharts <- renderChart({
# 			h1 <- makePlot(xData = xData(), yData = yData(), showColor = input$showColor,
# 										 showColorTissues = input$showColorTissues, dataSource = input$xDataset,
# 										 srcContent = srcContentReactive(), dom="rCharts")
# 		})
#   }

	# Alternative plotting
	output$rChartsAlternative <- renderPlotly({
		#-----[range check]----------------------------------------------------------
		# Note: Until a better solution can be found, these checks are needed.
		# The issue is that upon a data source change, there appears to be a moment 
		# when the xData() or yData() are updated, but the input$xAxisRange
		# or input$yAxisRange (from the sliderInput UI element) are not yet updated.
		# As such, the data value range can be out of synch with the invalidated
		# axis ranges. In the most extreme case, there are no points in the
		# specified range. The reactive code is quickly re-run with the proper
		# inputs, correcting the plot, but the error flashes briefly in a buggy 
		# looking way. 
		# Below we do a range check and quietly exit if something is amiss (knowing
		# that the reactivity will ensure that the code is re-run with a proper
		# set of inputs once thing settle down).
		#****************************************************************************
		xData <- xData()
		yData <- yData()
		shiny::validate(need(xData$uniqName != yData$uniqName, 
			"Please select distinct x and y axis variables."))
		
		xValRange <- range(xData$data, na.rm = TRUE)
		xLimits <- input$xAxisRange
		
		yValRange <- range(yData$data, na.rm = TRUE)
		yLimits <- input$yAxisRange
		
		# req(FALSE) causes immediate but quiet exit.
		req(!any(is.null(xValRange)), !any(is.null(xLimits)))
		req(!any(is.null(yValRange)), !any(is.null(yLimits)))
		req(!any(is.na(xValRange)), !any(is.na(xLimits)))
		req(!any(is.na(yValRange)), !any(is.na(yLimits)))
		# commented below
		# req((xLimits[1] <= xValRange[1]) && (xValRange[2] <= xLimits[2]))
		# req((yLimits[1] <= yValRange[1]) && (yValRange[2] <= yLimits[2]))
		
		cat("xAxis Limits: ", paste0(xLimits, collapse = " "), sep = "\n")
		cat("X_VAL_RANGE: ",  paste0(xValRange, collapse = " "), sep = "\n")
		
		cat("yAxis Limits: ", paste0(yLimits, collapse = " "), sep = "\n")
		cat("Y_VAL_RANGE: ",  paste0(yValRange, collapse = " "), sep = "\n")
		cat("-------------------------------------------", sep = "\n")
		#----------------------------------------------------------------------------
		
		p1 <- makePlotStatic(xData = xData, yData = yData, showColor = input$showColor, 
												 showColorTissues = showColorTissues(), dataSource = input$xDataset, 
												 xLimVals = xLimits, yLimVals = yLimits,
												 srcContent = srcContentReactive(),oncolor=oncolor, showCells=input$showCells)
		
		p1 <- p1 + theme(axis.text = element_text(size=16), plot.title = element_text(size = 16, hjust = 0.5), plot.subtitle = element_text(size = 16, hjust = 0.5), 
		                 axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) 
		p1 <- p1 + labs(subtitle = "salut")
		#theme_update(legend.position = c(0,0))
		g1 <- ggplotly(p1, width=plotWidth, height=plotHeight, tooltip=tooltipCol)
		#g1 <- layout(g1, margin=list(t = 75))
		g1 <- layout(g1, margin=list(t = 75), legend = list(font = list(size = 18)) )
		g2 <- config(p = g1,  cloud=FALSE, displaylogo=FALSE, displayModeBar=TRUE,
								 modeBarButtonsToRemove=c("select2d", "sendDataToCloud", "pan2d", "resetScale2d",
								 												 "hoverClosestCartesian", "hoverCompareCartesian",
								 												 "lasso2d", "zoomIn2d", "zoomOut2d"))
		g2
	})
	#--------------------------------------------------------------------------------------

  #----[Render Data Table in 'Download Data' Tab]----------------------------------------
  # Generate an HTML table view of the data
  output$table <- DT::renderDataTable({
  	# Column selection below is to restrict to cell line, x, y features,
  	# and tissue type information (source-provided + OncoTree).
  	dlDataTab <- getPlotData(xData = xData(), yData = yData(), showColor = input$showColor, 
  		showColorTissues = showColorTissues(), dataSource = input$xDataset, 
  		srcContent = srcContentReactive())
  
  	# cat(colnames(dlDataTab),"\n")
  	
  		shiny::validate(need(nrow(dlDataTab)>0, paste("ERROR:", " No common complete data found.")))
    # new stuff
  		dlDataTab[, 2] <- round(dlDataTab[, 2], 3)
  		dlDataTab[, 3] <- round(dlDataTab[, 3], 3)	
  		dlDataTabCols <- c(colnames(dlDataTab)[c(1,21,22,2:4)], paste0("OncoTree", 1:4)) # adding patient ID
  	#
  	## dlDataTabCols <- c(colnames(dlDataTab)[1:4], paste0("OncoTree", 1:4))
  	if ("EMT" %in% colnames(dlDataTab)) {
  	  dlDataTab[,"EMT"] <- gsub("EMT:","",dlDataTab[,"EMT"])
  		dlDataTabCols <- c(dlDataTabCols, "EMT")
  	}
  	
  	if ("priorTreatment" %in% colnames(dlDataTab)) {
  	  dlDataTab[,"priorTreatment"] <- gsub("PriorTreatment:","",dlDataTab[,"priorTreatment"])
  	  dlDataTabCols <- c(dlDataTabCols, "priorTreatment")
  	}
  	
  	if ("NAPY" %in% colnames(dlDataTab)) {
  	  dlDataTabCols <- c(dlDataTabCols, "NAPY")
  	}
  	# new stuff
  	if ("sampleType" %in% colnames(dlDataTab)) {
  	  dlDataTab[,"sampleType"] <- gsub("SampleType:","",dlDataTab[,"sampleType"])
  	  dlDataTabCols <- c(dlDataTabCols, "sampleType")
  	}
  	if ("biopsySite" %in% colnames(dlDataTab)) {
  	  dlDataTab[,"biopsySite"] <- gsub("BiopsySite:","",dlDataTab[,"biopsySite"])
  	  dlDataTabCols <- c(dlDataTabCols, "biopsySite")
  	}
  	if ("molecularSubtype" %in% colnames(dlDataTab)) {
  	  dlDataTab[,"molecularSubtype"] <- gsub("MolecularSubtype:","",dlDataTab[,"molecularSubtype"])
  	  dlDataTabCols <- c(dlDataTabCols, "molecularSubtype")
  	}
  	if ("dataSet" %in% colnames(dlDataTab)) {
  	  dlDataTab[,"dataSet"] <- gsub("DataSet:","",dlDataTab[,"dataSet"])
  	  dlDataTabCols <- c(dlDataTabCols, "dataSet")
  	}

  	dlDataTab <- dlDataTab[, dlDataTabCols]
  	# dlDataTab[, 2] <- round(dlDataTab[, 2], 3)
  	# dlDataTab[, 3] <- round(dlDataTab[, 3], 3)
  	colnames(dlDataTab)[1]="Patient_sample"
  	## new filtering
  	dlDataTab = dlDataTab[, -c(6,9,10,15)] # remove tissue, Onco3, Onco4, molecular subtype
  	DT::datatable(dlDataTab, rownames=FALSE, colnames=colnames(dlDataTab),
  								filter='top', style='bootstrap4', selection="none",
  								extensions = "FixedColumns",
  								## options=list(pageLength = nrow(dlDataTab), language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')))
  								options=list(scrollX = TRUE, fixedColumns = list(leftColumns = 1), lengthMenu = c(10, 25, 50, 100), pageLength = 10, language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')))
  								)
  })
	#--------------------------------------------------------------------------------------
	output$cortable <- DT::renderDataTable({
	  # Column selection below is to restrict to cell line, x, y features,
	  # and tissue type information (source-provided + OncoTree).
	  rescor= CorrelationTable(xData(),yData(),srcContentReactive())
	  
	  DT::datatable(rescor, rownames=FALSE, colnames=colnames(rescor),extensions='Buttons',
	                filter='top', style='bootstrap4', selection="none",
	                options=list(pageLength = nrow(rescor), language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipBt', buttons = list('copy', 'print', list(extend = 'collection',buttons = list(list(extend='csv',filename='tissue_correlation',title='Exported data from CellMinerCDB'), list(extend='excel',filename='tissue_correlation',title='Exported data from CellMinerCDB'), list(extend='pdf',filename='tissue_correlation',title='Exported data from CellMinerCDB')),text = 'Download'))))
	 
	  })
  #--------------------------------------------------------------------------------------
  #----[Render Data Table in 'Search IDs' Tab]-------------------------------------------
  # Generate an HTML table view of the data
  # Note: Searchable data is derived from the x-axis data source.
	output$ids <- DT::renderDataTable({
	  srcContent <- srcContentReactive()
	  drugIds   <- srcContent[[input$xDataset]][["drugInfo"]][, "ID"]
	  drugNames <- srcContent[[input$xDataset]][["drugInfo"]][, "NAME"]
	  moaNames  <- srcContent[[input$xDataset]][["drugInfo"]][, "MOA"]
	  exptype <- rep("act", length(drugIds))
	  
	  results <- data.frame(availableTypes=exptype, availableIds = drugIds, idNames=drugNames, properties=moaNames,
	                        stringsAsFactors=FALSE)
	  
	  # Make molecular data data.frame
	  molPharmData <- srcContent[[input$xDataset]][["molPharmData"]]
	  molData <- molPharmData[setdiff(names(molPharmData), "act")]
	  molDataIds <- as.vector(unlist(lapply(molData, function(x) { rownames(x) })))
	  
	  exptype2 <- substr(molDataIds,1,3)
	  molDataIds <- substr(molDataIds,4,nchar(molDataIds))
	  
	  molDataNames <- rep("", length(molDataIds))
	  moaNames <- rep("", length(molDataIds))
	  
	  tmp <- data.frame(availableTypes=exptype2,availableIds=molDataIds, idNames=molDataNames, properties=moaNames,
	                    stringsAsFactors=FALSE)
	  
	  # Join data.frames
	  results <- rbind(results, tmp)
	  
	  # Reverse Order/ no need for reverse for me
	  #results <- results[rev(rownames(results)),]
	  
	  colnames(results) <- c("Data type","ID ", "Drug Name", "Drug MOA")
	  selsource=metaConfig[[input$xDataset]][["fullName"]]
	  DT::datatable(results, rownames=FALSE, colnames=colnames(results),
	                filter='top', style='bootstrap4', selection = "none",
	                options=list(pageLength = 10, language=list(paginate = list(previous = 'Previous page', `next`= 'Next page'))), caption=htmltools::tags$caption(paste0("Ids table for ",selsource),style="color:dodgerblue; font-size: 18px"))
	})
		#--------------------------------------------------------------------------------------

	#----[Render Data Table in 'SearchIDs' Tab]-------------------------------------------
	# Generate an HTML table view of the data
	# Note: Searchable data is derived from the x-axis data source.
	output$ids_s <- DT::renderDataTable({
	  srcContent <- srcContentReactive()
	  if (input$dataTyp=="act") {
	    myframe=srcContent[[input$dataSrc]][["drugInfo"]]
	    
	  }
	  else
	  {
	    mytype=paste0(input$dataTyp,"A")
	    ## new 
	    shiny::validate(need(!is.null(srcContent[[input$dataSrc]][["molPharmData"]][[mytype]]), "Non valid data type")) ## for reactivity
	    if (ncol(srcContent[[input$dataSrc]][["molPharmData"]][[mytype]])==0) {
	    ## end
	    ## if (is.null(srcContent[[input$dataSrc]][["molPharmData"]][[mytype]]) | ncol(srcContent[[input$dataSrc]][["molPharmData"]][[mytype]])==0) {
	      ID=unlist(lapply(rownames(srcContent[[input$dataSrc]][["molPharmData"]][[input$dataTyp]]),function(x) {return(substr(x,4,nchar(x)))}))
	      myframe=data.frame(ID,stringsAsFactors = F)
	    }
	    else
	    { myframe=srcContent[[input$dataSrc]][["molPharmData"]][[mytype]]
	      #if (input$dataTyp=="swa") {
	      #   myframe=myframe[,c(2,1)]
	      # } else 
	        if (input$dataTyp=="exp" & input$dataSrc=="nciSclc")
	        { myframe=myframe[,c(4,3,5:13,1:2)] 
	        }
	      
	      colnames(myframe)[1]="ID" 
	    }
	  }
	  
	  
	  selsource=metaConfig[[input$dataSrc]][["fullName"]]
	  # DT::datatable(myframe, rownames=FALSE,extensions='Buttons',
	  #               filter='top', style='bootstrap4', selection = "none",
	  #               options=list(pageLength = 10,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipBt',buttons = list('copy', 'print', list(extend = 'collection',buttons = list(list(extend='csv',filename='search_id',title='Exported data from CellMinerCDB'), list(extend='excel',filename='search_id',title='Exported data from CellMinerCDB'), list(extend='pdf',filename='search_id',title='Exported data from CellMinerCDB')),text = 'Download')))
	  #               , caption=htmltools::tags$caption(paste0("Identifier search for ",selsource),style="color:dodgerblue; font-size: 18px")

	                DT::datatable(myframe, rownames=FALSE,extensions='Buttons',
	                              filter='top', style='bootstrap4', selection = "none",
	                              options=list(pageLength = 10,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipfBt',buttons = list('copy', 'print', list(extend = 'collection',buttons = list(list(extend='csv',filename='search_id',title='Exported data from CellMinerCDB'), list(extend='excel',filename='search_id',title='Exported data from CellMinerCDB'), list(extend='pdf',filename='search_id',title='Exported data from CellMinerCDB')),text = 'Download')))
	                              , caption=htmltools::tags$caption(paste0("Please use the boxes on top of the columns ","for specific search"),style="color:dodgerblue; font-size: 18px") 
	          
	                
	                
	   	)})

	## new version
	output$downloadPlotData2 <- DT::renderDataTable({
	  srcContent <- srcContentReactive()
	  drugIds   <- srcContent[[input$dataSrc]][["drugInfo"]][, "ID"]
	  drugNames <- srcContent[[input$dataSrc]][["drugInfo"]][, "NAME"]
	  moaNames  <- srcContent[[input$dataSrc]][["drugInfo"]][, "MOA"]
	  exptype <- rep("act", length(drugIds))
	  
	  results <- data.frame(availableTypes=exptype, availableIds = drugIds, idNames=drugNames, properties=moaNames,
	                        stringsAsFactors=FALSE)
	  
	  # Make molecular data data.frame
	  molPharmData <- srcContent[[input$dataSrc]][["molPharmData"]]
	  molData <- molPharmData[setdiff(names(molPharmData), "act")]
	  molDataIds <- as.vector(unlist(lapply(molData, function(x) { rownames(x) })))
	  
	  exptype2 <- substr(molDataIds,1,3)
	  molDataIds <- substr(molDataIds,4,nchar(molDataIds))
	  
	  molDataNames <- rep("", length(molDataIds))
	  moaNames <- rep("", length(molDataIds))
	  
	  tmp <- data.frame(availableTypes=exptype2,availableIds=molDataIds, idNames=molDataNames, properties=moaNames,
	                    stringsAsFactors=FALSE)
	  
	  # Join data.frames
	  results <- rbind(results, tmp)
	  
	  # Reverse Order/ no need for reverse for me
	  #results <- results[rev(rownames(results)),]
	  
	  colnames(results) <- c("Data type","ID ", "Drug Name", "Drug MOA")
	  selsource=metaConfig[[input$dataSrc]][["fullName"]]
	  DT::datatable(results, rownames=FALSE, colnames=colnames(results),extensions='Buttons',
	                filter='top', style='bootstrap4', selection = "none",
	                options=list(pageLength = 10,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipBt',buttons = list('copy', 'print', list(extend = 'collection',buttons = list(list(extend='csv',filename='search_id'), list(extend='excel',filename='search_id'), list(extend='pdf',filename='search_id')),text = 'Download')))
	                , caption=htmltools::tags$caption(paste0("Identifier search for ",selsource),style="color:dodgerblue; font-size: 18px")
	  )})
	
		#--------------------------------------------------------------------------------------
	
	
	#----[Render Data Table in 'Compare Patterns' Tab]-------------------------------------
	output$patternComparison <- DT::renderDataTable({
	# 	srcContent <- srcContentReactive()
	# 	
	# 	if (input$patternComparisonSeed == "xPattern"){
	# 		dat <- xData()
	# 		pcDataset <- input$xDataset
	# 	} else{
	# 		dat <- yData()
	# 		pcDataset <- input$yDataset
	# 	}
	# 	selectedLines <- names(dat$data)
	# 
	#   if(input$patternComparisonType == "drug") {
	#     results <- patternComparison(dat$data,
	#     														 srcContent[[pcDataset]][["molPharmData"]][["act"]][, selectedLines])
	#     results$ids <- rownames(results)
	#     results$NAME <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "NAME"]
	# 
	#     if ("MOA" %in% colnames(srcContent[[pcDataset]][["drugInfo"]])){
	#     	results$MOA <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "MOA"]
	#     	results <- results[, c("ids", "NAME", "MOA", "COR", "PVAL")]
	#     	colnames(results) <- c("ID", "Name", "MOA", "Correlation", "P-Value")
	#     } else{
	#     	results <- results[, c("ids", "NAME", "COR", "PVAL")]
	#     	colnames(results) <- c("ID", "Name", "Correlation", "P-Value")
	#     }
	#     results$FDR=p.adjust(results[,"P-Value"],method="BH",nrow(results))
	#     
	#   } else {
	#     molPharmData <- srcContent[[pcDataset]][["molPharmData"]]
	#     molData <- molPharmData[setdiff(names(molPharmData), c("act","copA","mutA","metA","expA","xaiA","proA","mirA","mdaA","swaA","xsqA"))]
	#     molData <- lapply(molData, function(X) X[, selectedLines])
	#     results <- patternComparison(dat$data, molData)
	#     results$ids <- rownames(results)
	# 
	#     results$molDataType <- getMolDataType(results$ids)
	#     results$gene <- removeMolDataType(results$ids)
	# 
	#     # Reorder columns
	#     results <- results[, c("ids", "molDataType", "gene", "COR", "PVAL")]
	#     colnames(results) <- c("ID", "Data Type", "Gene", "Correlation", "P-Value")
	# 
	#     if (require(rcellminerUtilsCDB)){
	#     	chromLocs <- character(nrow(results))
	#     	haveLoc <- results$Gene %in% names(geneToChromBand)
	#     	chromLocs[haveLoc] <- geneToChromBand[results$Gene[haveLoc]]
	# 
	#     	results$Location <- chromLocs
	#     	results <- results[, c("ID", "Data Type", "Gene", "Location", "Correlation", "P-Value")]
	#     }
	#     results$FDR=p.adjust(results[,"P-Value"],method="BH",nrow(results))
	#     
	#     if (require(geneSetPathwayAnalysis)){
	#     	results$Annotation <- geneSetPathwayAnalysis::geneAnnotTab[results$Gene, "SHORT_ANNOT"]
	# 			results$Annotation[is.na(results$Annotation)] <- ""
	#     }
	#     
	#     results$ID <- NULL
	#   }
	#   
	#   results[, "Correlation"] <- round(results[, "Correlation"], 3)
	#   results[, "P-Value"] <- signif(results[, "P-Value"], 3)
	#   results[, "FDR"] <- signif(results[, "FDR"], 3)
	# 	## sort by p-value
	#   results <- results[order(results[, "P-Value"]),]
    results= PatternCompTable()	
    
    
	  # DT::datatable(results, rownames=FALSE, colnames=colnames(results),extensions='Buttons',
	  # 							filter='top', style='bootstrap4', selection = "none",
	  # 							options=list(lengthMenu = c(10, 50, 100,500), pageLength = 100,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipBt', buttons = list('copy', 'print', list(extend = 'collection',buttons = list(list(extend='csv',filename='pattern_comp',title='Exported data from CellMinerCDB'), list(extend='excel',filename='pattern_comp',title='Exported data from CellMinerCDB'), list(extend='pdf',filename='pattern_comp',title='Exported data from CellMinerCDB')),text = 'Download'))))
	  DT::datatable(results, rownames=FALSE, colnames=colnames(results),
	                filter='top', style='bootstrap4', selection = "none",
	                options=list(lengthMenu = c(10, 50, 100,500), pageLength = 100,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipt'))
	  
	})
	

	#----[Render Data Table in 'Metadata' Tab]-------------------------------------------
	output$cellLineTable <- DT::renderDataTable({
		
		configSelect <- metaConfig[[input$mdataSource]][["packages"]][[1]][["MetaData"]]
		jsonFrame <- as.data.frame(configSelect)
		
		colnames(jsonFrame) <- c("Data Type", "Description", "Units", 
														 "Platform/Assay", "PubMed Ref. ID")
		
		DT::datatable(jsonFrame, rownames=FALSE, colnames=colnames(jsonFrame),
									filter='top', style='bootstrap4', selection = "none",
									options=list(pageLength = 10,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page'))),escape=F)
	})
	#---------------------------------------------------------------------------------------
	output$log <- renderText({
			paste(names(input), collapse=" ")
			query <- parseQueryString(session$clientData$url_search)
			sapply(names(input), function(x) { query[[x]] <- input[[x]] })
		})

# 	output$genUrl <- renderText({
# 		#query <- parseQueryString(session$clientData$url_search)
# 		query <- list()
#
# 		# Get current input values and combine with
# 		tmp <- sapply(names(input), function(x) { query[[x]] <- input[[x]] })
# 		query <- c(query, tmp)
#
# 		paramStr <- paste(sapply(names(query), function(x) { paste0(x, "=", query[[x]]) }), collapse="&")
#
# 		urlStr <- paste(session$clientData$url_protocol, "//",
# 										session$clientData$url_hostname, ":",
# 										session$clientData$url_port,
# 										session$clientData$url_pathname,
# 										"?", paramStr,
# 										sep="")
#
# 		paste0(a(paste("Shareable Link (Right-click then 'Copy' to share)"), href=urlStr), hr())
#
# 		#paste("Shareable Link:", urlStr)
# 	})

  #**********************************************************************************************
	output$tabsetPanel = renderUI({
		#verbatimTextOutput("log") can be used for debugging
		#tabPanel("Plot", verbatimTextOutput("genUrl"), showOutput("rCharts", "highcharts")),
	  tab4 <- tabPanel("Tissue Correlation", value=4,
	                   #downloadLink("downloadData", "Download selected x and y axis data as a Tab-Delimited File"),
	                   DT::dataTableOutput("cortable"))
		tab1 <- tabPanel("View Data", value=2,
                     downloadLink("downloadData", "Download selected x and y axis data as a Tab-Delimited File"),
                     DT::dataTableOutput("table"))
		tab2 <- tabPanel("Search IDs",
                     includeMarkdown("www/files/help.md"),
                     DT::dataTableOutput("ids"))
		tab3 <- tabPanel("Compare Patterns", value=3,
										 includeMarkdown("www/files/help.md"),
										 #br(),
										 HTML("<b>Pattern comparison results are computed with respect to that data defined and shared by both the x and y-axis inputs.</b>"),
										 br(),br(),
										 fluidRow(
                     	#column(3, selectInput("patternComparisonType", "Pattern Comparison",
                      #           						choices=c("Molecular Data"="molData", "Drug Data"="drug"), 
                     	#											selected="molData")),
                     	
                     	column(4, HTML(
                     	  paste("<label class='control-label' for='patternComparisonType'>Select molecular or activity data</label>","<select id='patternComparisonType'><option value='moldata' selected>Molecular Data</option><option value='drug'>Drug Data</option></select>")
                     	)),
										 	#column(3, selectInput("patternComparisonSeed", "With Respect to",
										 	#											choices=c("x-Axis Entry"="xPattern", 
										 	#																"y-Axis Entry"="yPattern"), 
										 	#											selected="xPattern"))
# removed
										 	# column(3, HTML(
										 	#   paste("<label class='control-label' for='patternComparisonSeed' id='lpcs'>With Respect to</label>","<select id='patternComparisonSeed'><option value='xPattern' selected>x-Axis Entry</option><option value='yPattern'>y-Axis Entry</option></select>")
										 	# )),
										 	## new staff
										 	# column(3, HTML(
										 	#   paste("<label class='control-label' for='crossdb'>Use y-Axis cell line set?</label>","<br><input type='radio' id='crossdb' value='No' checked> No  <input type='radio' id='crossdb' value='Yes'> Yes")
										 	# ))
										 	# column(3, radioButtons("crossdb", label = "Select type of comparison", choices = list("Compare x-Axis input to x-Axis molecular or activity data" = "No", "Compare x-Axis input to y-Axis molecular or activity data" = "Yes"), selected  = "No", inline=F)
										 	column(8, radioButtons("crossdb", label = NULL, choices = list("Compare x-Axis input to x-Axis molecular or activity data" = "No", "Compare x-Axis input to y-Axis molecular or activity data" = "Yes"), selected  = "No", inline=F, width="100%")       
										 	 )
										 	
										 ),
										 br(),
										 renderUI({
										   req(PatternCompTable())
										   downloadLink("downloadDataComp", "Download All as a Tab-Delimited File")
										 }),
										 ##downloadLink("downloadDataComp", "Download All as a Tab-Delimited File"),
                     withSpinner(DT::dataTableOutput("patternComparison")))
		
		                 # withSpinner(DT::dataTableOutput("patternComparison")),
		                 # downloadLink("downloadDataComp", "Download All as a Tab-Delimited File"))	

#if(input$hasRCharts == "TRUE") {
#		if (FALSE) {
#			tsPanel <- tabsetPanel(type="tabs",
#									#tabPanel("Plot Data", htmlOutput("genUrl"), showOutput("rCharts", "highcharts")),
#									tabPanel("Plot Data", showOutput("rCharts", "highcharts")),
#									tab1, tab2, tab3
#			)
#		} else {
			# plotPanel <- tabPanel("Plot Data", value=1, plotlyOutput("rChartsAlternative", width = plotWidth, height = plotHeight),
			# 											br(), br(), p("Plot point tooltips provide additional information."))
			plotPanel <- tabPanel("Plot Data", value=1, uiOutput("showCellsUi"), plotlyOutput("rChartsAlternative", width = plotWidth, height = plotHeight),
			                      br(), br(), p("Plot point tooltips provide additional information."))
			#tsPanel <- tabsetPanel(plotPanel, tab1, tab2, tab3)
			### tsPanel <- tabsetPanel(id="ts",plotPanel, tab1, tab3,tab4)
			tsPanel <- tabsetPanel(id="ts",plotPanel, tab1, tab3)
#		}

		return(tsPanel)
	})
	#**********************************************************************************************
	output$tabsetSurvival = renderUI({
	  
	  stab3 <- tabPanel("Heatmap", value=3,
	                   downloadLink("downloadData", "Download selected x and y axis data as a Tab-Delimited File"),
	                   DT::dataTableOutput("table"))
	  stab2 <- tabPanel("Cox Details",value=2,
	                   #includeMarkdown("www/files/help.md"),
	                   withSpinner(verbatimTextOutput("CoxDetails4") ))
	  stab1 <- tabPanel("Survival Plot", value=1,
	                   ## HTML("The survival plots are based on Cox proportional-hazards model. If you enter a set of genes, the survival is based on a risk score which is a weighted average of the genes across all patient samples. The weigth are the beta values of a Cox multivariate analysis using the selected genes as co-variates."),
	                   HTML("The survival plots are based on Cox proportional-hazards model. Multivariate analysis is only allowed for gene expression. If you enter a set of genes, the survival plot will be based on the risk score of the genes across all patient samples."),
	                   br(),br(),
	                    withSpinner(plotOutput("Survplot")))
	                   
	  stab1new <- tabPanel("Survival Plot", value=1,
	                    ## HTML("The survival plots are based on Cox proportional-hazards model. If you enter a set of genes, the survival is based on a risk score which is a weighted average of the genes across all patient samples. The weigth are the beta values of a Cox multivariate analysis using the selected genes as co-variates."),
	                    HTML("The survival plots are based on Cox proportional-hazards model. Multivariate analysis is only allowed for gene expression. If you enter a set of genes, the survival plot will be based on the average of the genes across all patient samples."),
	                    br(),br(),
	                    withSpinner(plotOutput("SurvplotNew")))
	  
	  stab1React <- tabPanel("Kaplan Meier Plot", value=1,
	                       ## old HTML("The survival plots are based on Cox proportional-hazards model. If you enter a set of genes, the survival is based on a risk score which is a weighted average of the genes across all patient samples. The weigth are the beta values of a Cox multivariate analysis using the selected genes as co-variates."),
	                       # HTML("The survival plots are based on Cox proportional-hazards model. Multivariate analysis is only allowed for gene expression. If you enter a set of genes, the survival plot will be based on the average of the genes across all patient samples."),
	                       ### HTML("The plot is based on Cox proportional-hazards model."), br(),
	                       ### HTML(paste0("<b>","Note:","</b>", " the  Multivariate analysis is only allowed for gene expression (not for mutations, metada/signatures). If you enter a set of genes, the survival plot will be based on the average of the genes across all samples (see 'Results Summary' TAB).")),
	                       HTML(paste0("<b>","Note:","</b>", " Multiple genes is only allowed for gene expression but we take the average to display the Kaplan Meier plot. Please look at Cox results TAB for multivariate analysis")),
	                       	                       br(),br(),
	                      ##  withSpinner(plotOutput("SurvplotReact")))
	                        withSpinner(plotOutput("SurvplotReact_km")))
	  stab2React <- tabPanel("Cox results",value=2,
	                    #includeMarkdown("www/files/help.md"),
	                    HTML(paste0("<b>","Note:","</b>", " Survival analysis is based on Cox proportional-hazards model. The multivariate analysis is only allowed for gene expression (not for mutations, metada/signatures).")),
	                    
	                    br(),br(),
	                    
	                    withSpinner(DT::dataTableOutput("CoxDetailsReact") ))

	  stab3React <- tabPanel("Cox Details",value=2,
	                    #includeMarkdown("www/files/help.md"),
	                    withSpinner(verbatimTextOutput("CoxPlusDetailsReact") ))
	  ## stsPanel <- tabsetPanel(id="sts",stab1, stab2)
	  
	  stabHeatReact <- tabPanel("Heatmap",value=2,
	                         #includeMarkdown("www/files/help.md"),
	                         HTML(paste0("<b>","Note:","</b>", " All rows are scaled to have range between 0 and 1.")),
	                         
	                         br(),br(),
	                         withSpinner(plotlyOutput("CoxHeatmapReact") ))
	  
	  stsPanel <- tabsetPanel(id="sts",stab1React, stab2React,stabHeatReact,  stab3React)
	  #		}
	  
	  return(stsPanel)
	})
	
	#**********************************************************************************************
	output$metadataPanel = renderUI({
		#verbatimTextOutput("log") can be used for debugging
		#tabPanel("Plot", verbatimTextOutput("genUrl"), showOutput("rCharts", "highcharts")),
		
		#mtab1 <- tabPanel("Features",
											#DT::dataTableOutput("featTable"))
		mtab2 <- tabPanel("Data Information",
		                  HTML("Here below are the data types that you can download. Please use left panel to proceed."),
		                  br(),br(),
											DT::dataTableOutput("cellLineTable"))
		#mtab3 <- tabPanel("Drug Information",
											#includeMarkdown("www/files/help.md"),
											#DT::dataTableOutput("drugTable"))
		
		dataFullname <- metaConfig[[input$mdataSource]][["fullName"]]
		
		# tabsetPanel(type="pills",
		# 						tabPanel(dataFullname, tags$hr(),
		# 											mtab2)
								tabsetPanel(
								            tabPanel(mtab2)
		)
	})
	##*********************************************************
	# output$searchPanel = renderUI({
	#   #verbatimTextOutput("log") can be used for debugging
	#   #tabPanel("Plot", verbatimTextOutput("genUrl"), showOutput("rCharts", "highcharts")),
	#   
	#   includeMarkdown("www/files/help.md")
	#   DT::dataTableOutput("ids2")
	#  
	# })
	
	output$searchPanel = renderUI({
	  #verbatimTextOutput("log") can be used for debugging
	  #tabPanel("Plot", verbatimTextOutput("genUrl"), showOutput("rCharts", "highcharts")),
	  
	  #includeMarkdown("www/files/help.md")
	  DT::dataTableOutput("ids_s")
	  
	})
	##*********************************************************
	#**************************************************************************************
	output$sourceLink <- renderUI({
		
		urlString <- metaConfig[[input$mdataSource]][["url"]]
		sourceName <- metaConfig[[input$mdataSource]][["displayName"]]
		visibleText <- paste("Select here to learn more about ", sourceName, sep="")
		if (input$mdataSource=="nci60")
		tags$div(
		tags$a(visibleText, href=paste(urlString), target = "_blank"),
		tags$a("   and the DTP",href='https://dtp.cancer.gov', target='_blank'))
		else tags$a(visibleText, href=paste(urlString), target = "_blank",id="tmd", class ="dm")
		
		#  
		# if (input$mdataSource=="nci60") { 
		#    a(visibleText, href=paste(urlString), target = "_blank")
		#    a("DTP",href='https://dtp.cancer.gov', target='_blank')
		# }
	})
	#**********************************************************************************************

  output$downloadData <- downloadHandler(
    filename = function() {
      query <- parseQueryString(session$clientData$url_search)
      
      if("filename" %in% names(query)) {
        filename <- query[["filename"]]
      } else {
        filename <- "dataset"
      }
      
      if("extension" %in% names(query)) {
        extension <- query[["extension"]]
      } else {
        extension <- "txt"
      }
      
      paste(filename, extension, sep=".")
    },
    content = function(file) {
      df <- getPlotData(xData = xData(), yData = yData(), showColor = input$showColor, 
                        showColorTissues = showColorTissues(), dataSource = input$xDataset, 
                        srcContent = srcContentReactive())
      
      # Column selection below is to restrict to cell line, x, y features,
      # and tissue type information (source-provided + OncoTree).
      
      dfCols <- c("Cell Line","patientID","paired_sample",colnames(df)[2:3], paste0("OncoTree", 1:2))
      ### ----------
      if ("EMT" %in% colnames(df)) {
        df[,"EMT"] <- gsub("EMT:","",df[,"EMT"])
        dfCols <- c(dfCols, "EMT")
      }
      
      if ("priorTreatment" %in% colnames(df)) {
        df[,"priorTreatment"] <- gsub("PriorTreatment:","",df[,"priorTreatment"])
        dfCols <- c(dfCols, "priorTreatment")
      }
      
      
      # new stuff
      if ("sampleType" %in% colnames(df)) {
        df[,"sampleType"] <- gsub("SampleType:","",df[,"sampleType"])
        dfCols <- c(dfCols, "sampleType")
      }
      if ("biopsySite" %in% colnames(df)) {
        df[,"biopsySite"] <- gsub("BiopsySite:","",df[,"biopsySite"])
        dfCols <- c(dfCols, "biopsySite")
      }
      # if ("molecularSubtype" %in% colnames(df)) {
      #   df[,"molecularSubtype"] <- gsub("MolecularSubtype:","",df[,"molecularSubtype"])
      #   dfCols <- c(dfCols, "molecularSubtype")
      # }
      if ("dataSet" %in% colnames(df)) {
        df[,"dataSet"] <- gsub("DataSet:","",df[,"dataSet"])
        dfCols <- c(dfCols, "dataSet")
      }
      
      df <- df[, dfCols]
      # df[, 2] <- round(df[, 2], 3)
      # df[, 3] <- round(df[, 3], 3)
      colnames(df)[1]="Patient_sample"
      
      ### ----------
      # dfCols <- c(colnames(df)[1:4], paste0("OncoTree", 1:4))
      # if ("EMT" %in% colnames(df)) {
      #   dfCols <- c(dfCols, "EMT")
      # }
      # df <- df[, dfCols]
      
      write.table(df, file, quote=FALSE, row.names=FALSE, sep="\t")
    }
  )
  
  ### onco3 for survival -------------------------------------------------------
  output$onco3Ui <- renderUI({
    srcContent <- srcContentReactive()
    
    # The last selected (data type) prefix is recorded in 
    # globalReactiveValues$xPrefix whenever xData() is updated. When the data set 
    # is changed, we try to use this same data type prefix, if it is available.
    prefixChoices <- srcContent[[input$survSource]][["featurePrefixes"]]
    selectedPrefix <- globalReactiveValues$xPrefix
    if ((is.null(selectedPrefix)) || (!(selectedPrefix %in% prefixChoices))){
      selectedPrefix <- srcContent[[input$mdataSource]][["defaultFeatureX"]]
      if (is.na(selectedPrefix)) selectedPrefix <- srcContent[[input$mdataSource]][["defaultFeatureY"]]
    }
    
    ##onco3choices = na.omit(unique(srcContent[[input$survSource]][["sampleData"]][["OncoTree3"]]))
    onco3choices = sort(na.omit(unique(srcContent[[input$survSource]][["sampleData"]][["dataSource"]])))
    
    ## new
    onco3choices = setdiff(onco3choices, OutOfSurvival) 
    
    opt = "";
    for(y in 1:length(onco3choices)){
     if  (!(onco3choices[y]  %in% c("GInorm_IOWA","PCnorm_TCGA", "Pnorm_IOWA","Tnorm_NIDDK")))  { 
      if (onco3choices[y]=="SCLC_TU")
      {
        opt =  paste0(opt,"<option value=",onco3choices[y]," selected>",onco3choices[y],"</option>;")
      }
      else
      {
        opt =  paste0(opt,"<option value=",onco3choices[y],">",onco3choices[y],"</option>;");
      }
     }
    }
    # selectInput("xPrefix", "x-Axis Type", choices = prefixChoices, selected = selectedPrefix)
    HTML(
      paste("<label class='control-label' for='onco3ch'>Select Data set</label>","<select id='onco3ch' size='4' style='word-wrap:break-word; width: 100%;' multiple>",opt,"</select>")
    )
  })
  
  ### data types for search ---------
  output$dataTypeUi <- renderUI({
    srcContent <- srcContentReactive()
    
    # The last selected (data type) prefix is recorded in 
    # globalReactiveValues$xPrefix whenever xData() is updated. When the data set 
    # is changed, we try to use this same data type prefix, if it is available.
    
    prefixChoices <- srcContent[[input$mdataSource]][["featurePrefixes"]]
    ## prefixChoices <- srcContent[["patient"]][["featurePrefixes"]] ## Dec 20, 2024
    ii = which(prefixChoices=="xsq")
    names(prefixChoices)[ii] = paste(names(prefixChoices)[ii], "with batch correction")
    ##
    selectedPrefix <- globalReactiveValues$xPrefix
    if ((is.null(selectedPrefix)) || (!(selectedPrefix %in% prefixChoices))){
      selectedPrefix <- srcContent[[input$mdataSource]][["defaultFeatureX"]]
      if (is.na(selectedPrefix)) selectedPrefix <- srcContent[[input$mdataSource]][["defaultFeatureY"]]
    }
    opt = "";
    for(y in 1:length(prefixChoices)){
      if (prefixChoices[y]==selectedPrefix)
      {
        opt =  paste0(opt,"<option value=",prefixChoices[y]," selected>",names(prefixChoices)[y],"</option>;")
      }
      else
      {
        opt =  paste0(opt,"<option value=",prefixChoices[y],">",names(prefixChoices)[y],"</option>;");
      }
    }
    # selectInput("xPrefix", "x-Axis Type", choices = prefixChoices, selected = selectedPrefix)
    HTML(
      # paste("<label class='control-label' for='dataType'>Select Data Type to Download</label>","<select id='dataType' style='word-wrap:break-word; width: 100%;'>",opt,"</select>")
      paste("<label class='control-label' for='dataType'>Please select Data Type and click on Download button</label>","<br><select id='dataType' style='word-wrap:break-word; width: 25%;'>",opt,"</select>")
    )
  })
 ## new version for search tab
  output$dataTypeUi_s <- renderUI({
    srcContent <- srcContentReactive()
    
    # The last selected (data type) prefix is recorded in 
    # globalReactiveValues$xPrefix whenever xData() is updated. When the data set 
    # is changed, we try to use this same data type prefix, if it is available.
    prefixChoices <- srcContent[[input$dataSrc]][["featurePrefixes"]]
    ## new
    ivar = match("var",prefixChoices)
    if (!is.na(ivar)) 	prefixChoices = prefixChoices[-ivar]
    ##
    selectedPrefix <- globalReactiveValues$xPrefix
    if ((is.null(selectedPrefix)) || (!(selectedPrefix %in% prefixChoices))){
      selectedPrefix <- srcContent[[input$dataSrc]][["defaultFeatureX"]]
      if (is.na(selectedPrefix)) selectedPrefix <- srcContent[[input$dataSrc]][["defaultFeatureY"]]
    }
    opt = "";
    for(y in 1:length(prefixChoices)){
      if (prefixChoices[y]==selectedPrefix)
      {
        opt =  paste0(opt,"<option value=",prefixChoices[y]," selected>",names(prefixChoices)[y],"</option>;")
      }
      else
      {
        opt =  paste0(opt,"<option value=",prefixChoices[y],">",names(prefixChoices)[y],"</option>;");
      }
    }
    # selectInput("xPrefix", "x-Axis Type", choices = prefixChoices, selected = selectedPrefix)
    HTML(
      paste("<label class='control-label' for='dataTyp'>Select Data Type</label>","<select id='dataTyp' style='word-wrap:break-word; width: 100%;'>",opt,"</select>")
    )
  })
  
  
  
  ### Download data
  output$downloadExp_orig <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      ## paste0("data_",input$mdataSource,"_",input$dataType,".txt")
      # paste0("data_",metaConfig[[input$mdataSource]][["displayName"]],"_",input$dataType,".txt")
      paste0("data_",metaConfig[[input$mdataSource]][["displayName"]],"_",input$dataType,".zip")
    },
    # filename = function() {
    #   paste0(input$mdataSource,"_",input$dataType,".zip")
    # },
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      file0 = paste0(tempdir(),"/data_",metaConfig[[input$mdataSource]][["displayName"]],"_",input$dataType,".txt")
      wdata=srcContent[[input$mdataSource]][["molPharmData"]][[input$dataType]]
      # new cancel next
      # rownames(wdata)=substr(rownames(wdata),4,nchar(rownames(wdata)))

      # Write to a file specified by the 'file' argument
      #write.table(srcContent[[input$mdataSource]][["molPharmData"]][[input$dataType]], file, sep = "\t",col.names = NA)      
      
      #new -------------------------------------------
      if (input$dataType!="act") {
      wdata.A=srcContent[[input$mdataSource]][["molPharmData"]][[paste0(input$dataType,"A")]]
          gene = substr(rownames(wdata),4,nchar(rownames(wdata)))
       if (nrow(wdata.A)==0) {
         
         # final = cbind(Gene=substr(rownames(wdata),4,nchar(rownames(wdata))),wdata)
         final = cbind(Gene=gene,wdata)
       }
       else{
          # stopifnot(identical(rownames(wdata),rownames(wdata.A)))
          stopifnot(identical(gene,as.character(wdata.A[,1])))
          final = cbind(wdata.A,wdata)
          ### rownames(final)=rownames(wdata)
       }
      }
      else {
        wdata.A=srcContent[[input$mdataSource]][["drugInfo"]]
        stopifnot(identical(rownames(wdata),rownames(wdata.A)))
        final = cbind(wdata.A,wdata)
      }
      ## write.table(final, file, sep = "\t", row.names = F)  
      write.table(final, file0, sep = "\t", row.names = F)  
      zip(zipfile=file, files=file0, flags = "-r9Xj")
      # ------------------------------------------------
      
      # new cancel next
      # write.table(wdata, file, sep = "\t", col.names = NA)  
      
      # myname=paste0(input$mdataSource,"_",input$dataType,".txt")
      # write.table(srcContent[[input$mdataSource]][["molPharmData"]][[input$dataType]], myname, sep = "\t",
      #             col.names = NA)
      # configSelect <- metaConfig[[input$mdataSource]][["packages"]][[1]][["MetaData"]]
      # jsonFrame <- as.data.frame(configSelect)
      # 
      # colnames(jsonFrame) <- c("DataType", "Description", "Units", 
      #                          "Platform/Assay", "PubMed Ref. ID")
      # write.table(t(jsonFrame[which(jsonFrame$DataType==input$dataType),]),"footnotes.txt",sep="\t",col.names=F)
      # zip(zipfile=file, files=c(myname,"footnotes.txt"))
    },
    contentType = "application/zip"
  )
##

  ### Download data
  output$downloadExp <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      ## paste0("data_",input$mdataSource,"_",input$dataType,".txt")
      # paste0("data_",metaConfig[[input$mdataSource]][["displayName"]],"_",input$dataType,".txt")
      paste0("sclc_",metaConfig[[input$mdataSource]][["displayName"]],"_",input$dataType,".zip")
    },
    content = function(file) {
      myfile = paste0(downloadpath,"/sclc_",metaConfig[[input$mdataSource]][["displayName"]],"_",input$dataType,".zip")
       
      # if (!file.exists(myfile)) showModal(modalDialog(title="Error", p("file does not exists")))
      #  file.copy(myfile, file)  
      
      if (file.exists(myfile)) file.copy(myfile, file)  
      else {
        file0 = paste0(tempdir(),"/sclc_",metaConfig[[input$mdataSource]][["displayName"]],"_",input$dataType,".txt")
        wdata=srcContent[[input$mdataSource]][["molPharmData"]][[input$dataType]]

        #new -------------------------------------------
        if (input$dataType!="act") {
          wdata.A=srcContent[[input$mdataSource]][["molPharmData"]][[paste0(input$dataType,"A")]]
          gene = substr(rownames(wdata),4,nchar(rownames(wdata)))
          if (nrow(wdata.A)==0) {
            
            # final = cbind(Gene=substr(rownames(wdata),4,nchar(rownames(wdata))),wdata)
            final = cbind(ID=gene,wdata)
          }
          else{
            # stopifnot(identical(rownames(wdata),rownames(wdata.A)))
            if (input$mdataSource == "nciSclc" & input$dataType =="exp" ) {
              stopifnot(identical(gene,as.character(wdata.A[,4]))) 
            } else 
            {
               stopifnot(identical(gene,as.character(wdata.A[,1]))) 
            }
            
            final = cbind(wdata.A,wdata)
          }
        }
        else {
          wdata.A=srcContent[[input$mdataSource]][["drugInfo"]]
          stopifnot(identical(rownames(wdata),rownames(wdata.A)))
          final = cbind(wdata.A,wdata)
        }
        ## write.table(final, file, sep = "\t", row.names = F)  
        write.table(final, file0, sep = "\t", row.names = F)  
        zip(zipfile=file, files=file0, flags = "-r9Xj")
        
      }
      
    },
    contentType = "application/zip"
  )
  ## ------------------------------
  output$downloadExpNOBER <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste0("sclc_Patient data_xsq_noBer",".zip")
    },
    content = function(file) {
      myfile = paste0(downloadpath,"/sclc_Patient data_xsq_noBer.zip")
      shiny::validate(need(length(analysisTissueTypes()) > 0, 
                           "There are no cell lines of the selected tissue type(s)."))
      ## shiny::validate(need(file.exists(myfile), "The file does not exist. Please contact the website administrator"))
      if (file.exists(myfile)) file.copy(myfile, file)  
       else
         shinyalert::shinyalert(
           # title = "Warning in Cox 2 groups",
           # text = paste("Message: ",w),
           title = "File not available.",
           text = " Please contact the website administrator",
           type = "error",
           imageHeight = 150,
           size ="m"
         ) 
      ## file.copy(myfile, file)  
    },
    contentType = "application/zip"
  )
  
  
  
## -----------Download Cell line info
  output$downloadCell <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste0("Sample_annotation_",input$mdataSource,".txt")
    },
     content = function(file) {
      
      wdata=srcContent[[input$mdataSource]][["sampleData"]]
      # rownames(wdata)=substr(rownames(wdata),4,nchar(rownames(wdata)))
      # Write to a file specified by the 'file' argument
      #write.table(srcContent[[input$mdataSource]][["molPharmData"]][[input$dataType]], file, sep = "\t",col.names = NA)      
      write.table(wdata, file, sep = "\t", row.names = F)  
      
    }
  )
##
  
### Download footnotes
  output$downloadFoot <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      # paste0("footnotes_",input$mdataSource,"_",input$dataType,".csv")
      paste0("footnotes_",metaConfig[[input$mdataSource]][["displayName"]],"_",input$dataType,".csv")
      
    },
    # filename = function() {
    #   paste0(input$mdataSource,"_",input$dataType,".zip")
    # },
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      
      # Write to a file specified by the 'file' argument
      #write.table(srcContent[[input$mdataSource]][["molPharmData"]][[input$dataType]], file, sep = "\t",
      #            col.names = NA)      
      
      # myname=paste0(input$mdataSource,"_",input$dataType,".txt")
      # write.table(srcContent[[input$mdataSource]][["molPharmData"]][[input$dataType]], myname, sep = "\t",
      #             col.names = NA)
      configSelect <- metaConfig[[input$mdataSource]][["packages"]][[1]][["MetaData"]]
      jsonFrame <- as.data.frame(configSelect)
       
      colnames(jsonFrame) <- c("DataType", "Description", "Units", 
                                "Platform/Assay", "PubMed Ref. ID(s)")
      mydata=t(jsonFrame[which(jsonFrame$DataType==input$dataType),])
      ##
      nr=nrow(mydata)
      mydata = rbind(" ",metaConfig[[input$mdataSource]][["displayName"]],mydata," ","For the 'Download Data' table for this data source, row one contains the sample names,","column one contains the 'Data Type' identifier, and the numerical values are ","the 'Units' for the 'Platform/Assay', as described in the 'PubMed Reference' above.","More details are in the help section.")
      rownames(mydata)[1:2]=c("Footnotes:","Data Set")
      rownames(mydata)[nr+4]="Remarks"
      ## colnames(mydata)=paste("footnotes for data source:",input$mdataSource)
       write.table(mydata,file,col.names = F,sep=",")
       #write.table(t(jsonFrame[which(jsonFrame$DataType==input$dataType),]),file,sep="\t",col.names=F)
       # zip(zipfile=file, files=c(myname,"footnotes.txt"))
    }
  )
  
  ### Download Synonym table
  output$downloadSyn <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste0("Table_","Drugs_Synonyms_cdb",".txt")
    },
    
     content = function(file) {
      
      write.table(findDrugIDs("*"), file, sep = "\t", row.names = F,quote=F)  
      
     }
  )
  ##
  output$downloadDataComp <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
#      if (input$patternComparisonSeed == "xPattern"){
        
        pcDataset <- input$xDataset
        pcType <- input$xPrefix
        pcId <- input$xId
      # } 
      # else{
      #   
      #   pcDataset <- input$yDataset
      #   pcType <- input$yPrefix
      #   pcId <- input$yId
      # 
      # }
        pcDataset2 <- input$yDataset
      if (input$crossdb == "No") {  
      ##paste0("Pattern_Comp_all_cdb_",input$patternComparisonSeed,"_",pcDataset,"_",pcType,"_",pcId,"_",input$patternComparisonType,".txt")
      paste0("Pattern_Comp_all_cdb_",pcDataset,"_",pcType,"_",pcId,"_",input$patternComparisonType,".txt") 
      }
        else {
          paste0("Pattern_Comp_all_cdb_",pcDataset,"_",pcType,"_",pcId,"_",input$patternComparisonType,"_from_",pcDataset2,".txt") 
        }
      },
    
    content = function(file) {
      
      write.table(PatternCompTable(), file, sep = "\t", row.names = F,quote=F)  
      
    }
  )
  
  
  ###
  
  output$xPrefixUi <- renderUI({
  	srcContent <- srcContentReactive()
  	
  	# cat("source: ", current_2var(), "\n")
  	# The last selected (data type) prefix is recorded in 
  	# globalReactiveValues$xPrefix whenever xData() is updated. When the data set 
  	# is changed, we try to use this same data type prefix, if it is available.
  	prefixChoices <- srcContent[[input$xDataset]][["featurePrefixes"]]
  	prefixChoices = sort(prefixChoices, decreasing = T) ## Dec20, 2024
  	
  	### new
  	ivar = match("var",prefixChoices)
  	if (!is.na(ivar)) 	prefixChoices = prefixChoices[-ivar]
  	###
  	selectedPrefix <- globalReactiveValues$xPrefix
  	if ((is.null(selectedPrefix)) || (!(selectedPrefix %in% prefixChoices))){
  		selectedPrefix <- srcContent[[input$xDataset]][["defaultFeatureX"]]
  		if (is.na(selectedPrefix)) selectedPrefix <- srcContent[[input$xDataset]][["defaultFeatureY"]]
  	}
  	opt = "";
  	for(y in 1:length(prefixChoices)){
  	  if (prefixChoices[y]==selectedPrefix)
  	  {
  	    opt =  paste0(opt,"<option value=",prefixChoices[y]," selected>",names(prefixChoices)[y],"</option>;")
  	  }
  	  else
  	  {
  	    opt =  paste0(opt,"<option value=",prefixChoices[y],">",names(prefixChoices)[y],"</option>;");
  	  }
  	}
  	selectInput("xPrefix", "Data Type", choices = prefixChoices, selected = selectedPrefix)
  ## 	selectInput("xPrefix", "x-Axis Type", choices = prefixChoices, selected = selectedPrefix)
  	# HTML(
  	#   paste("<label class='control-label' for='xPrefix'>x-Axis Data Type</label>","<select id='xPrefix' style='word-wrap:break-word; width: 100%;'>",opt,"</select>")
  	# )
  })

  ## new interface for xId: textInput("xId", "Identifier: (e.g. ASCL1)", "ASCL1") -----------
  
  output$xAdcUi <- renderUI({
   
    if ( length(input$xPrefix)!=0 )  {
      if (input$xPrefix == "xsq")
      {
        ## 11/07 ADC
        checkboxInput("xAdc", "ADC genes ?", value=FALSE)
        ## end ADC
      }
      }

  })
  
  varprevious <- debounce(reactive({
    if (is.null(input$varname)) "SLFN11" else input$varname
  }), 2000) # was 5000 
  
  output$varUi <- renderUI({
    
    if ( length(input$optsurv)!=0 )  {
      
        if (input$optsurv == "mda") {
          cellChoices2 <- sort(srcContent[[input$xDataset]][["molPharmData"]][["mdaA"]]$ID)
          selectInput("varname", "Please select a miscellaneous variable",choices = cellChoices2, selected = "TMB")
        } else 
        { 
          if (input$optsurv== "sig") { 
            cellChoices3 <- sort(srcContent[[input$xDataset]][["molPharmData"]][["sigA"]]$ID)
            selectInput("varname", "Please select a signature",choices = cellChoices3, selected = "NE_mean")
          }
          else 
          {
            if (input$optsurv == "nmf") { 
              cellChoices3 <- sort(srcContent[[input$xDataset]][["molPharmData"]][["nmfA"]]$ID)
              selectInput("varname", "Please select an NMF cluster",choices = cellChoices3, selected = "NMF1u")
            }
            else
            {
              # textInput("xId", "Identifier: (e.g. ASCL1)", xprevious()) 
              textInput("varname", "Gene(s) Symbol(s) separated by space(s) (e.g. SLFN11)", varprevious()) 
            }
          }
           
        }
      
    }
    else 
      
    {
      textInput("varname", "Gene(s) Symbol(s) separated by space(s) (e.g. SLFN11)", "SLFN11")
    }
    
  })
  
  xprevious <- debounce(reactive({
    if (is.null(input$xId)) "ASCL1" else input$xId
  }), 5000)
  
   output$xIdUi <- renderUI({
  
    ## textInput("xId", "Identifier: (e.g. ASCL1)", "ASCL1")
    
  if ( length(input$xPrefix)!=0 )  {
      if (input$xPrefix == "hll")
      {
      # srcContent <- srcContentReactive()
        cellChoices <- sort(srcContent[[input$xDataset]][["molPharmData"]][["hllA"]]$ID)
        selectInput("xId", "Please select a pathway",choices = cellChoices, selected = "DNA_REPAIR")
      } else 
      {
        if (input$xPrefix == "mda") {
          cellChoices2 <- sort(srcContent[[input$xDataset]][["molPharmData"]][["mdaA"]]$ID)
          selectInput("xId", "Please select a miscellaneous variable",choices = cellChoices2, selected = "TMB")
        } else 
        { 
          # 11/07 adc
          # if (input$xPrefix == "xsq" & length(input$xAdc) != 0)  {
          #     if (input$xAdc) {
          #       cellChoices3 <- sort(adcgenes$gene)
          #       selectInput("xId", "Please select a gene",choices = cellChoices3, selected = "DLK1")
          #       } else 
          #       {
          #         textInput("xId", "Identifier: (e.g. ASCL1)", xprevious())
          #       }
          #    }
          #   else
          #   { 
          
          if (input$xPrefix == "sig") { 
            cellChoices3 <- sort(srcContent[[input$xDataset]][["molPharmData"]][["sigA"]]$ID)
            selectInput("xId", "Please select a signature",choices = cellChoices3, selected = "NE_mean")
          }
          else 
          {
            if (input$xPrefix == "nmf") { 
              cellChoices3 <- sort(srcContent[[input$xDataset]][["molPharmData"]][["nmfA"]]$ID)
              selectInput("xId", "Please select an NMF cluster",choices = cellChoices3, selected = "NMF1u")
            }
             else
             {
                textInput("xId", "Identifier: (e.g. ASCL1)", xprevious()) 
             }
          }
          #    }
        }
      }
    }
   else 
     
     {
        textInput("xId", "Identifier: (e.g. ASCL1)", "ASCL1")
    }
      
  })
  
   output$yAdcUi <- renderUI({
     
     if ( length(input$yPrefix)!=0 )  {
       if (input$yPrefix == "xsq")
       {
         ## 11/07 ADC
         checkboxInput("yAdc", "ADC genes ?", value=FALSE)
         ## end ADC
       }
     }
     
   })
   
   yprevious <- debounce(reactive({
     if (is.null(input$yId)) "YAP1" else input$yId
   }), 5000)
   
   
  output$yIdUi <- renderUI({
    
    ## textInput("xId", "Identifier: (e.g. ASCL1)", "ASCL1")
    
    if ( length(input$yPrefix)!=0 )  {
      if (input$yPrefix == "hll")
      {
        # srcContent <- srcContentReactive()
        cellChoices <- sort(srcContent[[input$yDataset]][["molPharmData"]][["hllA"]]$ID)
        selectInput("yId", "Please select a pathway",choices = cellChoices, selected = "APOPTOSIS")
      } else 
      {
        if (input$yPrefix == "mda") {
          cellChoices2 <- sort(srcContent[[input$yDataset]][["molPharmData"]][["mdaA"]]$ID)
          selectInput("yId", "Please select a miscellaneous variable",choices = cellChoices2, selected = "Age")
        } else 
        {
          ##--------------------
          # 11/07 adc
        #   if (input$yPrefix == "xsq" & length(input$yAdc) != 0)  {
        #     if (input$yAdc) {
        #       cellChoices3 <- sort(adcgenes$gene)
        #       selectInput("yId", "Please select a gene",choices = cellChoices3, selected = "CD274")
        #     } else 
        #     {
        #       textInput("yId", "Identifier: (e.g. YAP1)", yprevious())
        #     }
        #   }
        # else
        # { 
        if (input$yPrefix == "sig") {
          cellChoices3 <- sort(srcContent[[input$yDataset]][["molPharmData"]][["sigA"]]$ID)
          selectInput("yId", "Please select a signature",choices = cellChoices3, selected = "APM_mean")
        }
          else
          {
            if (input$yPrefix == "nmf") {
              cellChoices3 <- sort(srcContent[[input$yDataset]][["molPharmData"]][["nmfA"]]$ID)
              selectInput("yId", "Please select an NMF cluster",choices = cellChoices3, selected = "NMF1r")
            }
            else {
             textInput("yId", "Identifier: (e.g. YAP1)", yprevious()) 
             }
          }
        # }
      }
    }
    }
  else {
      
      textInput("yId", "Identifier: (e.g. YAP1)", "YAP1")
      
    }
    
  })
  
  ## end new -----------------------------------------------------------------
  
  output$yPrefixUi <- renderUI({
  	srcContent <- srcContentReactive()
  	prefixChoices <- srcContent[[input$yDataset]][["featurePrefixes"]]
  	prefixChoices = sort(prefixChoices, decreasing = T) # Dec 20, 2024
  	### new
  	ivar = match("var",prefixChoices)
  	if (!is.na(ivar)) 	prefixChoices = prefixChoices[-ivar]
  	###
  	selectedPrefix <- globalReactiveValues$yPrefix
  	if ((is.null(selectedPrefix)) || (!(selectedPrefix %in% prefixChoices))){
  		selectedPrefix <- srcContent[[input$yDataset]][["defaultFeatureY"]]
  	}
  	opt = "";
  	for(y in 1:length(prefixChoices)){
  	  if (prefixChoices[y]==selectedPrefix)
  	  {
  	    opt =  paste0(opt,"<option value=",prefixChoices[y]," selected>",names(prefixChoices)[y],"</option>;")
  	  }
  	  else
  	  {
  	    opt =  paste0(opt,"<option value=",prefixChoices[y],">",names(prefixChoices)[y],"</option>;");
  	  }
  	}
  	## selectInput("yPrefix", "y-Axis Type", choices = prefixChoices, selected = selectedPrefix)
  	selectInput("yPrefix", "Data Type", choices = prefixChoices, selected = selectedPrefix)
  	
  	  	# HTML(
  	#   paste("<label class='control-label' for='yPrefix' id='lyp'>y-Axis Data Type</label>","<select id='yPrefix' style='word-wrap:break-word; width: 100%;'>",opt,"</select>")
  	# )
  })
  
  output$xAxisRangeUi <- renderUI({
  	srcContent <- srcContentReactive()
  	
  	# Note: req() ensures values are available or 'truthy' (not NULL, "", FALSE, empty, etc.),
  	# returning the value if so; otherwise the operation is stopped with a silent exception.
  	# The idea is to exit quietly if inputs are momentarily in an invalid state, as might
  	# occur when the app is first loading, etc.
  	
  	## new
  	shiny::validate(need(!is.na(match(input$xPrefix,srcContentReactive()[[input$xDataset]][["featurePrefixes"]])), "Non valid data type"))
  	
  	valRange <- srcContent[[req(input$xDataset)]][["featureValRanges"]][[req(input$xPrefix)]]
  	
  	xData <- NULL
  	try(xData <- xData())
  	if (is.null(xData)){
  		xInitSliderVals <- valRange
  	} else{
  		xDataRange <- range(xData$data, na.rm = TRUE)
  		delta <- max((0.05 * (xDataRange[2] - xDataRange[1])), 0.1)
  		xInitSliderVals <- c((xDataRange[1] - delta), (xDataRange[2] + delta))
  	}
  	
##  	sliderInput("xAxisRange", "x-Axis Range", min = valRange[1], max = valRange[2], value = xInitSliderVals, step = 0.5)
  	
  	sliderInput("xAxisRange", "Range", 
  	            min = valRange[1], max = valRange[2], value = xInitSliderVals, step = 0.5)
  })
  
  output$yAxisRangeUi <- renderUI({
  	srcContent <- srcContentReactive()
  	
 		# Note: see comment in output#xAxisRangeUi explaining the use of req().
  	## new
  	shiny::validate(need(!is.na(match(input$yPrefix,srcContentReactive()[[input$yDataset]][["featurePrefixes"]])), "Non valid data type"))
  	
  	valRange <- srcContent[[req(input$yDataset)]][["featureValRanges"]][[req(input$yPrefix)]]
  	
  	yData <- NULL
  	try(yData <- yData())
  	if (is.null(yData)){
  		yInitSliderVals <- valRange
  	} else{
  		yDataRange <- range(yData$data, na.rm = TRUE)
  		delta <- max((0.05 * (yDataRange[2] - yDataRange[1])), 0.1)
  		yInitSliderVals <- c((yDataRange[1] - delta), (yDataRange[2] + delta))
  	}
  	
  ## 	sliderInput("yAxisRange", "y-Axis Range", min = valRange[1], max = valRange[2], value = yInitSliderVals, step = 0.5)
  	sliderInput("yAxisRange", "Range", 
  	            min = valRange[1], max = valRange[2], value = yInitSliderVals, step = 0.5)
  })

  # output$xIdUi <- renderUI({
  # 	updateSelectizeInput(session, inputId='xId', choices=xIdChoices(), selected="SLFN11", server=TRUE)
  # 	selectizeInput('xId', label="ID: (e.g. 94600 or SLFN11); Case-Sensitive", choices=NULL, options=list(maxOptions=5))
  # })
  # 
  # output$yIdUi <- renderUI({
  # 	updateSelectizeInput(session, inputId='yId', choices=yIdChoices(), selected="94600", server=TRUE)
  # 	selectizeInput('yId', label="ID: (e.g. 94600 or SLFN11); Case-Sensitive", choices=NULL, options=list(maxOptions=5))
  # })
  
  convertListTree <- function(vect)
  {
    mylist = list()
    for (k in 1:length(vect))
    {
      rec = unlist(strsplit(vect[k], ":")) # take all levels, we assume 4 levels max
      stopifnot(length(rec)>=1)
      # first element
      if ( is.null(mylist[[rec[1]]]) ) { 
        # mylist = c(mylist, list())
        mylist = c(mylist, list(list()))
        names(mylist)[length(mylist)] = rec[1]
      }
      mylist
      # 2nd element, label is different LATER
      if (!is.na(rec[2])) { 
        temp = mylist[[rec[1]]]
        if ( is.null(temp[[rec[2]]]) ) { 
          if (length(temp)!=0) temp = c(temp, list(list())) else temp=list(list())
          names(temp)[length(temp)] = rec[2]
          mylist[[rec[1]]] = temp
        }
        
        # 3rd element, label is different LATER
        if (!is.na(rec[3])) { 
          temp3 = temp[[rec[2]]]
          if ( is.null(temp3[[rec[3]]]) ) { 
            if (length(temp3)!=0) temp3 = c(temp3, list(list())) else temp3=list(list())
            # temp3 = c(temp3, list()) 
            names(temp3)[length(temp3)] = rec[3]
            temp[[rec[2]]]= temp3
            mylist[[rec[1]]] = temp
          }
          
          if (!is.na(rec[4])) { 
            temp4 = temp3[[rec[3]]]
            if ( is.null(temp4[[rec[4]]]) ) { 
              if (length(temp4)!=0) temp4 = c(temp4, list(list())) else temp4=list(list())
              # temp4 = c(temp4, list())
              names(temp4)[length(temp4)] = rec[4]
              temp3[[rec[3]]]= temp4
              temp[[rec[2]]]= temp3
              mylist[[rec[1]]] = temp
            }
          }
          
        } #rec3
        
      } # rec2
      
    }
    return(mylist)
  }
  
  output$tree <- renderTree({
    srcContent <- srcContentReactive()
    tissueToSamplesMap <- srcContent[[input$xDataset]][["tissueToSamplesMap"]]
    ## tissueTypes <- sort(unique(names(tissueToSamplesMap)))
    tissueTypes <- unique(names(tissueToSamplesMap)) # no sorting
    ## very new ------------------------------------------------------------------------###
    mylist <- convertListTree(tissueTypes)
    if (input$tissueSelectionMode == "To include"){
      mylist = list(all=structure(mylist, stselected=TRUE))
      
      
    } else{ # input$tissueSelectionMode == "To exclude"
      mylist = list(none=structure(mylist, stselected=TRUE))
      
    }
    mylist
  })
  
  output$selectTissuesUi <- renderUI({
  	# srcContent <- srcContentReactive()
  	# tissueToSamplesMap <- srcContent[[input$xDataset]][["tissueToSamplesMap"]]
  	# tissueTypes <- sort(unique(names(tissueToSamplesMap)))
  	## very new ------------------------------------------------------------------------###
  	
  	# HTML("<p><b>Select Tissue/s of Origin</b></p>")
  	shinyTree("tree", checkbox = F, theme = "proton")
  	## end very new ------------------------------------------------------------------- ###
  	
  	# ## new code
  	# if (input$tissueSelectionMode == "To include"){
  	#        choices=c("all", tissueTypes); mysel="all"
  	# 
  	# } else{ # input$tissueSelectionMode == "To exclude"
  	#        choices=c("none", tissueTypes); mysel="none"
  	# }
  	# opt = "";
  	# for(y in 1:length(choices)){
  	#   # style works only for browser Chrome
  	#   if (choices[y]==mysel)
  	#   {
  	#     opt =  paste0(opt,"<option style='white-space: pre-wrap' selected>",choices[y],"</option>;");
  	#   }
  	#   else {
  	#   opt =  paste0(opt,"<option style='white-space: pre-wrap'>",choices[y],"</option>;");
  	#   }
  	# }
  	# HTML(
  	#   paste("<label class='control-label' for='selectedTissues'>Select Tissue/s of Origin</label>","<select id='selectedTissues' style='word-wrap:break-word; width: 100%;' multiple>",opt,"</select>")
  	# )
  	# # # end new
  	# 
  	# # if (input$tissueSelectionMode == "To include"){
  	# # selectInput("selectedTissues", label = "Select Tissue/s of Origin", choices=c("all", tissueTypes),
  	# # 						multiple=TRUE, selected="all")
  	# # } else { # input$tissueSelectionMode == "To exclude"
  	# # selectInput("selectedTissues", label = "Select Tissue/s of Origin", choices=c("none", tissueTypes),
  	# # 						multiple=TRUE, selected="none")
  	# # }
  	# # #
  })

  ## ------------------------------------------------------------------------
  showColorTissues <- reactive({
    tree2 <- input$tree2
    req(tree2)
    ## zz = get_selected(tree, format = "names")
    
    showColorTissues = names(unlist(get_selected(tree2, format = "slices")))
    showColorTissues = gsub("\\.",":", showColorTissues)
    showColorTissues = gsub("no_selection:","", showColorTissues)
    cat(showColorTissues,"\n")
    showColorTissues = setdiff(showColorTissues,"no_selection")
    return(sort(showColorTissues))
    
  })
  
  output$tree2 <- renderTree({
    
    shiny::validate(need(length(analysisTissueTypes()) > 0, 
                         "There are no cell lines of the selected tissue type(s)."))
    # end new
    matchedCellLinesTab = matchedCellLinesTab()
    
    tissueChoices <- getSampleSetTissueTypes(
      sampleSet = rownames(matchedCellLinesTab),
      ## sampleSet = rownames(req(matchedCellLinesTab())), 
      dataSource = input$xDataset, 
      srcContent = srcContentReactive()
    ) 
    # opt = "";
    # for(y in 1:length(tissueChoices)){
    #   # style works only for browser Chrome
    #   opt =  paste0(opt,"<option style='white-space: pre-wrap'>",tissueChoices[y],"</option>;");
    # }
    # HTML(
    #   paste("<label class='control-label' for='showColorTissues' id='lsc'>Select Tissues to Color</label>","<select id='showColorTissues' style='word-wrap:break-word; width: 100%;' multiple>",opt,"</select>")
    # )
    # print(tissueChoices)
    mylist2 <- convertListTree(tissueChoices)
    ## mylist2 = list(no_selection=structure(mylist2, stselected=TRUE))
    
    mylist2 = list(no_selection=structure(mylist2))
    if (! is.null(mylist2$no_selection$DataSet$NCI) ) attr(mylist2$no_selection$DataSet$NCI,"stselected")=TRUE
    if (! is.null(mylist2$no_selection$DataSet$TU) ) attr(mylist2$no_selection$DataSet$TU,"stselected")=TRUE
    if (! is.null(mylist2$no_selection$DataSet$UoC) ) attr(mylist2$no_selection$DataSet$UoC,"stselected")=TRUE
    if (! is.null(mylist2$no_selection$DataSet$UR) ) attr(mylist2$no_selection$DataSet$UR,"stselected")=TRUE
    #
    attr(mylist2$no_selection,"stopened")=TRUE
    attr(mylist2$no_selection$DataSet,"stopened")=TRUE
    mylist2
  })
    ## -------------------------------------------------------------------------
  output$showColorTissuesUi <- renderUI({
    
    shinyTree("tree2", checkbox = F, theme = "proton")
    
  # 	#tissueChoices <- analysisTissueTypes()
  #   #new
  #   shiny::validate(need(length(analysisTissueTypes()) > 0, 
  #                        "There are no cell lines of the selected tissue type(s)."))
  #   # end new
  #   matchedCellLinesTab = matchedCellLinesTab()
  #   
  # 	tissueChoices <- getSampleSetTissueTypes(
  # 	  sampleSet = rownames(matchedCellLinesTab),
  # 		## sampleSet = rownames(req(matchedCellLinesTab())), 
  # 		dataSource = input$xDataset, 
  # 		srcContent = srcContentReactive()
  # 		) 
  # 	opt = "";
  # 	for(y in 1:length(tissueChoices)){
  # 	    # style works only for browser Chrome
  # 	    opt =  paste0(opt,"<option style='white-space: pre-wrap'>",tissueChoices[y],"</option>;");
  # 	}
  # 	HTML(
  # 	  paste("<label class='control-label' for='showColorTissues' id='lsc'>Select Tissues to Color</label>","<select id='showColorTissues' style='word-wrap:break-word; width: 100%;' multiple>",opt,"</select>")
  # 	)
  # 	# selectInput("showColorTissues", "Select Tissues to Color",choices = tissueChoices, multiple = TRUE)
	})
  
  output$showCellsUi <- renderUI({
    #new
    #shiny::validate(need(length(analysisTissueTypes()) > 0, 
    #                     "There are no cell lines of the selected tissue type(s)."))
    # end new
    #tissueChoices <- analysisTissueTypes()
    # cellChoices <- getSampleSetTissueTypes(
    #   sampleSet = rownames(req(matchedCellLinesTab())), 
    #   dataSource = input$xDataset, 
    #   srcContent = srcContentReactive()
    # ) 
    
    dlDataTab <- getPlotData(xData = xData(), yData = yData(), showColor = input$showColor, 
                             showColorTissues = showColorTissues(), dataSource = input$xDataset, 
                             srcContent = srcContentReactive())
    
    
    
    shiny::validate(need(nrow(dlDataTab)>0, paste("ERROR:", " No common complete data found.")))
    cellChoices <- dlDataTab[,1]
    ### cellChoices <- rownames(req(matchedCellLinesTab()))
    opt = "";
    for(y in 1:length(cellChoices)){
      # style works only for browser Chrome
      opt =  paste0(opt,"<option style='white-space: pre-wrap'>",cellChoices[y],"</option>;");
    }
    # HTML(
    #   paste("<label class='control-label' for='showCells' id='lcells'>Select Patient sample to Color</label>","<select id='showCells' style='word-wrap:break-word; width: 100%;' multiple>",opt,"</select>")
    # )
    selectInput("showCells", "Select Patient sample to highlight",choices = cellChoices, multiple = TRUE, width='350px')
  })
  
  output$ipAddress <- renderText({
  	# debug
  	text <- readLines("http://api.ipify.org")
  })

  ## mutation variant ------------------------------------------------------------------
  output$mutationPanelV2 = renderUI({
    
    # plotPanel <- tabPanel("Plot Data", value=1, uiOutput("showCellsUi"), plotlyOutput("rChartsAlternative", width = plotWidth, height = plotHeight),
    #                       br(), br(), p("Plot point tooltips provide additional information."))
    # 
    muttab1 <- tabPanel("View variants", value=1 ,HTML("The following table displays the selected gene mutation variants for all patients with details such as variant allele frequency (VAF), mutation type (Exonic function) and Amino Acid changes (AAChange). The AAChange specfies the cDNA (c.) and the protein (p.) changes for all gene transcripts. All other details can be found by scrolling to the right of the screen."),
                        br(),br(),DT::dataTableOutput("mut_table_v2"))
    muttab2 <- tabPanel("Oncoplot", value=2,br(),HTML(paste0("<p style='text-align:center; font-size: 125%;'><b>Altered in ",mut_input()$nbpat," of total ",mut_input()$totpat, " samples</b></p>")), withSpinner(plotOutput("MutationPlot")))
    ##muttab2 <- tabPanel("Oncoplot", value=2, withSpinner(plotlyOutput("MutationPlotly", width = 1700, height = 1700)))
    
    mutPanel <- tabsetPanel(id="mp",muttab2, muttab1)
    
    
    return(mutPanel)
    
    
    
  })
  
  mut_input <- reactive({
    srcContent <- srcContentReactive()
    shiny::validate(need(!is.null(srcContent[[input$dataMut]][["molPharmData"]][["var"]]), "No mutation details found"))
    ## find gene
    # library(dplyr)
    # library(tidyr)
    # 
    vardat = srcContent[[input$dataMut]][["molPharmData"]][["var"]]
    dim(vardat) # 12997    89
    nbr = nrow(vardat)
    ina = apply(vardat,2,function(x) length(which(is.na(x))))
    j  = which(ina==nbr)
    
    vardat = data.frame(vardat[,-j], stringsAsFactors = F, check.names = F)
    # dim(vardat) # 12997    67
    ncol1 = ncol(vardat)
    
    varann = srcContent[[input$dataMut]][["molPharmData"]][["varA"]]
    # dim(varann)
    # 12997    55
    ncol2 = ncol(varann)
    # stopifnot(identical(rownames(vardat), rownames(varann)))
    # i= grep("slfn11", varann$Gene.refGene, ignore.case = T)
    ## new 
    selgenes = unlist(strsplit(trimws(input$geneId)," ")) ## many genes
    shiny::validate(need(length(selgenes)<6, "Please enter maximum 5 Gene symbols"))
    ##
    ## i = which(varann$Gene.refGene %in% toupper(input$geneId))
    i = which(varann$Gene.refGene %in% toupper(selgenes) )
    shiny::validate(need(length(i)>0, "Gene(s) not found"))
    res = cbind(varann[i,],vardat[i,])
    # cat(dim(res)) # 2 -122
    
    resu = pivot_longer(res, cols = (ncol2+1):(ncol1+ncol2), names_to ="sample_Id",  values_to = "VAF")
    # dim(resu) # 134 - 57
    resu = resu[,c(2:8,ncol2+1,ncol2+2,10:ncol2)]
    # dim(resu)
    
    k = which(resu$VAF==0)
    if (length(k)>0) resu = resu[-k,]
    
    nbpat = length(unique(resu$sample_Id))
    # dim(resu)
    # add Oncotree level 1 to 3
    sampleann = data.frame(srcContent[[input$dataMut]][["sampleData"]], stringsAsFactors = F, check.names = F)[,c(1,3,4,22,23)]
    colnames(sampleann)[1]="sample_Id"
    colnames(sampleann)[4]="dataSource"
    
    resu2 = merge(sampleann,resu, by.x="sample_Id", by.y="sample_Id")
    ### resu2 = resu2[,c(5:11,1:4,12:ncol(resu2))]
    # resu2 = resu2[,c(12,6:11,1:5,13:ncol(resu2))]
    resu2 = resu2[,c(14,6:13,1:5,15:ncol(resu2))] ## new on Dec17
    # create table with dplyr
    resu2$VAF = round(resu2$VAF,3)
    # return(resu2)
    return(list(resu2 =resu2, totpat = ncol1, nbpat = nbpat))
    
  })
  
  output$mut_table_v2 <- DT::renderDataTable({
    # ready
     resu2 = mut_input()$resu2
     # selsource=metaConfig[[input$dataMut]][["fullName"]]
  
    
    DT::datatable(resu2, rownames=TRUE,extensions='Buttons',
                  filter='top', style='bootstrap4', selection = "none",
                  options=list(
                               columnDefs = list(list(targets = c(0,8,12,15:ncol(resu2)), visible=FALSE )),
                               ## columnDefs = list(list(targets = c(0,9,12,15:ncol(resu2)), visible=FALSE )),
                               pageLength = 100,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipBt',scrollX = TRUE,buttons = list('colvis','copy', 'print', list(extend = 'collection',buttons = list(list(extend='csv',filename='TumorMiner_variants',title='Exported data from TumorMiner'), list(extend='excel',filename='TumorMiner_variants',title='Exported data from TumorMiner'), list(extend='pdf',filename='TumorMiner_variants',title='Exported data from TumorMiner')),text = 'Download') ))
               #   , caption=htmltools::tags$caption(paste0("Mutation variants details for ",selsource),style="color:dodgerblue; font-size: 18px")
    )
    
    
  })
  
  output$MutationPlot <- renderPlot({
    resu2 = mut_input()$resu2
    shiny::validate(need(length(unique(resu2$Gene.refGene))>1, "Need at least 2 genes to show oncoplot"))
    knowng = unique(resu2$Gene.refGene)
    # colnames(resu2)[1:14]
    # [1] "Gene.refGene"                  "Chr"                           "Start"                        
    # [4] "End"                           "Ref"                           "Alt"                          
    # [7] "Func.refGene"                  "sample_Id"                     "OncoTree1"                    
    # [10] "OncoTree2"                     "dataSource"                    "disease"                      
    # [13] "VAF"                           "ExonicFunc.refGene"            
    
  # colnames(resu2)[1:14] = c("Hugo_Symbol","Chromosome","Start_Position",
  #                       "End_Position","Reference_Allele","Tumor_Seq_Allele1", 
  #                       "Func.refGene", "Tumor_Sample_Barcode","OncoTree1",
  #                       "OncoTree2","dataSource", "disease" ,
  #                       "VAF", "Variant_Classification")
  #   
    colnames(resu2)[c(1:6,10)] = c("Hugo_Symbol","Chromosome","Start_Position",
                              "End_Position","Reference_Allele","Tumor_Seq_Allele1", 
                               "Tumor_Sample_Barcode") ## new dec17
 
     ## Variant_Classification can only be \{Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins,
    ###       Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, RNA, Targeted_Region}.
    ## variant type : Valid elements include: "SNP", "DNP", "TNP", "ONP", "DEL", "INS". Used to map frameshift_variant to more specific MAF columns (character)
    
    resu2$Tumor_Seq_Allele2 = NA
    
    ## Ignore -----------------
    # resu2$Variant_Type = ifelse(resu2$Variant_Classification %in% c("missense", "nonsense"), "SNP", "Indel") 
    # ## ?? read_through, Above should be updated !!!
    # 
    # resu2$Variant_Classification[which(resu2$Variant_Classification == "missense")] = "Missense_Mutation"
    # resu2$Variant_Classification[which(resu2$Variant_Classification == "nonsense")] = "Nonsense_Mutation"
    # resu2$Variant_Classification[which(resu2$Variant_Classification == "frameshift")] = "Frameshift"
    # resu2$Variant_Classification[which(resu2$Variant_Classification == "nonframeshift")] = "Nonframeshift"
    # resu2$Variant_Classification[which(resu2$Variant_Classification == "read_through")] = "Read_through"

    ## end ignore ------------
## above should be updated !!!
    
    # resuacc = resu2[which(resu2$disease == "ACC"),]
    # dim(resuacc) # 12234    18
    # clindata.acc  = unique(resuacc[,c("Tumor_Sample_Barcode", "disease")])
    # dim(clindata.acc) # 123 2
    
    # resu2 = resu2[order(resu2$disease),]
    
    ## clindata  = unique(resu2[,c("Tumor_Sample_Barcode", "disease")])
    clindata  = unique(resu2[,c("Tumor_Sample_Barcode", "dataSource")])
     
    dim(clindata) ## 332 2
    
    maf = read.maf(resu2, clinicalData = clindata, verbose = F)
    
    
    # mafacc = read.maf(resuacc, clinicalData = clindata.acc)
    # oncoplot(maf = mafacc, top = 10)  # Adjust 'top' to show the most frequent mutated genes
    
    # # Step 5: Optionally, customize the oncoplot with additional parameters
    # oncoplot(
    #   maf = maf,
    #   top = 10,  # Display the top 10 mutated genes
    #   mutationClassOrder = c("Missense_Mutation", "Frameshift_Del", "Nonsense_Mutation"),  # Order mutations by class
    #   drawBarPlot = TRUE,  # Draw a bar plot of mutation frequencies
    #   showTumorSampleBarcodes = FALSE  # Hide sample barcodes
    # )
    
    # oncoplot(maf = maf, clinicalFeatures = 'disease', genes = knowng, 
    #          fontSize = 1, legendFontSize = 1.5, annotationFontSize = 1.5, titleFontSize = 1.5) 
    # 
  # par(mar = c(2, 2, 1, 1))  # Adjust margins
  # oncoplot(maf = maf, clinicalFeatures = 'disease', top=5, sortByAnnotation = T,
  #                  fontSize = 1, legendFontSize = 2.5, annotationFontSize = 2.5, 
  #                  titleFontSize = 1.8, drawColBar = T
  #            )
  
  # HTML( paste0("Altered in ",nrow(clindata)," of total ",mut_input()$totpat, " samples") )

  par(mar = c(1, 1, 10, 1))  # Adjust margins
  oncoplot(maf = maf, clinicalFeatures = 'dataSource', top=5, sortByAnnotation = T, barcode_mar = 10,
           fontSize = 1.2, legendFontSize = 2.5, annotationFontSize = 2.5, gene_mar = 7, showTitle = F ,
           titleFontSize = 2, drawColBar = F , drawRowBar = T , titleText = paste0("Altered in ",nrow(clindata)," of total ",mut_input()$totpat, " samples")
  )
  

  
    
  } ,  width=1200, height=650 # was 1200 and 650
  )  
  
  output$MutationPlotly <- renderPlotly({
    resu2 = mut_input()
    ## shiny::validate(need(length(unique(resu2$Gene.refGene))>1, "Need at least 2 genes to show oncoplot"))
    knowng = unique(resu2$Gene.refGene)
    colnames(resu2)
    # [1] "Gene.refGene"                  "Chr"                           "Start"                        
    # [4] "End"                           "Ref"                           "Alt"                          
    # [7] "Func.refGene"                  "sample_Id"                     "OncoTree1"                    
    # [10] "OncoTree2"                     "dataSource"                    "disease"                      
    # [13] "VAF"                           "ExonicFunc.refGene"            
    
    colnames(resu2)[1:14] = c("Hugo_Symbol","Chromosome","Start_Position",
                              "End_Position","Reference_Allele","Tumor_Seq_Allele1", 
                              "Func.refGene", "Tumor_Sample_Barcode","OncoTree1",
                              "OncoTree2","dataSource", "disease" ,
                              "VAF", "Variant_Classification")
    
    ## Variant_Classification can only be \{Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins,
    ###       Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, RNA, Targeted_Region}.
    ## variant type : Valid elements include: "SNP", "DNP", "TNP", "ONP", "DEL", "INS". Used to map frameshift_variant to more specific MAF columns (character)
    
    resu2$Tumor_Seq_Allele2 = NA
    
    
    resu2$Variant_Type = ifelse(resu2$Variant_Classification %in% c("missense", "nonsense"), "SNP", "Indel") 
    ## ?? read_through, Above should be updated !!!
    
    resu2$Variant_Classification[which(resu2$Variant_Classification == "missense")] = "Missense_Mutation"
    resu2$Variant_Classification[which(resu2$Variant_Classification == "nonsense")] = "Nonsense_Mutation"
    resu2$Variant_Classification[which(resu2$Variant_Classification == "frameshift")] = "Frameshift"
    resu2$Variant_Classification[which(resu2$Variant_Classification == "nonframeshift")] = "Nonframeshift"
    resu2$Variant_Classification[which(resu2$Variant_Classification == "read_through")] = "Read_through"
    ## above should be updated !!!
    
    # resuacc = resu2[which(resu2$disease == "ACC"),]
    # dim(resuacc) # 12234    18
    # clindata.acc  = unique(resuacc[,c("Tumor_Sample_Barcode", "disease")])
    # dim(clindata.acc) # 123 2
    
    clindata  = unique(resu2[,c("Tumor_Sample_Barcode", "disease")])
    dim(clindata) ## 332 2
    
    maf = read.maf(resu2, clinicalData = clindata, verbose = F)
    
    
    # mafacc = read.maf(resuacc, clinicalData = clindata.acc)
    # oncoplot(maf = mafacc, top = 10)  # Adjust 'top' to show the most frequent mutated genes
    
    # # Step 5: Optionally, customize the oncoplot with additional parameters
    # oncoplot(
    #   maf = maf,
    #   top = 10,  # Display the top 10 mutated genes
    #   mutationClassOrder = c("Missense_Mutation", "Frameshift_Del", "Nonsense_Mutation"),  # Order mutations by class
    #   drawBarPlot = TRUE,  # Draw a bar plot of mutation frequencies
    #   showTumorSampleBarcodes = FALSE  # Hide sample barcodes
    # )
    
    # oncoplot(maf = maf, clinicalFeatures = 'disease', genes = knowng, 
    #          fontSize = 1, legendFontSize = 1.5, annotationFontSize = 1.5, titleFontSize = 1.5) 
    # 
    #par(mar = c(3, 3, 3, 1))  # Adjust margins
    gg <- oncoplot(maf = maf, clinicalFeatures = 'disease', top=5, 
                   fontSize = 1, legendFontSize = 1.5, annotationFontSize = 1.5, titleFontSize = 1.5)
    igg <- ggplotly(gg
                    # ,plotWidth=1200,plotHeight=650
    )
    igg
    
  } #,  width=1200, height=650
  )  
    
  output$mutationPanel = renderUI({
    DT::dataTableOutput("mut_table")
    
  })
  
  output$mut_table <- DT::renderDataTable({
    srcContent <- srcContentReactive()
    shiny::validate(need(!is.null(srcContent[[input$dataMut]][["molPharmData"]][["var"]]), "No mutation details found"))
    ## find gene
    library(dplyr)
    library(tidyr)
    
    vardat = srcContent[[input$dataMut]][["molPharmData"]][["var"]]
    dim(vardat) # 12997    89
    nbr = nrow(vardat)
    ina = apply(vardat,2,function(x) length(which(is.na(x))))
    j  = which(ina==nbr)
    
    vardat = data.frame(vardat[,-j], stringsAsFactors = F, check.names = F)
    # dim(vardat) # 12997    67
    ncol1 = ncol(vardat)
    
    varann = srcContent[[input$dataMut]][["molPharmData"]][["varA"]]
    # dim(varann)
    # 12997    55
    ncol2 = ncol(varann)
    # stopifnot(identical(rownames(vardat), rownames(varann)))
    # i= grep("slfn11", varann$Gene.refGene, ignore.case = T)
    
    i = which(varann$Gene.refGene %in% toupper(input$geneId))
    shiny::validate(need(length(i)>0, "Gene not found"))
    res = cbind(varann[i,],vardat[i,])
    # cat(dim(res)) # 2 -122
    
    resu = pivot_longer(res, cols = (ncol2+1):(ncol1+ncol2), names_to ="sample_Id",  values_to = "VAF")
    # dim(resu) # 134 - 57
    resu = resu[,c(2:8,ncol2+1,ncol2+2,10:ncol2)]
    # dim(resu)
    
    k = which(resu$VAF==0)
    if (length(k)>0) resu = resu[-k,]
    # dim(resu)
    # add Oncotree level 1 to 3
    sampleann = data.frame(srcContent[[input$dataMut]][["sampleData"]], stringsAsFactors = F, check.names = F)[,c(1,3,4,22)]
    colnames(sampleann)[1]="sample_Id"
    colnames(sampleann)[4]="dataSource"
    
    resu2 = merge(sampleann,resu, by.x="sample_Id", by.y="sample_Id")
    resu2 = resu2[,c(5:11,1:4,12:ncol(resu2))]
    # create table with dplyr
    resu2$VAF = round(resu2$VAF,3)
    # ready
    selsource=metaConfig[[input$dataMut]][["fullName"]]
    # DT::datatable(resu2, rownames=FALSE,extensions='Buttons',
    #               filter='top', style='bootstrap4', selection = "none",
    #               options=list(pageLength = 100,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipBt',buttons = list('copy', 'print', list(extend = 'collection',buttons = list(list(extend='csv',filename='Sarcoma_variants',title='Exported data from CellMinerCDB'), list(extend='excel',filename='Sarcoma_variants',title='Exported data from CellMinerCDB'), list(extend='pdf',filename='Sarcoma_variants',title='Exported data from CellMinerCDB')),text = 'Download')))
    #               , caption=htmltools::tags$caption(paste0("Mutation variants details for ",selsource),style="color:dodgerblue; font-size: 18px")
    # )

    DT::datatable(resu2, rownames=TRUE,extensions='Buttons',
                  filter='top', style='bootstrap4', selection = "none",
                  options=list(columnDefs = list(list(targets = c(0,6,7,14:ncol(resu2)), visible=FALSE )),pageLength = 100,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipBt',buttons = list('colvis','copy', 'print', list(extend = 'collection',buttons = list(list(extend='csv',filename='TumorMiner_variants',title='Exported data from TumorMiner'), list(extend='excel',filename='TumorMiner_variants',title='Exported data from TumorMiner'), list(extend='pdf',filename='TumorMiner_variants',title='Exported data from TumorMiner')),text = 'Download') ))
                  , caption=htmltools::tags$caption(paste0("Mutation variants details for ",selsource),style="color:dodgerblue; font-size: 18px")
    )
    
    
  })
  
  output$selectInputSurvLassoUI <- renderUI({
    # ns <- session$ns
    ## selectInput(ns("inputGeneSets"), "Select Gene Sets",
    #choices  = c(names(geneSetPathwayAnalysis::geneSets), "All Genes"),
    ##					choices = names(geneSetPathwayAnalysis::geneSets),
    ##					selected = "All Gene Sets",
    ##					multiple=TRUE)
    ## choices = names(geneSetPathwayAnalysis::geneSets)
    choices = list.genesets
    if (paste(input$onco3ch, collapse = "") == "ACC_TCGA") {
      choices = c("ACC variable genes (4997)", choices) ##  New ACC
      ## mygeneset = append(accgeneset, mygeneset)
    }
    mysel="All Gene Sets (8663)"
    opt = "";
    for(y in 1:length(choices)){
      # style works only for browser Chrome
      if (choices[y]==mysel)
      {
        opt =  paste0(opt,"<option style='white-space: pre-wrap' selected>",choices[y],"</option>;");
      }
      else {
        opt =  paste0(opt,"<option style='white-space: pre-wrap'>",choices[y],"</option>;");
      }
    }
    HTML(
      paste("<label class='control-label' for='inputGSsurv'>Select Gene Sets</label>","<select id='inputGSsurv', size='6', style='word-wrap:break-word; width: 100%;' multiple>",opt,"</select>")
      
       )
    
  })
  
  
  #----[observers]-----------------------------------------------------------------------

  # Observe reactive variable and send message to Javascript code
  observe({
  	if(isPackageLoadingComplete()) {
  if (is.null(appConfig$modal))
  	  session$sendCustomMessage(type='showLoading', list(show=FALSE))
  	 else 
  	  session$sendCustomMessage(type='showSkip', list(show=TRUE))
     }
  })
	
  #-----[NavBar Tab Server Code]---------------------------------------------------------
  rm <- callModule(regressionModels, "rm", srcContentReactive = srcContentReactive,
  								 appConfig = appConfig, oncolor=oncolor)
  
  #-----[NavBar Tab Server Code2]---------------------------------------------------------
  myPatient <- callModule(myModule, "myPatient", srcContent = srcContent)
  
  ## ******************************************************************************************
  ## --------------- Functions and UI specific to Predictive biomarkers -----------
  ## ******************************************************************************************
  
  ## main TAB ---------------
  output$tabsetPanelpb = renderUI({
    
    # tab4 <- tabPanel("Tissue Correlation", value=4,
    #                  #downloadLink("downloadData", "Download selected x and y axis data as a Tab-Delimited File"),
    #                  DT::dataTableOutput("cortable"))
    tab1 <- tabPanel("View Data", value=2,
                     downloadLink("downloadDatapb", "Download selected x and y axis data as a Tab-Delimited File"),
                     DT::dataTableOutput("tablepb"))
    # tab2 <- tabPanel("Search IDs",
    #                  includeMarkdown("www/files/help.md"),
    #                  DT::dataTableOutput("ids"))
    tab3 <- tabPanel("Compare Patterns", value=3,
                     includeMarkdown("www/files/help.md"),
                     #br(),
                     HTML("<b>Pattern comparison results are computed with respect to that data defined and shared by both the x and y-axis inputs.</b>"),
                     br(),br(),
                     fluidRow(
                       #column(3, selectInput("patternComparisonType", "Pattern Comparison",
                       #           						choices=c("Molecular Data"="molData", "Drug Data"="drug"), 
                       #											selected="molData")),
                       
                       column(4, HTML(
                         paste("<label class='control-label' for='patternComparisonTypepb'>Select molecular or activity data</label>","<select id='patternComparisonTypepb'><option value='moldata' selected>Molecular Data</option><option value='drug'>Drug Data</option></select>")
                       )),
            
                       column(8, radioButtons("crossdbpb", label = NULL, choices = list("Compare x-Axis input to x-Axis molecular or activity data" = "No", "Compare x-Axis input to y-Axis molecular or activity data" = "Yes"), selected  = "No", inline=F, width="100%")       
                       )
                       
                     ),
                     br(),
                     renderUI({
                       req(PatternCompTablepb())
                       downloadLink("downloadDataComppb", "Download All as a Tab-Delimited File")
                     }),
                     ##downloadLink("downloadDataComp", "Download All as a Tab-Delimited File"),
                     withSpinner(DT::dataTableOutput("patternComparisonpb")))
    
     
    
    plotPanel <- tabPanel("Plot Data", value=1, uiOutput("showCellsUipb"), plotlyOutput("rChartsAlternativepb", width = plotWidth, height = plotHeight),
                          br(), br(), p("Plot point tooltips provide additional information."))
    
    tsPanel <- tabsetPanel(id="tspb",plotPanel, tab1, tab3)
    #		}
    
    return(tsPanel)
  })
  
  ## end main TAB ------------------------------------------------
  
  
  output$downloadDataComppb <- downloadHandler(
    
    
    filename = function() {
      
      pcDataset <- input$xDatasetpb
      pcType <- input$xPrefixpb
      pcId <- input$xIdpb
     
      pcDataset2 <- input$yDatasetpb
      if (input$crossdbpb == "No") {  
        ##paste0("Pattern_Comp_all_cdb_",input$patternComparisonSeed,"_",pcDataset,"_",pcType,"_",pcId,"_",input$patternComparisonType,".txt")
        paste0("Pattern_Comp_all_biomarker_",pcDataset,"_",pcType,"_",pcId,"_",input$patternComparisonTypepb,".txt") 
      }
      else {
        paste0("Pattern_Comp_all_biomarker_",pcDataset,"_",pcType,"_",pcId,"_",input$patternComparisonTypepb,"_from_",pcDataset2,".txt") 
      }
    },
    
    content = function(file) {
      
      write.table(PatternCompTablepb(), file, sep = "\t", row.names = F,quote=F)  
      
    }
  )
  
  
  output$patternComparisonpb <- DT::renderDataTable({
    results= PatternCompTablepb()	
    
    DT::datatable(results, rownames=FALSE, colnames=colnames(results),
                  filter='top', style='bootstrap4', selection = "none",
                  options=list(lengthMenu = c(10, 50, 100,500), pageLength = 100,language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')) ,dom='lipt'))
    
  })
  
  
  ##------------
  PatternCompTablepb <- reactive({
    srcContent <- srcContentReactive()
    ## new
    if (input$crossdbpb == "Yes")
    {
      # if (input$patternComparisonSeed == "xPattern"){
      dat <- xDatapb()
      pcDataset <- input$yDatasetpb
      # selectedLines <- names(yData()$data)
      selectedLines <- as.character(matchedCellLinesTabpb()[, "yDatasetpb"])
      # } else{
      #   dat <- yData()
      #   pcDataset <- input$xDataset
      #   # selectedLines <- names(xData()$data)
      #   selectedLines <- as.character(matchedCellLinesTab()[, "xDataset"])
      # }
      names(dat$data) = selectedLines
    }
    
    else
    {
      
      # if (input$patternComparisonSeed == "xPattern"){
      dat <- xDatapb()
      pcDataset <- input$xDatasetpb
      # } else{
      #   dat <- yData()
      #   pcDataset <- input$yDataset
      # }
      selectedLines <- names(dat$data)
    } # end new
    
    
    shiny::validate(need(length(selectedLines)>0, paste("ERROR:", " No common complete data found.")))
    shiny::validate(need(length(selectedLines)>2, paste("ERROR:", " No display for less than 3 observations.")))
    
    if(input$patternComparisonTypepb == "drug") {
      shiny::validate(need(srcContent[[pcDataset]][["molPharmData"]][["act"]], "No drug available for this cell line set"))
      #if (is.null(srcContent[[pcDataset]][["molPharmData"]][["act"]])) stop("No drug available for this cell line set")
      results <- patternComparison(dat$data,
                                   srcContent[[pcDataset]][["molPharmData"]][["act"]][, selectedLines, drop=FALSE])
      results$ids <- rownames(results)
      results$NAME <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "NAME"]
      # OLD
      # if ("MOA" %in% colnames(srcContent[[pcDataset]][["drugInfo"]])){
      #   results$MOA <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "MOA"]
      #   results <- results[, c("ids", "NAME", "MOA", "COR", "PVAL")]
      #   colnames(results) <- c("ID", "Name", "MOA", "Correlation", "P-Value")
      # } else{
      #   results <- results[, c("ids", "NAME", "COR", "PVAL")]
      #   colnames(results) <- c("ID", "Name", "Correlation", "P-Value")
      # }
      
      #NEW
      results$MOA <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "MOA"]
      results$CLINICAL.STATUS <- srcContent[[pcDataset]][["drugInfo"]][rownames(results), "CLINICAL.STATUS"]
      results <- results[, c("ids", "NAME", "MOA", "CLINICAL.STATUS", "COR", "PVAL")]
      colnames(results) <- c("ID", "Name", "MOA", "CLINICAL.STATUS", "Correlation", "P-Value")
      # end new
      
      DataType <- getMolDataType(results$ID)
      DrugID <- removeMolDataType(results$ID)
      results$ID <- NULL
      results$FDR=p.adjust(results[,"P-Value"],method="BH",nrow(results))
      results=cbind(DataType,DrugID,results)
      colnames(results)[1:2]=c("Data Type", "Drug ID")
      
    } else {
      molPharmData <- srcContent[[pcDataset]][["molPharmData"]]
      molData <- molPharmData[setdiff(names(molPharmData), c("act","copA","mutA","metA","expA","xaiA","proA","mirA","mdaA","swaA","xsqA","mthA","hisA","criA","mtbA","rrbA","tmmA","tpmA","varA","var","mucA","hllA","sigA","nmfA"))]
      shiny::validate(need(length(molData)>0, "No molecular data available for this cell line set"))
      ##if (length(molData)==0) stop("No molecular data available for this cell line set")
      ## old: molData <- lapply(molData, function(X) X[, selectedLines])
      ## new selection in case a dataset has only one row
      ## one solution: molData <- lapply(molData, function(X) subset(X, select=selectedLines))
      molData <- lapply(molData, function(X) X[, selectedLines,drop=FALSE])
      results <- patternComparison(dat$data, molData)
      results$ids <- rownames(results)
      
      results$molDataType <- getMolDataType(results$ids)
      results$gene <- removeMolDataType(results$ids)
      
      # Reorder columns
      results <- results[, c("ids", "molDataType", "gene", "COR", "PVAL")]
      colnames(results) <- c("ID", "Data Type", "Gene", "Correlation", "P-Value")
      
      if (require(rcellminerUtilsCDB)){
        chromLocs <- character(nrow(results))
        haveLoc <- results$Gene %in% names(geneToChromBand)
        chromLocs[haveLoc] <- geneToChromBand[results$Gene[haveLoc]]
        
        results$Location <- chromLocs
        results <- results[, c("ID", "Data Type", "Gene", "Location", "Correlation", "P-Value")]
      }
      results$FDR=p.adjust(results[,"P-Value"],method="BH",nrow(results))
      
      if (require(geneSetPathwayAnalysis)){
        # old :results$Annotation <- geneSetPathwayAnalysis::geneAnnotTab[results$Gene, "SHORT_ANNOT"]
        # issue related to prefix subsetting ex: "age" >> gives "ager"
        #st=Sys.time()
        results$Annotation <- geneSetPathwayAnalysis::geneAnnotTab[match(results$Gene,rownames(geneSetPathwayAnalysis::geneAnnotTab)), "SHORT_ANNOT"]
        #et=Sys.time()
        #cat(et-st,"\n")
        results$Annotation[is.na(results$Annotation)] <- ""
      }
      
      results$ID <- NULL
      colnames(results)[2]="ID"
    }
    
    results[, "Correlation"] <- round(results[, "Correlation"], 3)
    results[, "P-Value"] <- signif(results[, "P-Value"], 3)
    results[, "FDR"] <- signif(results[, "FDR"], 3)
    ## sort by p-value
    results <- results[order(results[, "P-Value"]),]
    
    return(results)
  })
  ##--------------
  
  output$downloadDatapb <- downloadHandler(
    filename = function() {
      query <- parseQueryString(session$clientData$url_search)
      
      if("filename" %in% names(query)) {
        filename <- query[["filename"]]
      } else {
        filename <- "Biomarker_dataset"
      }
      
      if("extension" %in% names(query)) {
        extension <- query[["extension"]]
      } else {
        extension <- "txt"
      }
      
      paste(filename, extension, sep=".")
    },
    content = function(file) {
      df <- getPlotData(xData = xDatapb(), yData = yDatapb(), showColor = input$showColorpb, 
                        showColorTissues = showColorTissuespb(), dataSource = input$xDatasetpb, 
                        srcContent = srcContentReactive())
      
      # Column selection below is to restrict to cell line, x, y features,
      # and tissue type information (source-provided + OncoTree).
      
      dfCols <- c("Cell Line","patientID","paired_sample",colnames(df)[2:3], paste0("OncoTree", 1:2))
      ### ----------
      if ("EMT" %in% colnames(df)) {
        df[,"EMT"] <- gsub("EMT:","",df[,"EMT"])
        dfCols <- c(dfCols, "EMT")
      }
      
      if ("priorTreatment" %in% colnames(df)) {
        df[,"priorTreatment"] <- gsub("PriorTreatment:","",df[,"priorTreatment"])
        dfCols <- c(dfCols, "priorTreatment")
      }
      
      
      # new stuff
      if ("sampleType" %in% colnames(df)) {
        df[,"sampleType"] <- gsub("SampleType:","",df[,"sampleType"])
        dfCols <- c(dfCols, "sampleType")
      }
      if ("biopsySite" %in% colnames(df)) {
        df[,"biopsySite"] <- gsub("BiopsySite:","",df[,"biopsySite"])
        dfCols <- c(dfCols, "biopsySite")
      }
      # if ("molecularSubtype" %in% colnames(df)) {
      #   df[,"molecularSubtype"] <- gsub("MolecularSubtype:","",df[,"molecularSubtype"])
      #   dfCols <- c(dfCols, "molecularSubtype")
      # }
      if ("dataSet" %in% colnames(df)) {
        df[,"dataSet"] <- gsub("DataSet:","",df[,"dataSet"])
        dfCols <- c(dfCols, "dataSet")
      }
      
      df <- df[, dfCols]
      # df[, 2] <- round(df[, 2], 3)
      # df[, 3] <- round(df[, 3], 3)
      colnames(df)[1]="Patient_sample"
      
      ### ----------
      # dfCols <- c(colnames(df)[1:4], paste0("OncoTree", 1:4))
      # if ("EMT" %in% colnames(df)) {
      #   dfCols <- c(dfCols, "EMT")
      # }
      # df <- df[, dfCols]
      
      write.table(df, file, quote=FALSE, row.names=FALSE, sep="\t")
    }
  )
  
  output$tablepb <- DT::renderDataTable({
    # Column selection below is to restrict to cell line, x, y features,
    # and tissue type information (source-provided + OncoTree).
    dlDataTab <- getPlotData(xData = xDatapb(), yData = yDatapb(), showColor = input$showColorpb, 
                             showColorTissues = showColorTissuespb(), dataSource = input$xDatasetpb, 
                             srcContent = srcContentReactive())
    
    # cat(colnames(dlDataTab),"\n")
    
    shiny::validate(need(nrow(dlDataTab)>0, paste("ERROR:", " No common complete data found.")))
    # new stuff
    dlDataTab[, 2] <- round(dlDataTab[, 2], 3)
    dlDataTab[, 3] <- round(dlDataTab[, 3], 3)	
    dlDataTabCols <- c(colnames(dlDataTab)[c(1,21,22,2:4)], paste0("OncoTree", 1:4)) # adding patient ID
    #
    ## dlDataTabCols <- c(colnames(dlDataTab)[1:4], paste0("OncoTree", 1:4))
    if ("EMT" %in% colnames(dlDataTab)) {
      dlDataTab[,"EMT"] <- gsub("EMT:","",dlDataTab[,"EMT"])
      dlDataTabCols <- c(dlDataTabCols, "EMT")
    }
    
    if ("priorTreatment" %in% colnames(dlDataTab)) {
      dlDataTab[,"priorTreatment"] <- gsub("PriorTreatment:","",dlDataTab[,"priorTreatment"])
      dlDataTabCols <- c(dlDataTabCols, "priorTreatment")
    }
    
    if ("NAPY" %in% colnames(dlDataTab)) {
      dlDataTabCols <- c(dlDataTabCols, "NAPY")
    }
    # new stuff
    if ("sampleType" %in% colnames(dlDataTab)) {
      dlDataTab[,"sampleType"] <- gsub("SampleType:","",dlDataTab[,"sampleType"])
      dlDataTabCols <- c(dlDataTabCols, "sampleType")
    }
    if ("biopsySite" %in% colnames(dlDataTab)) {
      dlDataTab[,"biopsySite"] <- gsub("BiopsySite:","",dlDataTab[,"biopsySite"])
      dlDataTabCols <- c(dlDataTabCols, "biopsySite")
    }
    if ("molecularSubtype" %in% colnames(dlDataTab)) {
      dlDataTab[,"molecularSubtype"] <- gsub("MolecularSubtype:","",dlDataTab[,"molecularSubtype"])
      dlDataTabCols <- c(dlDataTabCols, "molecularSubtype")
    }
    if ("dataSet" %in% colnames(dlDataTab)) {
      dlDataTab[,"dataSet"] <- gsub("DataSet:","",dlDataTab[,"dataSet"])
      dlDataTabCols <- c(dlDataTabCols, "dataSet")
    }
    
    dlDataTab <- dlDataTab[, dlDataTabCols]
    # dlDataTab[, 2] <- round(dlDataTab[, 2], 3)
    # dlDataTab[, 3] <- round(dlDataTab[, 3], 3)
    colnames(dlDataTab)[1]="Patient_sample"
    ## new filtering
    dlDataTab = dlDataTab[, -c(6,9,10,15)] # remove tissue, Onco3, Onco4, molecular subtype
    DT::datatable(dlDataTab, rownames=FALSE, colnames=colnames(dlDataTab),
                  filter='top', style='bootstrap4', selection="none",
                  extensions = "FixedColumns",
                  ## options=list(pageLength = nrow(dlDataTab), language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')))
                  options=list(scrollX = TRUE, fixedColumns = list(leftColumns = 1), lengthMenu = c(10, 25, 50, 100), pageLength = 10, language=list(paginate = list(previous = 'Previous page', `next`= 'Next page')))
    )
  })
  
  
  
  output$showCellsUipb <- renderUI({
    #new
   
    
    dlDataTab <- getPlotData(xData = xDatapb(), yData = yDatapb(), showColor = input$showColorpb, 
                             showColorTissues = showColorTissuespb(), dataSource = input$xDatasetpb, 
                             srcContent = srcContentReactive())
    
    
    
    shiny::validate(need(nrow(dlDataTab)>0, paste("ERROR:", " No common complete data found.")))
    cellChoices <- dlDataTab[,1]
    ### cellChoices <- rownames(req(matchedCellLinesTab()))
    opt = "";
    for(y in 1:length(cellChoices)){
      # style works only for browser Chrome
      opt =  paste0(opt,"<option style='white-space: pre-wrap'>",cellChoices[y],"</option>;");
    }
    # HTML(
    #   paste("<label class='control-label' for='showCells' id='lcells'>Select Patient sample to Color</label>","<select id='showCells' style='word-wrap:break-word; width: 100%;' multiple>",opt,"</select>")
    # )
    selectInput("showCellspb", "Select Patient sample to highlight",choices = cellChoices, multiple = TRUE, width='350px')
  })
  
  
  output$rChartsAlternativepb <- renderPlotly({
    #-----[range check]----------------------------------------------------------
    # Note: Until a better solution can be found, these checks are needed.
    # The issue is that upon a data source change, there appears to be a moment 
    # when the xData() or yData() are updated, but the input$xAxisRange
    # or input$yAxisRange (from the sliderInput UI element) are not yet updated.
    # As such, the data value range can be out of synch with the invalidated
    # axis ranges. In the most extreme case, there are no points in the
    # specified range. The reactive code is quickly re-run with the proper
    # inputs, correcting the plot, but the error flashes briefly in a buggy 
    # looking way. 
    # Below we do a range check and quietly exit if something is amiss (knowing
    # that the reactivity will ensure that the code is re-run with a proper
    # set of inputs once thing settle down).
    #****************************************************************************
    xData <- xDatapb()
    yData <- yDatapb()
    shiny::validate(need(xData$uniqName != yData$uniqName, 
                         "Please select distinct x and y axis variables."))
    
    xValRange <- range(xData$data, na.rm = TRUE)
    xLimits <- input$xAxisRangepb
    
    yValRange <- range(yData$data, na.rm = TRUE)
    yLimits <- input$yAxisRangepb
    
    # req(FALSE) causes immediate but quiet exit.
    req(!any(is.null(xValRange)), !any(is.null(xLimits)))
    req(!any(is.null(yValRange)), !any(is.null(yLimits)))
    req(!any(is.na(xValRange)), !any(is.na(xLimits)))
    req(!any(is.na(yValRange)), !any(is.na(yLimits)))
    # commented below
    # req((xLimits[1] <= xValRange[1]) && (xValRange[2] <= xLimits[2]))
    # req((yLimits[1] <= yValRange[1]) && (yValRange[2] <= yLimits[2]))
    
    cat("xAxis Limits: ", paste0(xLimits, collapse = " "), sep = "\n")
    cat("X_VAL_RANGE: ",  paste0(xValRange, collapse = " "), sep = "\n")
    
    cat("yAxis Limits: ", paste0(yLimits, collapse = " "), sep = "\n")
    cat("Y_VAL_RANGE: ",  paste0(yValRange, collapse = " "), sep = "\n")
    cat("-------------------------------------------", sep = "\n")
    #----------------------------------------------------------------------------
    
    p1 <- makePlotStatic(xData = xData, yData = yData, showColor = input$showColorpb, 
                         showColorTissues = showColorTissuespb(), dataSource = input$xDatasetpb, 
                         xLimVals = xLimits, yLimVals = yLimits,
                         srcContent = srcContentReactive(),oncolor=oncolor, showCells=input$showCellspb)
    p1 <- p1 + theme(axis.text = element_text(size=16), plot.title = element_text(size = 16, hjust = 0.5), 
                     axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
    #theme_update(legend.position = c(0,0))
    g1 <- ggplotly(p1, width=plotWidth, height=plotHeight, tooltip=tooltipCol)
    #g1 <- layout(g1, margin=list(t = 75))
    g1 <- layout(g1, margin=list(t = 75), legend = list(font = list(size = 18)))
    g2 <- config(p = g1,  cloud=FALSE, displaylogo=FALSE, displayModeBar=TRUE,
                 modeBarButtonsToRemove=c("select2d", "sendDataToCloud", "pan2d", "resetScale2d",
                                          "hoverClosestCartesian", "hoverCompareCartesian",
                                          "lasso2d", "zoomIn2d", "zoomOut2d"))
    g2
  })
  
  
  showColorTissuespb <- reactive({
    tree2pb <- input$tree2pb
    req(tree2pb)
    ## zz = get_selected(tree, format = "names")
    
    showColorTissues = names(unlist(get_selected(tree2pb, format = "slices")))
    showColorTissues = gsub("\\.",":", showColorTissues)
    showColorTissues = gsub("no_selection:","", showColorTissues)
    cat(showColorTissues,"\n")
    showColorTissues = setdiff(showColorTissues,"no_selection")
    return(sort(showColorTissues))
    
  })
  
  output$tree2pb <- renderTree({
    
    shiny::validate(need(length(analysisTissueTypespb()) > 0, 
                         "There are no cell lines of the selected tissue type(s)."))
    # end new
    matchedCellLinesTab = matchedCellLinesTabpb()
    
    tissueChoices <- getSampleSetTissueTypes(
      sampleSet = rownames(matchedCellLinesTab),
      ## sampleSet = rownames(req(matchedCellLinesTab())), 
      dataSource = input$xDatasetpb, 
      srcContent = srcContentReactive()
    ) 
   
    mylist2 <- convertListTree(tissueChoices)
    # mylist2 = list(no_selection=structure(mylist2, stselected=TRUE))
    # mylist2
    
    mylist2 = list(no_selection=structure(mylist2))
    if (! is.null(mylist2$no_selection$DataSet$NCI) ) attr(mylist2$no_selection$DataSet$NCI,"stselected")=TRUE
    if (! is.null(mylist2$no_selection$DataSet$TU) ) attr(mylist2$no_selection$DataSet$TU,"stselected")=TRUE
    if (! is.null(mylist2$no_selection$DataSet$UoC) ) attr(mylist2$no_selection$DataSet$UoC,"stselected")=TRUE
    if (! is.null(mylist2$no_selection$DataSet$UR) ) attr(mylist2$no_selection$DataSet$UR,"stselected")=TRUE
    
    # attr(mylist2$no_selection$DataSet$NCI,"stselected")=TRUE
    # attr(mylist2$no_selection$DataSet$TU,"stselected")=TRUE
    # attr(mylist2$no_selection$DataSet$UoC,"stselected")=TRUE
    # attr(mylist2$no_selection$DataSet$UR,"stselected")=TRUE
    #
    attr(mylist2$no_selection,"stopened")=TRUE
    attr(mylist2$no_selection$DataSet,"stopened")=TRUE
    mylist2
  })
  
  
  output$showColorTissuesUipb <- renderUI({
    
    shinyTree("tree2pb", checkbox = F, theme = "proton")
    
    #
  })
  
  
  output$selectTissuesUipb <- renderUI({
   
    
    # HTML("<p><b>Select Tissue/s of Origin</b></p>")
    shinyTree("treepb", checkbox = F, theme = "proton")
    ## end very new ------------------------------------------------------------------- ###
    
   
    # # #
  })
  
  
  output$treepb <- renderTree({
    srcContent <- srcContentReactive()
    tissueToSamplesMap <- srcContent[[input$xDatasetpb]][["tissueToSamplesMap"]]
    ## tissueTypes <- sort(unique(names(tissueToSamplesMap)))
    tissueTypes <- unique(names(tissueToSamplesMap)) # no sorting
    ## very new ------------------------------------------------------------------------###
    mylist <- convertListTree(tissueTypes)
    if (input$tissueSelectionModepb == "To include"){
      mylist = list(all=structure(mylist, stselected=TRUE))
      
      
    } else{ # input$tissueSelectionMode == "To exclude"
      mylist = list(none=structure(mylist, stselected=TRUE))
      
    }
    mylist
  })
  
  yDatapb <- reactive({
    ## shiny::validate(need(length(input$selectedTissues) > 0, "Please select tissue types."))
    shiny::validate(need(length(selectedTissuespb()) > 0, "Please select tissue types."))
    ## new
    shiny::validate(need(!is.na(match(input$yPrefixpb,srcContentReactive()[[input$yDatasetpb]][["featurePrefixes"]])), "Non valid data type"))
    ##
    
    yPrefix <- input$yPrefixpb
    if (!is.character(yPrefix)){
      yPrefix <- srcContentReactive()[[input$yDatasetpb]][["defaultFeatureY"]]
    }
    #
    originalId <- trimws(input$yIdpb)
    # 
    
    yId <- getMatchedIds(yPrefix, trimws(input$yIdpb), input$yDatasetpb, srcContent = srcContentReactive())
    
    if (length(yId) == 0){
      shiny::validate(need(FALSE, paste("ERROR:", paste0(yPrefix, input$yIdpb), "not found. Please use the Search IDs tab to find available IDs for each dataset.")))
    } else{
      globalReactiveValues$yPrefixpb <- yPrefix
      if (length(yId) > 1){
        warningMsg <- paste0("Other identifiers matching y-axis ID: ",
                             paste0(yId[-1], collapse = ", "), ".")
        showNotification(warningMsg, duration = 10, type = "message")
        yId <- yId[1]
      }
      # yData <- getFeatureData(yPrefix, yId, input$yDataset, srcContent = srcContentReactive())
      yData <- getFeatureData(yPrefix, yId, input$yDatasetpb, srcContent = srcContentReactive(), originalId)
      
      matchedLinesTab <- matchedCellLinesTabpb()
      # yData$data <- yData$data[matchedLinesTab[, "yDataset"]]
      yData$data <- yData$data[as.character(matchedLinesTab[, "yDatasetpb"])]
    }
    
    return(yData)
  })
  
  
  
  output$yAxisRangeUipb <- renderUI({
    srcContent <- srcContentReactive()
    
    # Note: see comment in output#xAxisRangeUi explaining the use of req().
    ## new
    shiny::validate(need(!is.na(match(input$yPrefixpb,srcContentReactive()[[input$yDatasetpb]][["featurePrefixes"]])), "Non valid data type"))
    
    valRange <- srcContent[[req(input$yDatasetpb)]][["featureValRanges"]][[req(input$yPrefixpb)]]
    
    yData <- NULL
    try(yData <- yDatapb())
    if (is.null(yData)){
      yInitSliderVals <- valRange
    } else{
      yDataRange <- range(yData$data, na.rm = TRUE)
      delta <- max((0.05 * (yDataRange[2] - yDataRange[1])), 0.1)
      yInitSliderVals <- c((yDataRange[1] - delta), (yDataRange[2] + delta))
    }
    
    ## 	sliderInput("yAxisRange", "y-Axis Range", min = valRange[1], max = valRange[2], value = yInitSliderVals, step = 0.5)
    sliderInput("yAxisRangepb", "Range", 
                min = valRange[1], max = valRange[2], value = yInitSliderVals, step = 0.5)
  })
  
  
  ypreviouspb <- debounce(reactive({
    if (is.null(input$yIdpb)) "YAP1" else input$yIdpb
  }), 5000)
  
  output$yIdUipb <- renderUI({
    
    ## textInput("xId", "Identifier: (e.g. ASCL1)", "ASCL1")
    
    if ( length(input$yPrefixpb)!=0 )  {
      if (input$yPrefixpb == "hll")
      {
        # srcContent <- srcContentReactive()
        cellChoices <- sort(srcContent[[input$yDatasetpb]][["molPharmData"]][["hllA"]]$ID)
        selectInput("yIdpb", "Please select a pathway",choices = cellChoices, selected = "APOPTOSIS")
      } else 
      {
        if (input$yPrefixpb == "mda") {
          cellChoices2 <- sort(srcContent[[input$yDatasetpb]][["molPharmData"]][["mdaA"]]$ID)
          selectInput("yIdpb", "Please select a miscellaneous variable",choices = cellChoices2, selected = "TMB")
        } else 
        {
          ##--------------------
          # 11/07 adc
          # if (input$yPrefixpb == "xsq" & length(input$yAdcpb) != 0)  {
          #   if (input$yAdcpb) {
          #     cellChoices3 <- sort(adcgenes$gene)
          #     selectInput("yIdpb", "Please select a gene",choices = cellChoices3, selected = "CD274")
          #   } else 
          #   {
          #     textInput("yIdpb", "Identifier: (e.g. YAP1)", ypreviouspb())
          #   }
          # }
         # else
          #{ 
          if (input$yPrefixpb == "sig") {
            cellChoices3 <- sort(srcContent[[input$yDatasetpb]][["molPharmData"]][["sigA"]]$ID)
            selectInput("yIdpb", "Please select a signature",choices = cellChoices3, selected = "APM_mean")
          }
          else {
            textInput("yIdpb", "Identifier: (e.g. YAP1)", ypreviouspb()) 
          }
          
          # }
        }
      }
    }
    else {
      
      textInput("yIdpb", "Identifier: (e.g. YAP1)", "YAP1")
      
    }
    
  })
  
  
  output$yAdcUipb <- renderUI({
    
    if ( length(input$yPrefixpb)!=0 )  {
      if (input$yPrefixpb == "xsq")
      {
        ## 11/07 ADC
        checkboxInput("yAdcpb", "ADC genes ?", value=FALSE)
        ## end ADC
      }
    }
    
  })
  
  output$yPrefixUipb <- renderUI({
    srcContent <- srcContentReactive()
    prefixChoices <- srcContent[[input$yDatasetpb]][["featurePrefixes"]]
    prefixChoices = sort(prefixChoices, decreasing = T) # Dec 20, 2024
    ### new
    ivar = match("var",prefixChoices)
    if (!is.na(ivar)) 	prefixChoices = prefixChoices[-ivar]
    ###
    selectedPrefix <- globalReactiveValues$yPrefixpb
    if ((is.null(selectedPrefix)) || (!(selectedPrefix %in% prefixChoices))){
      selectedPrefix <- srcContent[[input$yDatasetpb]][["defaultFeatureY"]]
    }
    opt = "";
    for(y in 1:length(prefixChoices)){
      if (prefixChoices[y]==selectedPrefix)
      {
        opt =  paste0(opt,"<option value=",prefixChoices[y]," selected>",names(prefixChoices)[y],"</option>;")
      }
      else
      {
        opt =  paste0(opt,"<option value=",prefixChoices[y],">",names(prefixChoices)[y],"</option>;");
      }
    }
    ## selectInput("yPrefix", "y-Axis Type", choices = prefixChoices, selected = selectedPrefix)
    selectInput("yPrefixpb", "Data Type", choices = prefixChoices, selected = selectedPrefix)
    
    # HTML(
    #   paste("<label class='control-label' for='yPrefix' id='lyp'>y-Axis Data Type</label>","<select id='yPrefix' style='word-wrap:break-word; width: 100%;'>",opt,"</select>")
    # )
  })
  
  analysisTissueTypespb <- reactive({
        
    tissueSelectionMode <- isolate(input$tissueSelectionModepb)
    
    
    selectedTissues = selectedTissuespb()
  
    
    srcContent <- srcContentReactive()
    tissueToSamplesMap <- srcContent[[input$xDatasetpb]][["tissueToSamplesMap"]]
    tissueTypes <- names(tissueToSamplesMap)
    
    if (tissueSelectionMode == "To include"){
      if (!("all" %in% selectedTissues)){
        selectedLines <- unique(c(tissueToSamplesMap[selectedTissues], recursive = TRUE))
        # For which tissue types are ALL lines in selectedLines?
        allInSelectedLines <- vapply(tissueToSamplesMap, function(x){
          all(x %in% selectedLines)
        }, logical(1))
        tissueTypes <- names(tissueToSamplesMap[allInSelectedLines])
      }
    } else{ # tissueSelectionMode == "Exclude"
      if (!("none" %in% selectedTissues)){
        selectedLines <- unique(c(tissueToSamplesMap[selectedTissues], recursive = TRUE))
        # For which tissue types are NO lines in selectedLines?
        notInSelectedLines <- vapply(tissueToSamplesMap, function(x){
          length(intersect(x, selectedLines)) == 0
        }, logical(1))
        tissueTypes <- names(tissueToSamplesMap[notInSelectedLines])
      }
    }
    
    #cat("--- LEAVING analysisTissueTypes()", sep = "\n")
    
    return(sort(unique(tissueTypes)))
  })
  
  
  matchedCellLinesTabpb <- reactive({
    srcContent <- srcContentReactive()
    analysisTissueTypes <- analysisTissueTypespb() ### to add
    
    if (input$xDatasetpb == input$yDatasetpb){
      matchedCellLinesTab <- data.frame(
        xDatasetpb = srcContent[[input$xDatasetpb]]$sampleData[, "Name"],
        stringsAsFactors = FALSE
      )
      matchedCellLinesTab$yDatasetpb <- matchedCellLinesTab$xDatasetpb
    } else{
      shiny::validate(need(require(rcellminerUtilsCDB),
                           "ERROR: x and y axis data sets must be the same."))
      matchedCellLinesTab <- getMatchedCellLines(c(input$xDatasetpb, input$yDatasetpb))
      shiny::validate(need(nrow(matchedCellLinesTab) > 0, 
                           "There are no shared cell lines between the selected datasets."))
      colnames(matchedCellLinesTab) <- c("xDatasetpb", "yDatasetpb")
    }
    stopifnot(all(!duplicated(matchedCellLinesTab$xDatasetpb)))
    rownames(matchedCellLinesTab) <- matchedCellLinesTab$xDatasetpb
    
    tissueMatchedLines <- getTissueTypeSamples(analysisTissueTypes, input$xDatasetpb, srcContent)
    tissueMatchedLines <- intersect(tissueMatchedLines, rownames(matchedCellLinesTab))
    
    shiny::validate(need(length(tissueMatchedLines) > 0, 
                         "There are no cell lines of the selected tissue type(s)."))
    
    matchedCellLinesTab <- matchedCellLinesTab[tissueMatchedLines, ]
    
    return(matchedCellLinesTab)
  })
  
  selectedTissuespb <- reactive({
    treepb <- input$treepb
    req(treepb)
    ## zz = get_selected(tree, format = "names")
    selectedTissuespb = names(unlist(get_selected(treepb, format = "slices")))
    selectedTissuespb = gsub("\\.",":", selectedTissuespb)
    selectedTissuespb = gsub("all:","", selectedTissuespb)
    selectedTissuespb = gsub("none:","", selectedTissuespb)
    cat(selectedTissuespb,"\n")
    return(sort(selectedTissuespb))
    
  })
  
  
  xDatapb <- reactive({
    ## shiny::validate(need(length(input$selectedTissues) > 0, "Please select tissue types."))
    shiny::validate(need(length(selectedTissuespb()) > 0, "Please select tissue types."))
    ## new
    cat(input$xPrefixpb," xprefix pb \n")
    shiny::validate(need(!is.na(match(input$xPrefixpb,srcContentReactive()[[input$xDatasetpb]][["featurePrefixes"]])), "Non valid data type"))
    ##
    xPrefix <- input$xPrefixpb
    if (!is.character(xPrefix)){
      xPrefix <- srcContentReactive()[[input$xDatasetpb]][["defaultFeatureX"]]
    }
    #
    originalId <- trimws(input$xIdpb)
    # 
    xId <- getMatchedIds(xPrefix, trimws(input$xIdpb), input$xDatasetpb, srcContent = srcContentReactive())
    
    if (length(xId) == 0){
      shiny::validate(need(FALSE, paste("ERROR:", paste0(xPrefix, input$xIdpb), "not found. Please use the Search IDs tab to find available IDs for each dataset.")))
    } else{
      globalReactiveValues$xPrefixpb <- xPrefix
      if (length(xId) > 1){
        warningMsg <- paste0("Other identifiers matching x-axis ID: ",
                             paste0(xId[-1], collapse = ", "), ".")
        showNotification(warningMsg, duration = 10, type = "message")
        xId <- xId[1]
      }
      # xData <- getFeatureData(xPrefix, xId, input$xDataset, srcContent = srcContentReactive())
      xDatapb <- getFeatureData(xPrefix, xId, input$xDatasetpb, srcContent = srcContentReactive(), originalId)
      
      matchedLinesTab <- matchedCellLinesTabpb()
      # xData$data <- xData$data[matchedLinesTab[, "xDataset"]]
      xDatapb$data <- xDatapb$data[as.character(matchedLinesTab[, "xDatasetpb"])]
      
    }
    
    return(xDatapb)
  })
  
  
  output$xPrefixUipb <- renderUI({
    srcContent <- srcContentReactive()
    
    # The last selected (data type) prefix is recorded in 
    # globalReactiveValues$xPrefix whenever xData() is updated. When the data set 
    # is changed, we try to use this same data type prefix, if it is available.
    prefixChoices <- srcContent[[input$xDatasetpb]][["featurePrefixes"]]  ## xsq
    
    ## selectedPrefix = srcContent[[input$xDatasetpb]][["defaultFeatureX"]]
    selectedPrefix= "xsq"
   ##  prefixChoices = sort(prefixChoices, decreasing = T) ## Dec20, 2024
    
    # ### new
    # ivar = match("var",prefixChoices)
    # if (!is.na(ivar)) 	prefixChoices = prefixChoices[-ivar]
    # ###
    # selectedPrefix <- globalReactiveValues$xPrefix
    # if ((is.null(selectedPrefix)) || (!(selectedPrefix %in% prefixChoices))){
    #   selectedPrefix <- srcContent[[input$xDataset]][["defaultFeatureX"]]
    #   if (is.na(selectedPrefix)) selectedPrefix <- srcContent[[input$xDataset]][["defaultFeatureY"]]
    # }
    opt = "";
    for(y in 1:length(prefixChoices)){
      if (prefixChoices[y]==selectedPrefix)
      {
        opt =  paste0(opt,"<option value=",prefixChoices[y]," selected>",names(prefixChoices)[y],"</option>;")
      }
      else
      {
        opt =  paste0(opt,"<option value=",prefixChoices[y],">",names(prefixChoices)[y],"</option>;");
      }
    }
    
    # selectInput("xPrefixpb", "Data Type", choices = prefixChoices, selected = selectedPrefix)
    

    # HTML(
    #   paste("<label class='control-label' for='xPrefixpb'>x-Axis Gene Expression</label>","<select hidden id='xPrefixpb' style='word-wrap:break-word; width: 100%;'>",opt,"</select>")
    # )
    HTML(
      paste("<label class='control-label' for='xPrefixpb'></label>","<select hidden id='xPrefixpb' style='word-wrap:break-word; width: 100%;'>",opt,"</select>")
    )
    
  })
  
 
  
  output$xIdUipb <- renderUI({
              # cellChoices3 <- sort(adcgenes$gene)
              cellChoices3 <- sort(unique(adcgenes$gene))
              selectInput("xIdpb", "Please select a predictive biomarker gene", choices = cellChoices3 , selected = "DLK1" ) 
              })
     
  output$xAxisRangeUipb <- renderUI({
    srcContent <- srcContentReactive()
    
    # Note: req() ensures values are available or 'truthy' (not NULL, "", FALSE, empty, etc.),
    # returning the value if so; otherwise the operation is stopped with a silent exception.
    # The idea is to exit quietly if inputs are momentarily in an invalid state, as might
    # occur when the app is first loading, etc.
    
    ## new
    shiny::validate(need(!is.na(match(input$xPrefixpb,srcContentReactive()[[input$xDataset]][["featurePrefixes"]])), "Non valid data type"))
    
    valRange <- srcContent[[req(input$xDatasetpb)]][["featureValRanges"]][[req(input$xPrefixpb)]]
    
    xData <- NULL
    try(xData <- xDatapb())
    if (is.null(xData)){
      xInitSliderVals <- valRange
    } else{
      xDataRange <- range(xData$data, na.rm = TRUE)
      delta <- max((0.05 * (xDataRange[2] - xDataRange[1])), 0.1)
      xInitSliderVals <- c((xDataRange[1] - delta), (xDataRange[2] + delta))
    }
    
    ##  	sliderInput("xAxisRange", "x-Axis Range", min = valRange[1], max = valRange[2], value = xInitSliderVals, step = 0.5)
    
    sliderInput("xAxisRangepb", "Range", 
                min = valRange[1], max = valRange[2], value = xInitSliderVals, step = 0.5)
  })
  
  # current_2var<-reactiveVal(NULL)
  # current_test<-reactiveVal(NULL)
  # 
  # observeEvent(input$nv,{
  #   cat(input$nv, "..nv\n")
  #   if(input$nv=="test"){
  #     cat("OK\n")
  #     ## updateNavbarPage(session,"2 variable Analyses")
  #     updateTabsetPanel(session,"nv", selected = "2 variable Analyses")
  #     # cat(input$nv, "..nv\n")
  #     current_test(1)
  #     current_2var(0) 
  #     cat("OK1 2var", current_2var(),"test ",current_test(), "\n")
  #   }
  #   else {
  #     current_2var(1)
  #     current_test(0)
  #     cat("OK2 2var", current_2var(),"test ",current_test(), "\n")
  #     
  #   }
  #   
  # })
          
  ## end---------------------------------------------------------------------------
})
#-----[end of shinyServer()]-----------------------------------------------------------------------
