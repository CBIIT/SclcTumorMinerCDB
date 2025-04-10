source("regressionModels.R")
source("myPatientModule.R")
metaConfig <- jsonlite::fromJSON("configMeta.json")
#--------------------------------------------------------------------------------------------------
# Helper functions.
#--------------------------------------------------------------------------------------------------
getMatchedIds <- function(prefix, id, dataSource, srcContent){
	dat <- srcContent[[dataSource]][["molPharmData"]][[prefix]]
	idSet <- rcellminer::removeMolDataType(rownames(dat))
	
	if(id %in% idSet){
		return(id)
	}
	
	matchedIds <- NULL
	
	i <- which(toupper(id) == toupper(idSet))
	if (length(i) == 1){
		# Straightforward case-insensitive match.
		matchedIds <- unname(idSet[i])
	} else{
		# For drugs: try to match synonyms to source-specific identifiers.
		#if (require(rcellminerUtilsCDB) && isDrugActivityDataType(prefix)){
			# matchedIds <- rcellminerUtilsCDB::getDbDrugIds(drugName = id, dbName = dataSource)
			# matchedIds <- intersect(matchedIds, idSet)
			# ## new final search with drug name
			# if (length(matchedIds)==0 & nchar(id)>=3) {
			#   vnames <- srcContent[[dataSource]][["drugInfo"]]
			#   matchedIds <- vnames[grep(toupper(id),toupper(vnames[,2])),1]
			# }
		#}
	  #
	  if (require(rcellminerUtilsCDB) && isGeneID(prefix)){
	    synonyms <- trimws(rcellminerUtilsCDB::getGeneSynonyms(id))
	    # cat(synonyms,"\n")
	    if (!is.null(synonyms) & length(synonyms)>1) {
	       matchedIndex <- which(!is.na(match(synonyms, idSet)))
         # cat(synonyms,":",matchedIndex,"\n")
	       if (length(matchedIndex)==1) {

	           matchedIds <- synonyms[matchedIndex]
	       }
	    }
	  }
	  #
	}
	
	return(matchedIds)
}

validateEntry <- function(prefix, id, dataSource, srcContent) {
	molPharmData <- srcContent[[dataSource]][["molPharmData"]]
	
	if(paste0(prefix, id) %in% rownames(molPharmData[[prefix]])) {
		return(TRUE)
	}
	
	return(FALSE)
}

getFeatureData <- function(prefix, id, dataSource, srcContent, originalId) {
	molPharmData <- srcContent[[dataSource]][["molPharmData"]]
	
	name <- paste0(prefix, id)
	data <- as.numeric(molPharmData[[prefix]][name, ])
	names(data) <- names(molPharmData[[prefix]][name, ])
	
	results <- list(name=name, data=data)
	
	# e.g., expTOP1 with dataSource=nci60 becomes TOP1 (exp, nci60)
	labs=metaConfig[[dataSource]][["displayName"]]
	if (toupper(id)==toupper(originalId))  {
	   ## results$plotLabel <- paste0(id, " (", prefix, ", ", labs, ")")
	   results$plotLabel <- paste0(id, " (", prefix, ")")
	}
	else {
	  ## if (prefix!="act") results$plotLabel <- paste0(toupper(originalId),": ",id, " (", prefix, ", ", labs, ")") else results$plotLabel <- paste0(id, " (", prefix, ", ", labs, ")")
	  if (prefix!="act") results$plotLabel <- paste0(toupper(originalId),": ",id, " (", prefix, ")") else results$plotLabel <- paste0(id, " (", prefix, ")")
	  
	}
	  
	
	# e.g., expTOP1 with dataSource=nci60 becomes expTOP1_nci60; needed for 
	# getPlotData() results (data.frame) with data for same feature from different sources.
	results$uniqName <- paste0(results$name, "_", dataSource)
	
	results$dataSource <- dataSource
	
	return(results)
}

getTissueTypeSamples <- function(tissueTypes, dataSource, srcContent) {
	matchedSamples <- c(lapply(tissueTypes, function(tissue){
		srcContent[[dataSource]]$tissueToSamplesMap[[tissue]]
	}), recursive=TRUE)
	return(unique(matchedSamples))
}

# Returns all tissue types associated with one or more samples in sampleSet.
getSampleSetTissueTypes <- function(sampleSet, dataSource, srcContent) {
	tissueToSamples <- srcContent[[dataSource]]$tissueToSamplesMap
	
	isMatchedType <- vapply(names(tissueToSamples), function(tissueType) {
		length(intersect(sampleSet, tissueToSamples[[tissueType]])) > 0
	}, logical(1))
	
	matchedTypes <- character(0)
	if (any(isMatchedType)) {
		## matchedTypes <- sort(unique(names(tissueToSamples[isMatchedType])))
		matchedTypes <-  unique(names(tissueToSamples[isMatchedType])) 
	}
	
	return(matchedTypes)
}

getPlotData <- function(xData, yData, showColor, showColorTissues, dataSource=NULL,
												srcContent){
	if (is.null(dataSource)){
		dataSource <- xData$dataSource
	}
	
	#-----[make sure x and y data cell lines are matched]----------------------------------
	if (xData$dataSource == yData$dataSource){
		stopifnot(identical(names(xData$data), names(yData$data)))
	}
	#--------------------------------------------------------------------------------------
	
	df <- data.frame(x=names(xData$data), y=xData$data, z=yData$data, stringsAsFactors = FALSE)
	rownames(df) <- df$x
	colnames(df) <- c("Cell Line", xData$uniqName, yData$uniqName)
	
	# HighCharts series name
	df$tissues <- srcContent[[dataSource]]$sampleData[rownames(df), "TissueType"]
	
	# HighCharts point name
	df$name <- srcContent[[dataSource]]$sampleData[rownames(df), "Name"]
	
	# HighCharts point x, y for scatter plot
	df$x <- df[,xData$uniqName]
	df$y <- df[,yData$uniqName]
	
	# RESTRICT (to rows with no NAs in either column x or column y).
	notNaData <- (!is.na(df[, xData$uniqName])) & (!is.na(df[, yData$uniqName]))
	df <- df[notNaData, ]
	
	if (nrow(df) > 0){
		cellLineSet <- rownames(df)
		# ADD COLOR COLUMN --------------------------------------------------------------------
		if (showColor){
			# Are there any cell lines to highlight in red?
			highlightedLineSet <- character(0)
			if (length(showColorTissues) > 0){
				highlightedLineSet <- unname(intersect(cellLineSet,
					getTissueTypeSamples(showColorTissues, dataSource, srcContent)))
			}
			#cat(highlightedLineSet,"22 \n")
			if (length(highlightedLineSet) > 0){
				## colorsToUse <- rep("rgba(0,0,255,0.3)", nrow(df)) #blue
				colorsToUse <- rep("rgba(192,192,192,0.3)", nrow(df)) # Grey
				
				names(colorsToUse) <- rownames(df)
				# colorsToUse[highlightedLineSet] <- "rgba(255,0,0,0.7)" # red
				   ### sel4=c("red","green","blue","orange")
				# 
				# sel7=c("red","green","blue","orange","cyan","yellow","violet")
				## sel7=c("red1","green1","blue1","brown","darkturquoise","khaki2","violet")
				# sel7=c("green1","blue1","brown","red","darkturquoise","khaki2","violet")
				sel7=c("green","brown","red","blue","darkturquoise","khaki2","violet")
				
				# sel7=rainbow(7)
				# library(RColorBrewer)
				# sel7 = brewer.pal(7, "Paired")
				
				sampleTissueTypes <- rep("unselected", nrow(df))
				names(sampleTissueTypes) <- rownames(df)
				for (k in 1:length(showColorTissues)) { 
				  ind=intersect(rownames(df),getTissueTypeSamples(showColorTissues[k], dataSource, srcContent))
				  ## sampleTissueTypes[ind]=showColorTissues[k] 
				  sampleTissueTypes[ind]=gsub(":","\n",showColorTissues[k])
				  ##
				    ### if (k >=1 & k<=4) colorsToUse[ind] <- sel4[k] else colorsToUse[ind] <-  "rgba(255,160,122,0.7)" # lightsalmon
				  if (k >=1 & k<=7) colorsToUse[ind] <- sel7[k] else colorsToUse[ind] <-  "rgba(255,160,122,0.7)" # lightsalmon
				  ##
				  }
			
				#cat(sampleTissueTypes,length(sampleTissueTypes),"\n")
				##---
			} else{
				# sampleTissueTypes <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree1"]
				# colorsToUse <- srcContent[[dataSource]]$tissueColorMap[sampleTissueTypes]
				
				sampleTissueTypes <- rep("unselected", nrow(df))
				## colorsToUse <- rep("rgba(0,0,255,0.3)", nrow(df)) #blue
				colorsToUse <- rep("rgba(192,192,192,0.3)", nrow(df)) # Grey
			}
			## ------
			df$sampleTissueTypes <- sampleTissueTypes
			## ------
		} else{
		  df$sampleTissueTypes <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree1"]
			## colorsToUse <- rep("rgba(0,0,255,0.3)", nrow(df)) #blue
			colorsToUse <- rep("rgba(192,192,192,0.3)", nrow(df)) # Grey
		}
		df$color <- colorsToUse
		
		# ADD ONCOTREE TISSUE TYPE COLUMNS ----------------------------------------------------
		df$OncoTree1 <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree1"]
		df$OncoTree2 <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree2"]
		df$OncoTree3 <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree3"]
		df$OncoTree4 <- srcContent[[dataSource]]$sampleData[rownames(df), "OncoTree4"]
		df$PlotTissueType <- ifelse(is.na(df$tissues),
																paste0(df$OncoTree1, ifelse(is.na(df$OncoTree2), "", 
																														paste0(":", df$OncoTree2))), 
																df$tissues)
		if (any(is.na(df$PlotTissueType))){
			df$PlotTissueType[which(is.na(df$PlotTissueType))] <- "TISSUE_TYPE_NA"
		}
		
		if ("EMT" %in% colnames(srcContent[[dataSource]]$sampleData)) {
			df$EMT <- srcContent[[dataSource]]$sampleData[rownames(df), "EMT"]
		}
		
		if ("sampleType" %in% colnames(srcContent[[dataSource]]$sampleData)) {
		  df$sampleType <- srcContent[[dataSource]]$sampleData[rownames(df), "sampleType"]
		}
		if ("biopsySite" %in% colnames(srcContent[[dataSource]]$sampleData)) {
		  df$biopsySite <- srcContent[[dataSource]]$sampleData[rownames(df), "biopsySite"]
		}
		
		if ("NAPY" %in% colnames(srcContent[[dataSource]]$sampleData)) {
		  df$NAPY <- srcContent[[dataSource]]$sampleData[rownames(df), "NAPY"]
		}
		
		if ("priorTreatment" %in% colnames(srcContent[[dataSource]]$sampleData)) {
		  df$priorTreatment <- srcContent[[dataSource]]$sampleData[rownames(df), "priorTreatment"]
		}
		if ("molecularSubtype" %in% colnames(srcContent[[dataSource]]$sampleData)) {
		  df$molecularSubtype <- srcContent[[dataSource]]$sampleData[rownames(df), "molecularSubtype"]
		}
		if ("dataSet" %in% colnames(srcContent[[dataSource]]$sampleData)) {
		  df$dataSet <- srcContent[[dataSource]]$sampleData[rownames(df), "dataSet"]
		}
		if ("patientID" %in% colnames(srcContent[[dataSource]]$sampleData)) {
		  df$patientID <- srcContent[[dataSource]]$sampleData[rownames(df), "patientID"]
		  kk = which(duplicated(df$patientID))
		  paired = which(df$patientID %in% unique(df$patientID[kk]))
		  df$paired_sample = "No"
		  df$paired_sample[paired] = "Yes"
		}
		
		
	}
	# cat(dim(df),"\n")
	return(df)
}


# makePlot <- function(xData, yData, showColor, showColorTissues, dataSource,
# 										 srcContent, dom="rCharts", showPValue = TRUE) {
# 	df <- getPlotData(xData, yData, showColor, showColorTissues, dataSource, srcContent)
# 	
# 	# Scatter plot
# 	h1 <- rCharts::Highcharts$new()
# 	
# 	# Divide the dataset, split by category and put into list() format
# 	# From: http://rcharts.io/viewer/?5735146#.VF6NS4W1Fy4
# 	series <- lapply(split(df, df$PlotTissueType), function(x) {
# 		res <- lapply(split(x, rownames(x)), as.list)
# 		names(res) <- NULL
# 		return(res)
# 	})
# 	
# 	invisible(sapply(series, function(x) {
# 		h1$series(data=x, type="scatter", name=x[[1]]$PlotTissueType)
# 	}
# 	))
# 	
# 	# Regression Line
# 	fit <- lm('y~x', data=df)
# 	x1 <- min(df[[xData$uniqName]])
# 	y1 <- fit$coefficients[[2]]*x1 + fit$coefficients[[1]]
# 	x2 <- max(df[[xData$uniqName]])
# 	y2 <- fit$coefficients[[2]]*x2 + fit$coefficients[[1]]
# 	
# 	h1$series(data=list(c(x1,y1), c(x2,y2)), type="line", color="#FF0000",
# 						marker=list(enabled=FALSE), enableMouseTracking=FALSE)
# 	
# 	# corResults <-cor.test(df[,xData$uniqName], df[,yData$uniqName], use="pairwise.complete.obs")
# 	# title <- paste0(paste(yData$plotLabel, '~', xData$plotLabel),
# 	# 								', r=', round(corResults$estimate, 2),
# 	# 								' p=', signif(corResults$p.value, 2))
# 	
# 	# rbind(y, y) forces more precise pvalue computation.
# 	corResults <- crossCors(df[,xData$uniqName], 
# 													rbind(df[,yData$uniqName], df[,yData$uniqName]))
# 	title <- paste0(paste(yData$plotLabel, '~', xData$plotLabel),
# 									', r=', signif(corResults$cor[1], digits=3))
# 	if (showPValue){
# 		title <- paste0(title, ' p=', signif(corResults$pval[1], digits=3))
# 	}
# 	
# 	h1$title(text=title)
# 	
# 	xAxisMin <- min(xData$data, na.rm = TRUE) - 0.25
# 	xAxisMax <- max(xData$data, na.rm = TRUE) + 0.25
# 	
# 	yAxisMin <- min(yData$data, na.rm = TRUE) - 0.25
# 	yAxisMax <- max(yData$data, na.rm = TRUE) + 0.25
# 	
# 	h1$xAxis(title=list(enabled=TRUE, text=xData$plotLabel, style=list(fontSize="24px", fontWeight="bold")),
# 					 min=xAxisMin, max=xAxisMax, labels=list(style=list(fontSize="20px")))
# 	h1$yAxis(title=list(enabled=TRUE, text=yData$plotLabel, style=list(fontSize="24px", fontWeight="bold")),
# 					 min=yAxisMin, max=yAxisMax, labels=list(style=list(fontSize="20px")))
# 	
# 	h1$legend(enabled=FALSE)
# 	
# 	# Force circle markers, set default size, hover color (otherwise color unpredictable)
# 	h1$plotOptions(series=list(animation=50),
# 								 scatter=list(marker=list(symbol='circle', radius=6,
# 								 												 states=list(hover=list(fillColor='white')))))
# 	
# 	tooltipFormat <- paste0("#! function() { return 'Cell: ' + this.point.name +
# 													'<br/>Tissue: ' + this.series.name +
# 													'<br/>", xData$uniqName, ": ' + Math.round(this.x * 100) / 100 +
# 													'<br/>", yData$uniqName, ": ' + Math.round(this.y * 100) / 100; } !#")
# 	
# 	h1$tooltip(backgroundColor="rgba(255,255,255,1)", formatter=tooltipFormat)
# 	
# 	h1$chart(zoomType="xy", style=list(fontFamily="Helvetica Neue"))
# 	
# 	# Enable exporting
# 	h1$exporting(enabled=TRUE)
# 	
# 	# Set name
# 	h1$set(dom=dom)
# 	
# 	# Print chart
# 	return(h1)
# }
# 

makePlotStatic <- function(xData, yData, showColor, showColorTissues, dataSource, 
													 srcContent, xLimVals = NULL, yLimVals = NULL,oncolor, showCells) {
	df <- getPlotData(xData, yData, showColor, showColorTissues, dataSource, srcContent)
	# contains column color
	shiny::validate(need(nrow(df)>0, paste("ERROR:", " No common complete data found.")))
	shiny::validate(need(nrow(df)>2, paste("ERROR:", " No display for less than 3 observations.")))

	df$biopsySite[which(is.na(df$biopsySite))] = "BiopsySite:NA"
	df$priorTreatment[which(is.na(df$priorTreatment))] = "PriorTreatment:NA"
	
		df$tooltip <- paste0(
		"Cell: ", df$name, "\n",
		"Tissue: ", df$PlotTissueType, "\n",
		"OncoTree1: ", df$OncoTree1, "\n",
		"OncoTree2: ", df$OncoTree2, "\n",
		df$sampleType, "\n",
		df$biopsySite, "\n",
		df$priorTreatment, "\n",
		df$dataSet, "\n",
		xData$uniqName, ": ", round(df$x, 2), "\n",
		yData$uniqName, ": ", round(df$y, 2), "\n"
	)
	
	# colorTab <- loadNciColorSet(returnDf=TRUE)
	# tissueColorTab <- unique(colorTab[, c("tissues", "colors")])
	tooltipCol <- "tooltip"
	
	# Plot parameters 
	# classCol <- "color"
	# colorPalette <- df[, "color"]
	# names(colorPalette) <- df[, classCol]
	# ##names(colorPalette) <- df[, "OncoTree1"]
	# # Merge data
	# df[, classCol] <- as.factor(df[, classCol])
#	if (showColor)
	   colorPalette <- oncolor[toupper(df[,"OncoTree1"]),]
	   classCol <- "OncoTree1"
	   leg <- TRUE
	# else 
	   if ( (length(showColorTissues) > 0) | (!showColor))
	         colorPalette <- df[, "color"]
	   
	   if (length(showColorTissues) > 0)  {
	     # classCol <- "color"
	     # leg <- FALSE
	     classCol <- "sampleTissueTypes"
	     colorPalette <- df[, "color"]
	     
	     	}
	   ## new -----------------------------
	   ## cat(showCells,  "1 \n")
	   if (length(showCells) > 0 & showColor)  {
	     # classCol <- "color"
	     # leg <- FALSE
	     write.csv(df,"df_showcells.csv")
	     df[, "sampleTissueTypes"] <- "unselected"
	     ## new stuff ---------------------------------
	     k = match(showCells, df$name)
	     patset = unique(df$patientID[k])
	     kpairs = which(df$patient %in% patset)
	     ipairs = setdiff(kpairs,k)
	     df[k, "sampleTissueTypes"] <- showCells
	     print(k)
	     print(patset)
	     print(kpairs)
	     print(ipairs)
	     if (length(ipairs) > 0 ) df[ipairs, "sampleTissueTypes"] <- df$name[ipairs]
	     ## end new stuff -----------------------------

	      ## old df[match(showCells, df$name), "sampleTissueTypes"] <- showCells
	     
	     # colorsToUse <- rep("rgba(0,0,255,0.3)", nrow(df)) #blue
	     # names(colorsToUse) <- rownames(df)
	     # colorsToUse[match(showCells, df$name)] <- "rgba(255,0,0,0.7)"
	     
	     ## df[, "color"] <- "rgba(0,0,255,0.3)" # blue
	     df[, "color"] <- "rgba(192,192,192,0.3)" # Grey
	     
	     df[match(showCells, df$name), "color"] <- "rgba(255,0,0,0.7)"
	     ## if (length(ipairs) > 0 ) df[ipairs, "color"] <- "rgba(0,128,0,0.7)"   ## green
	     if (length(ipairs) > 0 ) df[ipairs, "color"] <- "rgba(255,0,0,0.7)"   ## red
	     
	     classCol <- "sampleTissueTypes"
	     colorPalette <- df[, "color"]
	     
	   }
	   ## cat(df[,"color"], "2 \n")
	   ## end new -------------------------
	   
	df[, classCol] <- as.factor(df[, classCol])
	
	#colorPalette <- df[, "color"]
	names(colorPalette) <- df[, classCol]
	print(length(colorPalette))
  # write.csv(colorPalette,"colorpalette.csv")
	# write.csv(df,"df_showcells2.csv")
	# vtitle <- paste(yData$plotLabel, 'vs.', xData$plotLabel,"#samples:",nrow(df))

		# p1 <- rcellminer::plotCellMiner2D(df, xCol="x", yCol="y", xLabel = xData$plotLabel, yLabel = yData$plotLabel,
		# 										colorPalette=colorPalette, classCol=classCol, tooltipCol=tooltipCol,
		# 										xLimVal = xLimVals, yLimVal = yLimVals, showLegend = leg,pointSize = 4)
	 
		p1 <- plotCellMiner2Dv2(df, xCol="x", yCol="y", xLabel = xData$plotLabel, yLabel = yData$plotLabel,
		                                  colorPalette=colorPalette, classCol=classCol, tooltipCol=tooltipCol,
		                                  xLimVal = xLimVals, yLimVal = yLimVals, showLegend = leg,pointSize = 4)
		
		
	return(p1)
}

##-------------
## make correlation table by oncotype1 (or tissue of origin)
# library(dplyr)
# CorrelationTable <- function(xData, yData, showColor, showColorTissues=NULL, dataSource, 
#                            srcContent, xLimVals = NULL, yLimVals = NULL,oncolor) {
  
CorrelationTable_old <- function(xData, yData, srcContent) {
  
  df <- getPlotData(xData, yData, showColor=FALSE, showColorTissues=NULL, dataSource=NULL, srcContent)
  
  # for each oncotype 1 compute correlation
  # cor=df %>% group_by(OncoTree1) %>% summarize(cor.test(x,y,use="pairwise.complete.obs")$estimate) %>% as.data.frame
  # pval=df %>% group_by(OncoTree1) %>% summarize(cor.test(x,y,use="pairwise.complete.obs")$p.value) %>% as.data.frame
  # res=cbind(cor,pval[,2])
  # colnames(res)=c("Tissue_of_origin","Correlation","P-value")
  # 
  res=df %>% group_by(OncoTree1) %>% summarize(cor(x,y,use="pairwise.complete.obs")) %>% as.data.frame
  
  colnames(res)=c("Tissue_of_origin","Correlation")
  ## add counts 
  ## filter by number of rows / na should be more than 2 to do t.test
  return(res)
}

CorrelationTable <- function(xData, yData, srcContent) {
  
  df <- getPlotData(xData, yData, showColor=FALSE, showColorTissues=NULL, dataSource=NULL, srcContent)
  onco=unique(df$OncoTree1); onco=c("ALL",onco)
  nb=length(onco)
  res=matrix(NA,nb,3)
  rownames(res)=onco
  colnames(res)=c("Cell lines with complete observations","Correlation","P.value")
  for (k in 1:nb) {
    if (k==1) temp=df else temp=df[which(df$OncoTree1==onco[k]),]
    n=length(which(!is.na(temp$x) & !is.na(temp$x)))
    res[k,1]=n
    if (n>2)
    {
      tt=cor.test(temp$x,temp$y,use="pairwise.complete.obs")
      res[k,2]=tt$estimate
      res[k,3]=tt$p.value
    }
  }
  
  res=data.frame(res)
  res=cbind(rownames(res),res)
  res[, "P.value"] <- signif(res[, "P.value"], 2)
  res[, "Correlation"] <- signif(res[, "Correlation"], 2)
  res=res[order(res[,"P.value"]),]
  colnames(res)=c("Tissue of origin","Patients with complete observations","Pearson correlation","P-value")
  
  return(res)
}



##---------------



getLmEquationString <- function(predictorWts, orderByDecrAbsVal = TRUE, numSigDigits = 3){
	if (length(predictorWts) == 0){
		return("")
	}
	if (is.null(names(predictorWts))){
		names(predictorWts) <- paste0("predictor_", 1:length(predictorWts))
	}
	if (orderByDecrAbsVal){
		predictorWts <- predictorWts[order(abs(predictorWts), decreasing = TRUE)]
	}
	if ((length(predictorWts) > 0) && ("(Intercept)" %in% names(predictorWts))){
		i <- which(names(predictorWts) == "(Intercept)")
		predictorWts <- c(predictorWts[i], predictorWts[-i])
	}
	
	predictorWts <- signif(predictorWts, numSigDigits)
	if (names(predictorWts)[1] == "(Intercept)"){
		eqStr <- paste0("Y = ", predictorWts[1])
	} else{
		eqStr <- paste0("Y = (", predictorWts[1], "*", names(predictorWts)[1], ")")
	}
	
	if (length(predictorWts) > 1){
		for (predName in names(predictorWts[-1])){
			eqStr <- paste0(eqStr, " + (", predictorWts[predName], "*", predName, ")")
		}
	}
		
	return(eqStr)
}

#--------------------------------------------------------------------------------------------------
# searching drug IDs across all data sources using the synonyms list
# return a dataframe
## call : findDrugIDs("topo")
## find all : findDrugIDs("*")
#----------------------------------------------------------------------____________________________
findDrugIDs <- function(drugname) {
  tmp <- rcellminerUtilsCDB::drugSynonymTab
  y = unlist(lapply(tmp$NAME_SET, function(x) length(grep(drugname,x,ignore.case=T))))
  res=tmp[which(y!=0),]
  #found=which(y!=0)
  nb=dim(res)[2]
  nr=dim(res)[1]
  matres=matrix("",nr,nb)
  colnames(matres)=colnames(res)
  colnames(matres)[1]="Drug_Synonyms"
  for (i in 1:nb)
  {
    matres[,i]= unlist(lapply(res[,i], function(x) paste(x,collapse=";")))
  }
  dmatres=as.data.frame(matres)
  dmatres = dmatres[,c("Drug_Synonyms",intersect(names(metaConfig),colnames(matres)))]
  nbd=dim(dmatres)[2]
  # colnames(matres)[2:nb]=paste0(vapply(colnames(matres)[2:nb],function(x) {metaConfig[[x]][["displayName"]]},character(1)),"_IDs")
  # return(matres)
  colnames(dmatres)[2:nbd]=paste0(vapply(colnames(dmatres)[2:nbd],function(x) {metaConfig[[x]][["displayName"]]},character(1)),"_IDs")
  return(dmatres)
}
#---------------------------------------------------------------------------------------------------
# query TCGA data 
#---------------------------------------------------------

geneExpTcga <- function(varName, aproject) {
  #
  # First we're going to build the string representing the BigQuery #
  #
  # q <- paste(
  #   "SELECT project_short_name, sample_barcode, HGNC_gene_symbol, normalized_count,platform FROM `isb-cgc.TCGA_hg19_data_v0.RNAseq_Gene_Expression_UNC_RSEM` WHERE platform = 'IlluminaHiSeq' AND HGNC_gene_symbol = '", varName, "'",sep="")
 
  q <- paste(
    "SELECT project_short_name, Safe.SUBSTR(sample_barcode,14,2) AS ttype, HGNC_gene_symbol, normalized_count,platform FROM `isb-cgc.TCGA_hg19_data_v0.RNAseq_Gene_Expression_UNC_RSEM` WHERE HGNC_gene_symbol = '", varName, "'",sep="")
  
  ## new
  service_token <- set_service_token("./isb-cgc-fathi-b2c2573d2b53.json")
  ##
  response <- query_exec(q, project = aproject, use_legacy_sql =FALSE,max_pages = Inf)
  
  ##   dat <- data.frame()
  dat <- response
  if(!is.null(response))
  {
    
    # then we need to do a little data-cleaning to get ready for our survival model
    dat$project_short_name <- as.factor(dat$project_short_name)
    dat$normalized_count <- log2(dat$normalized_count+1)
  }
  return(dat)
} ## end expression query

### survival functions --------------------------------------------------------------

findInfoExp2old <- function(varName, cohort, geneexp, clin,sampletype,priortrt) {
  # find gene and survival information
  # sampleinfo is the oncotree3 
  
  if (sampletype) { 
    geneexp = geneexp[, which(clin$sampleType=="SampleType: primary tumor")] 
    clin = clin[which(clin$sampleType=="SampleType: primary tumor"),]
    # sampleinfo = sampleinfo[which(clin$sampleType=="SampleType: primary tumor")]
  }
  cat("1", dim(geneexp), " ", dim(clin), "\n")
  if (priortrt) {
    geneexp = geneexp[, which(clin$priorTreatment=="PriorTreatment: none")]
    clin    = clin[which(clin$priorTreatment=="PriorTreatment: none"),]
    # sampleinfo = sampleinfo[which(clin$priorTreatemnt=="PriorTreatment: none")]
  }
  cat("2", dim(geneexp), " ", dim(clin), "\n")
  dat =  data.frame()
  cat("varName: ", varName, "\n")
  if(!is.na(match(varName, rownames(geneexp))))
  { 
    ge <- geneexp[varName, which(clin$OncoTree3==cohort)]
    dat = data.frame(patient = clin[names(ge),"Name"], stringsAsFactors = F)
    ## dat = data.frame(patient = clin[names(ge),"Patient ID"], stringsAsFactors = F)
    dat$normalized_count <- ge
    # then we need to do a little data-cleaning to get ready for our survival model
    dat$days_to_last_known_alive <- as.numeric(as.character(clin[names(ge),"osDays"]))
    dat$vital_status <- clin[names(ge),"vitalStatus"]
    dat$vital_status <- ifelse(dat$vital_status == "alive", 0, 1)
    ### # dat$platform <- as.factor(dat$platform)
    # remove normal and controls and metastatic samples
    #vind = which(as.numeric(dat$ttype) %in% c(1,3,9))
    #
  }
  # return(dat[vind,])
  return(dat)
} 


drawPlotExp2old <- function(cohort, varname, geneexp, clin, sampletype,priortrt) {
  #
  # first make a call to BigQuery, and build the data frame
  ## dat <- findInfoExp(toupper(varname), cohort, geneexp,clin,sampleinfo)
  dat <- findInfoExp2(paste0("xsq",toupper(varname)), cohort, geneexp,clin,sampletype,priortrt)
  # write.csv(dat,"temp.csv")
  cat(dim(dat), "\n")
  shiny::validate(need(nrow(dat)>0, 
                       "data not found"))
  # unicity of patients? 
  # shiny::validate(need(!any(duplicated(dat$case_barcode)),
  #                      "duplicate patients found"))
  
  # if (any(duplicated(dat$case_barcode))) showNotification("Duplicated patients found!", duration=10, type="error")
  # if (length(levels(dat$platform))>1) showNotification("Platform is not unique", duration=10, type="warning")
  #
  # need to deal with duplicates / add cox model !!!!
  #
  #print(dat$normalized_count)
  response.M25=ifelse( dat$normalized_count <=quantile(dat$normalized_count)[2],0,NA); nb1 = length(which(response.M25==0))
  response.M25[which(dat$normalized_count>quantile(dat$normalized_count)[4])]=1; nb2 = length(which(response.M25==1))
  # response.M25=as.factor(response.M25)
  #
  # then we fit our survival model
  # fit <- survfit(Surv(days_to_last_known_alive, vital_status) ~ response.M25, data=dat)
  # lrt <- survdiff(Surv(days_to_last_known_alive, vital_status) ~ response.M25, data=dat)
  # plot(fit,main=paste("Survival prediction in ",cohort," for ",varname," gene expression",sep=""),lwd=2,xlab="Time (days)",ylab="Survival probability", col=c("blue","red"),lty=1)
  # mtext(paste("Log-rank P = ",signif(1-pchisq(lrt$chisq,1),4)," ",sep=""), side=1,line=-1,adj=0,cex=1)
  ## legend('topright', c("Q1","Q3"),lty=1, col=c("blue","red"), bty='n', cex=.75,title=varname)
  
  ## LSCC - slfn11 .00524
  ## now using cox model
  c0=coxph(Surv(days_to_last_known_alive, vital_status)~normalized_count,data=dat)
  c=coxph(Surv(days_to_last_known_alive, vital_status)~response.M25,data=dat)
  mfit.subtype1=survfit(c,newdata=data.frame(response.M25=c(0,1)))
  plot(mfit.subtype1,main=paste("Survival prediction in ",cohort," for ",toupper(varname)," gene expression using ",nrow(dat), " samples.",sep=""),lwd=2,xlab="Time (days)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5, col.main="blue", cex.axis=1, cex.lab=1.2)
  mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
  mtext(paste("HR 2groups: ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep=""),side=1,line=-2,adj=0,cex=1.2)

  # mtext(paste("HR 2groups: ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4)," ,Log-rank test p-val = ",round(summary(c)$logtest["pvalue"],4),sep=""),side=1,line=-2,adj=0,cex=1.2)
  
  mtext(paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep=""),side=1,line=-3,adj=0,cex=1.2)
  legend('topright', c(paste0("Exp <= 25% percentile: ",nb1," samples"),paste0("Exp  > 75%  percentile: ",nb2," samples")),lty=1, col=c("blue","red"), bty='n', cex=1,title=toupper(varname))
}

findInfoExp2 <- function(varName, cohort, geneexp, clin,sampletype,priortrt,optsurv) {
  # find gene and survival information
  # sampleinfo is the oncotree3 
  
  if (sampletype) { 
    geneexp = geneexp[, which(clin$sampleType=="SampleType: primary tumor")] 
    clin = clin[which(clin$sampleType=="SampleType: primary tumor"),]
    # sampleinfo = sampleinfo[which(clin$sampleType=="SampleType: primary tumor")]
  }
  # cat("1", dim(geneexp), " ", dim(clin), "\n")
  if (priortrt) {
    geneexp = geneexp[, which(clin$priorTreatment=="PriorTreatment: none")]
    clin    = clin[which(clin$priorTreatment=="PriorTreatment: none"),]
    # sampleinfo = sampleinfo[which(clin$priorTreatemnt=="PriorTreatment: none")]
  }
  # cat("2", dim(geneexp), " ", dim(clin), "\n")

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
  # predIds <- stringr::str_split(stringr::str_trim(input$predIds), pattern = "\\s+")[[1]] 
  nb=length(gl)
  comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp))
  if (length(comgenes)> 0) 
    ## 
    ## if(!is.na(match(toupper(varName), rownames(geneexp))))
  { 
    ##   ge <- geneexp[toupper(varName), which(clin$OncoTree3==cohort)]
    ge <- geneexp[comgenes, which(clin$OncoTree3==cohort), drop=F]
    dat = data.frame(patient = clin[colnames(ge),"Name"], stringsAsFactors = F)
    # dat = data.frame(patient = clin[names(ge),"Patient ID"], stringsAsFactors = F)
    ## dat$normalized_count <- ge
    dat$normalized_count <- colMeans(ge)
    # then we need to do a little data-cleaning to get ready for our survival model
    dat$days_to_last_known_alive <- as.numeric(as.character(clin[colnames(ge),"osDays"]))
    dat$vital_status <- clin[colnames(ge),"vitalStatus"]
    dat$vital_status <- ifelse(dat$vital_status == "alive", 0, 1)
    ### # dat$platform <- as.factor(dat$platform)
    # remove normal and controls and metastatic samples
    #vind = which(as.numeric(dat$ttype) %in% c(1,3,9))
    # cat("4 dat", dim(dat), "\n") 
  }
  # return(dat[vind,])
  return(list(data=dat,comgenes=comgenes))
} 


drawPlotExp2 <- function(cohort, varname, geneexp, clin, sampletype,priortrt,optsurv) {
  #
  # first make a call to BigQuery, and build the data frame
  ## dat <- findInfoExp(toupper(varname), cohort, geneexp,clin,sampleinfo)
  resu <- findInfoExp2(toupper(varname), cohort, geneexp,clin,sampletype,priortrt,optsurv)
  dat = resu$data
  comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
  # write.csv(dat,"temp.csv")
  # cat(dim(dat), "\n")
  shiny::validate(need(nrow(dat)>0, 
                       "data not found"))
  
  shiny::validate(need(length(na.omit(dat$normalized_count))>0, 
                       "all values are NA"))
  #
  #print(dat$normalized_count)
  tlabchange=F
  response.M25 = dat$normalized_count
  if (!(comgenes %in% c("CGA.STAIN","HORMONE.PROD","SYP.POSITIVE", "INSM1","MALEGENDER"))) {
  response.M25=ifelse( dat$normalized_count <=quantile(dat$normalized_count,na.rm=T)[2],0,NA)
  response.M25[which(dat$normalized_count>quantile(dat$normalized_count,na.rm=T)[4])]=1
  tlabchange = T
  } 
  


  nb1 = length(which(response.M25==0))
  nb2 = length(which(response.M25==1))
  ## now using cox model
  c0=coxph(Surv(days_to_last_known_alive, vital_status)~normalized_count,data=dat)
  c=coxph(Surv(days_to_last_known_alive, vital_status)~response.M25,data=dat)
  mfit.subtype1=survfit(c,newdata=data.frame(response.M25=c(0,1)))
  
  
  if (optsurv=="xsq") stitle = " gene expression using " else stitle = " metadata using "
  
  plot(mfit.subtype1,main=paste("Survival prediction in ",cohort," for ",paste(comgenes,collapse=","),stitle,nrow(dat), " samples.",sep=""),lwd=2,xlab="Time (days)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
  ## mtext(paste("Log-rank P = ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1)
  mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
  mtext(paste("HR 2groups: ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep=""),side=1,line=-2,adj=0,cex=1.2)
  mtext(paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep=""),side=1,line=-3,adj=0,cex=1.2)
  # legend('topright', c("Q1","Q4"),lty=1, col=c("blue","red"), bty='n', cex=.75,title=paste(comgenes,collapse=","))
  if (tlabchange==F) 
    {
      stitle21 = "No   : "
      stitle22 = "Yes  : "
        } else { 
      stitle21 = "Exp <= 25% percentile: "
      stitle22 = "Exp  > 75%  percentile: "
    }
    legend('topright', c(paste0(stitle21,nb1," samples"),paste0(stitle22,nb2," samples")),lty=1, col=c("blue","red"), bty='n', cex=1.2,title=paste(comgenes,collapse=","))
  
}

## new ver 3  -----------------------------------------------------------------------------------------------------------------------
findInfoExp3 <- function(varName, cohort, geneexp, clin,sampletype,priortrt,optsurv) {
  # find gene and survival information
  # sampleinfo is the oncotree3 
  
  if (sampletype) { 
    geneexp = geneexp[, which(clin$sampleType=="SampleType: primary tumor")] 
    clin = clin[which(clin$sampleType=="SampleType: primary tumor"),]
    # sampleinfo = sampleinfo[which(clin$sampleType=="SampleType: primary tumor")]
  }
  # cat("1", dim(geneexp), " ", dim(clin), "\n")
  if (priortrt) {
    geneexp = geneexp[, which(clin$priorTreatment=="PriorTreatment: none")]
    clin    = clin[which(clin$priorTreatment=="PriorTreatment: none"),]
    # sampleinfo = sampleinfo[which(clin$priorTreatemnt=="PriorTreatment: none")]
  }
  # cat("2", dim(geneexp), " ", dim(clin), "\n")
  
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
  
  # mutation1
  if (optsurv=="mut") geneexp[geneexp>0] = 1
  
  nb=length(gl)
  shiny::validate(need(nb<21, 
                       "Error: Max number of genes is 20"))
  comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp))
  if (length(comgenes)> 0) 
    ## 
    ## if(!is.na(match(toupper(varName), rownames(geneexp))))
  { 
    ##   ge <- geneexp[toupper(varName), which(clin$OncoTree3==cohort)]
    ge <- geneexp[comgenes, which(clin$OncoTree3==cohort), drop=F]
    # dat = data.frame(patient = clin[colnames(ge),"Name"], stringsAsFactors = F)
    
    dat = data.frame(t(ge))
    colnames(dat) = paste0(tolower(substr(colnames(dat),1,3)), substr(colnames(dat),4,nchar(colnames(dat))))
    comgenes = colnames(dat)
        ## dat$normalized_count <- colMeans(ge)
    # then we need to do a little data-cleaning to get ready for our survival model
    dat$days_to_last_known_alive <- as.numeric(as.character(clin[colnames(ge),"osDays"]))
    dat$vital_status <- clin[colnames(ge),"vitalStatus"]
    dat$vital_status <- ifelse(dat$vital_status == "alive", 0, 1)
   
    # cat("4 dat", dim(dat), class(dat), "\n") 
  }
  #
  shiny::validate(need(nrow(dat)>1, 
                       "data not found or only 1 sample"))
  
  shiny::validate(need(length(na.omit(dat[,-match(c("days_to_last_known_alive","vital_status"),colnames(dat))]))>0, 
                       "all values are NA"))
  #
  #print(dat$normalized_count)
  # first compute cox
  
  nbc = length(comgenes)
  
  if (nbc==1) {
    dat$variable = dat[,comgenes]
    cesp = NULL
  } else {
    cesp=coxph(Surv(days_to_last_known_alive, vital_status)~.,data=dat)
    # print(summary(cesp))
    beta = unlist(summary(cesp)[7])[1:nbc]
    rscore = as.matrix(dat[,comgenes]) %*% beta ## risk score
    dat$variable = rscore
  }
  
  
  #
  return(list(data=dat,comgenes=comgenes,coxmulti = cesp))
} 

drawPlotExp3 <- function(cohort, varname, geneexp, clin, sampletype,priortrt,optsurv) {
  #
  # first make a call to BigQuery, and build the data frame
  ## dat <- findInfoExp(toupper(varname), cohort, geneexp,clin,sampleinfo)
  resu <- findInfoExp3(toupper(varname), cohort, geneexp,clin,sampletype,priortrt,optsurv)
  dat = resu$data
  # covariates = resu$comgenes
  comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
  # write.csv(dat,"temp.csv")
  # cat(dim(dat), "\n")
  # shiny::validate(need(nrow(dat)>0, 
  #                      "data not found"))
  # 
  # shiny::validate(need(length(na.omit(dat[,-c("days_to_last_known_alive","vital_status")]))>0, 
  #                      "all values are NA"))
  # #
  # #print(dat$normalized_count)
  # # first compute cox
  # 
  # nbc = length(covariates)
  # 
  # if (nbc==1) {
  #   dat$variable = dat[,covariates]
  # } else {
  #   cesp=coxph(Surv(days_to_last_known_alive, vital_status)~.,data=dat)
  #   beta = unlist(summary(cesp)[7])[1:nbc]
  #   rscore = t(dat[,covariates]) %*% beta ## risk score
  #   dat$variable = rscore
  # }
  # 
  
  
  tlabchange=F
  response.M25 = dat$variable
  # mutation2
  if (!(comgenes %in% c("CGA.STAIN","HORMONE.PROD","SYP.POSITIVE", "INSM1","MALEGENDER")) & (optsurv!="mut") ) {
    response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[2],0,NA)
    response.M25[which(dat$variable>quantile(dat$variable,na.rm=T)[4])]=1
    tlabchange = T
  } 
  
  # cat(response.M25, "\n")
  
  
  nb1 = length(which(response.M25==0))
  nb2 = length(which(response.M25==1))
  ## now using cox model
  c0=coxph(Surv(days_to_last_known_alive, vital_status)~variable,data=dat)
  nbobs= unlist(summary(c0)[4])
  c=coxph(Surv(days_to_last_known_alive, vital_status)~response.M25,data=dat)
  
  # print(summary(c))
  # print(table(response.M25, dat$vital_status))
  # write.csv(cbind(dat,response.M25),paste0("tcga_acc_data_",comgenes,".csv"))
  
  mfit.subtype1=survfit(c,newdata=data.frame(response.M25=c(0,1)))
  
  
  if (optsurv=="xsq") {
    
    if (length(comgenes)==1) stitle = " gene expression using "  else stitle = " Risk score using "
    }
  else
  {
    if (optsurv=="mda") 
    {stitle = " metadata using " } else 
    {stitle = " mutation using " }
  }
  plot(mfit.subtype1,main=paste("Survival prediction in ",cohort," based on",stitle, paste(comgenes,collapse=",")," with ",nbobs, " samples.",sep=""),lwd=2,xlab="Time (days)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
  ## mtext(paste("Log-rank P = ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1)
  mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
  mtext(paste("HR 2groups: ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep=""),side=1,line=-2,adj=0,cex=1.2)
  mtext(paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep=""),side=1,line=-3,adj=0,cex=1.2)
  # legend('topright', c("Q1","Q4"),lty=1, col=c("blue","red"), bty='n', cex=.75,title=paste(comgenes,collapse=","))
  if (tlabchange==F) 
  {  if(optsurv!="mut") 
    {
    stitle21 = "No   : "
    stitle22 = "Yes  : " } 
    else 
    {
      stitle21 = "Wt   : "
      stitle22 = "Mut  : " }
  } else { 
    stitle21 = "Exp <= 25% percentile: "
    stitle22 = "Exp  > 75%  percentile: "
  }
  legend('topright', c(paste0(stitle21,nb1," samples"),paste0(stitle22,nb2," samples")),lty=1, col=c("blue","red"), bty='n', cex=1.2,title=paste(comgenes,collapse=","))
  
}

## new ver 4-------------------------------------------------------------------------------------------------------------------------
findInfoExp4 <- function(varName, cohort, geneexp, clin,sampletype,priortrt,optsurv,chkLasso,gsets,nbpred) {
  # find gene and survival information
  # sampleinfo is the oncotree3 >> dataSource
  
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
    shiny::validate(need(nbpred>=1 & nbpred <=20, "Error: The number of genes should be between 1 and 20"))
    shiny::validate(need(length(gsets) > 0, "Please select one or more gene sets."))
    ## cat(gsets, "...gsets..\n")
    # gl = union(gl, geneSetPathwayAnalysis::geneSets[[gsets]])
    gl = union(gl, unique(c(geneSetPathwayAnalysis::geneSets[gsets], recursive = TRUE)))
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
    beta = lassoPredictorWts0
    selgenes = names(lassoPredictorWts0)
    ## cat(selgenes, "selected genes .....\n")
    ## dat = data.frame(t(ge[selgenes,]))
    dat = data.frame(t(ge[selgenes,]))
    ## cat(dim(dat), "  dat  ..\n")
    rscore = t(ge[selgenes,]) %*% beta ## risk score
    
    dat$variable = rowMeans(t(ge[selgenes,]))
    dat$rscore = rscore
    
    ## cat(rscore, " ..rscore  ..\n")
    ## dat$variable = rscore
    dat$months_to_last_known_alive = months_to_last_known_alive
    dat$vital_status = vital_status
    ## cat(dim(dat), "  ..dat\n")
    ## cat(dim(res), "  ..cox multi\n")
    return(list(data=dat,comgenes=selgenes,coxmulti = res))
  }
  
  else
  {  
  
  # mutation1
  if (optsurv=="mut") geneexp[geneexp>0] = 1
  
  nb=length(gl)
  shiny::validate(need(nb<21, 
                       "Error: Max number of genes is 20"))
  comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp))
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
                       "data not found or only 1 sample"))
  
  shiny::validate(need(length(na.omit(dat[,-match(c("months_to_last_known_alive","vital_status"),colnames(dat))]))>0, 
                       "all values are NA"))
  #
  #print(dat$normalized_count)
  # first compute cox
  
  nbc = length(comgenes)
  
  if (nbc==1) {
    dat$variable = dat[,comgenes]
    cesp = NULL
  } else {
    cesp=coxph(Surv(months_to_last_known_alive, vital_status)~.,data=dat)
    # print(summary(cesp))
    beta = unlist(summary(cesp)[7])[1:nbc]
    rscore = as.matrix(dat[,comgenes]) %*% beta ## risk score
    
    dat$variable = rowMeans(as.matrix(dat[,comgenes]))
    dat$rscore = rscore
  }
  return(list(data=dat,comgenes=comgenes,coxmulti = cesp))
  # new for Lasso
  }
  
  #
 
} 

drawPlotExp4 <- function(cohort, varname, geneexp, clin, sampletype,priortrt,optsurv,chkLasso,gsets,nbpred) {
  #
  # first make a call to BigQuery, and build the data frame
  ## dat <- findInfoExp(toupper(varname), cohort, geneexp,clin,sampleinfo)
  resu <- findInfoExp4(toupper(varname), cohort, geneexp,clin,sampletype,priortrt,optsurv,chkLasso,gsets,nbpred)
  dat = resu$data
  # covariates = resu$comgenes
  comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq

    # write.csv(dat,"temp.csv")
  # cat(dim(dat), "\n")
  # shiny::validate(need(nrow(dat)>0, 
  #                      "data not found"))
  # 
  # shiny::validate(need(length(na.omit(dat[,-c("months_to_last_known_alive","vital_status")]))>0, 
  #                      "all values are NA"))
  # #
  # #print(dat$normalized_count)
  # # first compute cox
  # 
  # nbc = length(covariates)
  # 
  # if (nbc==1) {
  #   dat$variable = dat[,covariates]
  # } else {
  #   cesp=coxph(Surv(months_to_last_known_alive, vital_status)~.,data=dat)
  #   beta = unlist(summary(cesp)[7])[1:nbc]
  #   rscore = t(dat[,covariates]) %*% beta ## risk score
  #   dat$variable = rscore
  # }
  # 
  
  
  tlabchange=F
  response.M25 = dat$variable
  # mutation2
  if (!((comgenes %in% c("CGA.STAIN","HORMONE.PROD","SYP.POSITIVE", "INSM1","MALEGENDER"))[1]) & (optsurv!="mut") ) {
    response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[2],0,NA)
    response.M25[which(dat$variable>quantile(dat$variable,na.rm=T)[4])]=1
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
  c=coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)
  mfit.subtype1=survfit(c,newdata=data.frame(response.M25=c(0,1)))
  
  if (length(comgenes)<=10) genlab= paste(comgenes,collapse=",") else genlab = paste(paste(comgenes[1:10],collapse=","),paste(comgenes[11:length(comgenes)],collapse=","), sep="\n")
  labcohort = paste(cohort,collapse=",")
  
  if (optsurv=="xsq") {
    if (chkLasso) { stitle = " Lasso Avg genes of " } else {
    if (length(comgenes)==1) stitle = " gene expression of "  else stitle = " Multi-genes Average of " }
  }
  else
  {
    if (optsurv=="mda") 
    {stitle = " metadata using " } else 
    {stitle = " mutation of " }
  }
  ##plot(mfit.subtype1,main=paste("Survival prediction in ",cohort," based on \n",stitle, paste(comgenes,collapse=",")," with ",nbobs, " samples.",sep=""),lwd=2,xlab="Time (days)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
  plot(mfit.subtype1,main=paste("Survival prediction in ",labcohort," based on \n",stitle, genlab," with ",nbobs, " samples.",sep=""),lwd=2,xlab="Time (months)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
  
  ## mtext(paste("Log-rank P = ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1)
  
  # mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
  # mtext(paste("HR 2groups: ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep=""),side=1,line=-2,adj=0,cex=1.2)
  # mtext(paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep=""),side=1,line=-3,adj=0,cex=1.2)

#  mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
  mtext((paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep="")),side=1,line=-1,adj=0,cex=1.3,col="red")
  mtext((paste("HR 2 groups:    ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep="")),side=1,line=-2,adj=0,cex=1.3,col="red")
  mtext((paste("Hazard Ratio (HR): ")),side=1,line=-3,adj=0,cex=1.3,col="red")
  
  
    # legend('topright', c("Q1","Q4"),lty=1, col=c("blue","red"), bty='n', cex=.75,title=paste(comgenes,collapse=","))
  if (tlabchange==F) 
  {  if(optsurv!="mut") 
  {
    stitle21 = "No   : "
    stitle22 = "Yes  : " } 
    else 
    {
      stitle21 = "Wt   : "
      stitle22 = "Mut  : " }
  } else { 
    stitle21 = "Exp <= 25% percentile: "
    stitle22 = "Exp  > 75%  percentile: "
  }
  ## legend('topright', c(paste0(stitle21,nb1," samples"),paste0(stitle22,nb2," samples")),lty=1, col=c("blue","red"), bty='n', cex=1.2,title=paste(comgenes,collapse=","))
  legend('topright', c(paste0(stitle21,nb1," samples"),paste0(stitle22,nb2," samples")),lty=1, col=c("blue","red"), bty='n', cex=1.3,title=genlab , inset=0.03)
  
}

## new ver 5-------------------------------------------------------------------------------------------------------------------------
findInfoExp5 <- function(varName, cohort, geneexp, clin,sampletype,priortrt,optsurv,chkLasso,gsets,nbpred,gender) {
  # find gene and survival information
  # sampleinfo is the oncotree3 >> dataSource
  
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
    shiny::validate(need(nbpred>=1 & nbpred <=20, "Error: The number of genes should be between 1 and 20"))
    shiny::validate(need(length(gsets) > 0, "Please select one or more gene sets."))
    ## cat(gsets, "...gsets..\n")
    # gl = union(gl, geneSetPathwayAnalysis::geneSets[[gsets]])
    gl = union(gl, unique(c(geneSetPathwayAnalysis::geneSets[gsets], recursive = TRUE)))
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
    beta = lassoPredictorWts0
    selgenes = names(lassoPredictorWts0)
    ## cat(selgenes, "selected genes .....\n")
    ## dat = data.frame(t(ge[selgenes,]))
    dat = data.frame(t(ge[selgenes,]))
    ## cat(dim(dat), "  dat  ..\n")
    rscore = t(ge[selgenes,]) %*% beta ## risk score
    
    dat$variable = rowMeans(t(ge[selgenes,]))
    dat$rscore = rscore
    
    ## cat(rscore, " ..rscore  ..\n")
    ## dat$variable = rscore
    dat$months_to_last_known_alive = months_to_last_known_alive
    dat$vital_status = vital_status
    ## cat(dim(dat), "  ..dat\n")
    ## cat(dim(res), "  ..cox multi\n")
    return(list(data=dat,comgenes=selgenes,coxmulti = res))
  }
  
  else
  {  
    
    # mutation1
    if (optsurv=="mut") geneexp[geneexp>0] = 1
    
    nb=length(gl)
    shiny::validate(need(nb<21, 
                         "Error: Max number of genes is 20"))
    comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp))
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
                         "data not found or only 1 sample"))
    
    shiny::validate(need(length(na.omit(dat[,-match(c("months_to_last_known_alive","vital_status"),colnames(dat))]))>0, 
                         "all values are NA"))
    #
    #print(dat$normalized_count)
    # first compute cox
    
    nbc = length(comgenes)
    
    if (nbc==1) {
      dat$variable = dat[,comgenes]
      cesp = NULL
    } else {
      cesp=coxph(Surv(months_to_last_known_alive, vital_status)~.,data=dat)
      # print(summary(cesp))
      beta = unlist(summary(cesp)[7])[1:nbc]
      rscore = as.matrix(dat[,comgenes]) %*% beta ## risk score
      
      dat$variable = rowMeans(as.matrix(dat[,comgenes]))
      dat$rscore = rscore
    }
    return(list(data=dat,comgenes=comgenes,coxmulti = cesp))
    # new for Lasso
  }
  
  #
  
} 

drawPlotExp5 <- function(cohort, varname, geneexp, clin, sampletype,priortrt,optsurv,chkLasso,gsets,nbpred, gender) {
  #
  # first make a call to BigQuery, and build the data frame
  ## dat <- findInfoExp(toupper(varname), cohort, geneexp,clin,sampleinfo)
  resu <- findInfoExp5(toupper(varname), cohort, geneexp,clin,sampletype,priortrt,optsurv,chkLasso,gsets,nbpred,gender)
  dat = resu$data
  # covariates = resu$comgenes
  comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
  
  # write.csv(dat,"temp.csv")
  # cat(dim(dat), "\n")
  # shiny::validate(need(nrow(dat)>0, 
  #                      "data not found"))
  # 
  # shiny::validate(need(length(na.omit(dat[,-c("months_to_last_known_alive","vital_status")]))>0, 
  #                      "all values are NA"))
  # #
  # #print(dat$normalized_count)
  # # first compute cox
  # 
  # nbc = length(covariates)
  # 
  # if (nbc==1) {
  #   dat$variable = dat[,covariates]
  # } else {
  #   cesp=coxph(Surv(months_to_last_known_alive, vital_status)~.,data=dat)
  #   beta = unlist(summary(cesp)[7])[1:nbc]
  #   rscore = t(dat[,covariates]) %*% beta ## risk score
  #   dat$variable = rscore
  # }
  # 
  
  
  tlabchange=F
  response.M25 = dat$variable
  # mutation2
  if (!((comgenes %in% c("CGA.STAIN","HORMONE.PROD","SYP.POSITIVE", "INSM1","MALEGENDER"))[1]) & (optsurv!="mut") ) {
    response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[2],0,NA)
    response.M25[which(dat$variable>quantile(dat$variable,na.rm=T)[4])]=1
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
  c=coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)
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
    if (chkLasso) { stitle = " Lasso Avg genes of " } else {
      if (length(comgenes)==1) stitle = " gene expression of "  else stitle = " Multi-genes Average of " }
  }
  else
  {
    if (optsurv=="mda") 
    {stitle = " metadata using " } else 
    {stitle = " mutation of " }
  }
  
  if (gender=="Both") mesgender = "All" else mesgender = gender
  
  ##plot(mfit.subtype1,main=paste("Survival prediction in ",cohort," based on \n",stitle, paste(comgenes,collapse=",")," with ",nbobs, " samples.",sep=""),lwd=2,xlab="Time (days)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
  par(mar = c(5.1,4.1,7.1,2.1))
  plot(mfit.subtype1,main=paste("Survival prediction in ",labcohort," based on \n",stitle, genlab," with ",mesgender," ",nbobs, " patients.",sep=""),lwd=2,xlab="Time (months)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
  
  ## mtext(paste("Log-rank P = ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1)
  
  # mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
  # mtext(paste("HR 2groups: ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep=""),side=1,line=-2,adj=0,cex=1.2)
  # mtext(paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep=""),side=1,line=-3,adj=0,cex=1.2)
  
  #  mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
  mtext((paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep="")),side=1,line=-1,adj=0,cex=1.3,col="red")
  mtext((paste("HR 2 groups:    ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep="")),side=1,line=-2,adj=0,cex=1.3,col="red")
  mtext((paste("Hazard Ratio (HR): ")),side=1,line=-3,adj=0,cex=1.3,col="red")
  
  
  # legend('topright', c("Q1","Q4"),lty=1, col=c("blue","red"), bty='n', cex=.75,title=paste(comgenes,collapse=","))
  if (tlabchange==F) 
  {  if(optsurv!="mut") 
  {
    stitle21 = "No   : "
    stitle22 = "Yes  : " } 
    else 
    {
      stitle21 = "Wt   : "
      stitle22 = "Mut  : " }
  } else { 
    stitle21 = "Exp <= 25% percentile: "
    stitle22 = "Exp  > 75%  percentile: "
  }
  ## legend('topright', c(paste0(stitle21,nb1," patients"),paste0(stitle22,nb2," patients")),lty=1, col=c("blue","red"), bty='n', cex=1.2,title=paste(comgenes,collapse=","))
  legend('topright', c(paste0(stitle21,nb1," patients"),paste0(stitle22,nb2," patients")),lty=1, col=c("blue","red"), bty='n', cex=1.3,title=genlab , inset=0.03)
  
}


findInfoExpNew <- function(varName, cohort, geneexp, clin,sampletype,priortrt,optsurv,chkLasso,gsets,nbpred,gender) {
  # find gene and survival information
  # sampleinfo is the oncotree3 >> dataSource
  ## new: remove samples with all NA for all genes
  vna0 = apply(geneexp,2, function(x) length(which(is.na(x))))
  kna0 = which(vna0==nrow(geneexp))
  if (length(kna0)!=0) {
    geneexp = geneexp[,-kna0]
    clin = clin[-kna0,]
  }
  ## end new
  
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
    shiny::validate(need(nbpred>=1 & nbpred <=20, "Error: The number of genes should be between 1 and 20"))
    shiny::validate(need(length(gsets) > 0, "Please select one or more gene sets."))
    ## cat(gsets, "...gsets..\n")
    # gl = union(gl, geneSetPathwayAnalysis::geneSets[[gsets]])
    gl = union(gl, unique(c(geneSetPathwayAnalysis::geneSets[gsets], recursive = TRUE)))
    comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp)) ## adding prefix
    ## cat(length(gl), ".. \n")
    ## ge <- geneexp[comgenes, which(clin$OncoTree3==cohort), drop=F]
    # ge <- geneexp[comgenes, which(clin$OncoTree3 %in% cohort), drop=F]
    ge <- geneexp[comgenes, which(clin$dataSource %in% cohort), drop=F]
    rownames(ge) =  paste0(tolower(substr(rownames(ge),1,3)), substr(rownames(ge),4,nchar(rownames(ge))))
    ## cat(dim(ge), "..\n")
    
    # check if samples have NA for all genes so remove them for Lasso -------------- TO CHECK ??
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
    ## rm 1 res = data.frame(beta = lassoPredictorWts0, HR = exp(lassoPredictorWts0), stringsAsFactors = F)
    beta = lassoPredictorWts0
    selgenes = names(lassoPredictorWts0)
    res = data.frame(beta = lassoPredictorWts0, HR = exp(lassoPredictorWts0), index = 1:length(selgenes), gene = selgenes, stringsAsFactors = F)
    ## cat(selgenes, "selected genes .....\n")
    ## dat = data.frame(t(ge[selgenes,]))
    dat = data.frame(t(ge[selgenes,]))
    ## cat(dim(dat), "  dat  ..\n")
    rscore = t(ge[selgenes,]) %*% beta ## risk score
    
    dat$variable = rowMeans(t(ge[selgenes,]))
    dat$rscore = rscore
    
    ## cat(rscore, " ..rscore  ..\n")
    ## dat$variable = rscore
    dat$months_to_last_known_alive = months_to_last_known_alive
    dat$vital_status = vital_status
    ## cat(dim(dat), "  ..dat\n")
    ## cat(dim(res), "  ..cox multi\n")
    return(list(data=dat,comgenes=selgenes,coxmulti = res))
  }
  
  else ## NO LASSO
  {  
    
    # mutation1
    if (optsurv=="mut") geneexp[geneexp>0] = 1
    
    nb=length(gl)
    shiny::validate(need(nb<21, 
                         "Error: Max number of genes is 20"))
    comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp))
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
                         "data not found or only 1 sample"))
    
    shiny::validate(need(length(na.omit(dat[,-match(c("months_to_last_known_alive","vital_status"),colnames(dat))]))>0, 
                         "all values are NA"))
    #
    #print(dat$normalized_count)
    # first compute cox
    
    nbc = length(comgenes)
    
    if (nbc==1) {
      dat$variable = dat[,comgenes]
      cesp = NULL
    } else {
      cesp=coxph(Surv(months_to_last_known_alive, vital_status)~.,data=dat)
      # print(summary(cesp))
      beta = unlist(summary(cesp)[7])[1:nbc]
      rscore = as.matrix(dat[,comgenes]) %*% beta ## risk score
      
      dat$variable = rowMeans(as.matrix(dat[,comgenes]))
      dat$rscore = rscore
    }
    return(list(data=dat,comgenes=comgenes,coxmulti = cesp))
    # new for Lasso
  }
  
  #
  
} 

drawPlotExpNew <- function(cohort, varname, geneexp, clin, sampletype,priortrt,optsurv,chkLasso,gsets,nbpred, gender) {
  #
  # first make a call to BigQuery, and build the data frame
  ## dat <- findInfoExp(toupper(varname), cohort, geneexp,clin,sampleinfo)
  resu <- findInfoExpNew(toupper(varname), cohort, geneexp,clin,sampletype,priortrt,optsurv,chkLasso,gsets,nbpred,gender)
  dat = resu$data
  # covariates = resu$comgenes
  comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
  
  if (is.null(resu$coxmulti)) {
  ## one gene
  
  tlabchange=F
  response.M25 = dat$variable
  # mutation2
  if (!((comgenes %in% c("CGA.STAIN","HORMONE.PROD","SYP.POSITIVE", "INSM1","MALEGENDER"))[1]) & (optsurv!="mut") ) {
    response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[2],0,NA)
    response.M25[which(dat$variable>quantile(dat$variable,na.rm=T)[4])]=1
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
  c=coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)
  mfit.subtype1=survfit(c,newdata=data.frame(response.M25=c(0,1)))

  # survminer::ggsurvplot(fit=mfit.subtype1, data=dat, pval=T, risk.table=T, conf.int=T)   
 ##rm1  if (length(comgenes)<=10) genlab= paste(comgenes,collapse=",") else genlab = paste(paste(comgenes[1:10],collapse=","),paste(comgenes[11:length(comgenes)],collapse=","), sep="\n")
    genlab = comgenes
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
    # if (chkLasso) { stitle = " Lasso Avg genes of " } else {
    #   if (length(comgenes)==1) stitle = " gene expression of "  else stitle = " Multi-genes Average of " }
         stitle = " gene expression of "
  }
  else
  {
    if (optsurv=="mda") 
    {stitle = " metadata using " } else 
    {stitle = " mutation of " }
  }
  
  if (gender=="Both") mesgender = "All" else mesgender = gender
  
  ##plot(mfit.subtype1,main=paste("Survival prediction in ",cohort," based on \n",stitle, paste(comgenes,collapse=",")," with ",nbobs, " samples.",sep=""),lwd=2,xlab="Time (days)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
  par(mar = c(5.1,4.1,7.1,2.1))
  plot(mfit.subtype1,main=paste("Survival prediction in ",labcohort," based on \n",stitle, genlab," with ",mesgender," ",nbobs, " patients.",sep=""),lwd=3,xlab="Time (months)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=2,col.main="blue", cex.axis=1, cex.lab=1.2)
  
  ## mtext(paste("Log-rank P = ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1)
  
  # mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
  # mtext(paste("HR 2groups: ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep=""),side=1,line=-2,adj=0,cex=1.2)
  # mtext(paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep=""),side=1,line=-3,adj=0,cex=1.2)
  
  #  mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
  mtext((paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep="")),side=1,line=-1,adj=0,cex=1.5,col="red")
  mtext((paste("HR 2 groups:     ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep="")),side=1,line=-2,adj=0,cex=1.5,col="red")
  mtext((paste("Hazard Ratio (HR): ")),side=1,line=-3,adj=0,cex=1.5,col="red")
  
  
  # legend('topright', c("Q1","Q4"),lty=1, col=c("blue","red"), bty='n', cex=.75,title=paste(comgenes,collapse=","))
  if (tlabchange==F) 
  {  if(optsurv!="mut") 
  {
    stitle21 = "No   : "
    stitle22 = "Yes  : " } 
    else 
    {
      stitle21 = "Wt   : "
      stitle22 = "Mut  : " }
  } else { 
    stitle21 = "Exp <= 25% percentile: "
    stitle22 = "Exp  > 75%  percentile: "
  }
  ## legend('topright', c(paste0(stitle21,nb1," patients"),paste0(stitle22,nb2," patients")),lty=1, col=c("blue","red"), bty='n', cex=1.2,title=paste(comgenes,collapse=","))
  legend('topright', c(paste0(stitle21,nb1," patients"),paste0(stitle22,nb2," patients")),lty=1, col=c("blue","red"), bty='n', cex=1.3,title=genlab , inset=0.03)
  }
  else {
    if (chkLasso & optsurv=="xsq" ) ## Lasso HR
    {
      ggplot(data=resu$coxmulti, aes(y=index, x=HR)) +
        geom_point(shape=18, size=5) + 
        # geom_hline(aes(yintercept = index, colour = dff$colour), size=7) +
        scale_y_continuous(breaks=1:nrow(resu$coxmulti), labels=resu$coxmulti$gene) +
        labs(title='Lasso Results', x='Hazard Ratio (HR)', y = 'Gene') +
        geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
        ## theme(text = element_text(size=20)) +
        theme_classic(base_size = 20) 
      
      
    }
    else ## Cox multivariate HR
    {
      ## ggforest(resu$coxmulti, data= resu$dat)
      ggforest(resu$coxmulti, data= resu$dat,fontsize =1.8)
    }
  }
  
  
}

## new ver 6: cut perc , risk score ------------------------------------------------------------------------------------------------------------------------

# findInfoExp6 <- function(varName, cohort, geneexp, clin,sampletype,priortrt,optsurv,chkLasso,gsets,nbpred,gender) {
#   # find gene and survival information
#   # sampleinfo is the oncotree3 >> dataSource
#   
#   if (sampletype) { 
#     geneexp = geneexp[, which(clin$sampleType=="SampleType:Primary")] 
#     clin = clin[which(clin$sampleType=="SampleType:Primary"),]
#     # sampleinfo = sampleinfo[which(clin$sampleType=="SampleType: Primary")]
#   }
#   # cat("1", dim(geneexp), " ", dim(clin), "\n")
#   if (priortrt) {
#     geneexp = geneexp[, which(clin$priorTreatment=="PriorTreatment:none")]
#     clin    = clin[which(clin$priorTreatment=="PriorTreatment:none"),]
#     # sampleinfo = sampleinfo[which(clin$priorTreatemnt=="PriorTreatment: none")]
#   }
#   # cat("2", dim(geneexp), " ", dim(clin), "\n")
#   
#   if (gender == "Female") {
#     geneexp = geneexp[, which(clin$gender=="Gender:Female")]
#     clin    = clin[which(clin$gender=="Gender:Female"),]
#   } else
#   {
#     if (gender == "Male") {
#       geneexp = geneexp[, which(clin$gender=="Gender:Male")]
#       clin    = clin[which(clin$gender=="Gender:Male"),]
#     }
#   }
#   # cat("3", dim(geneexp), " ", dim(clin), "\n")
#   # Focus on non redundant samples
#   geneexp = geneexp[, which(clin$survivalStatus==TRUE)]
#   clin    = clin[which(clin$survivalStatus==TRUE),]
#   # cat("3", dim(geneexp), " ", dim(clin), "\n") 
#   #
#   dat =  data.frame()
#   ## new stuff
#   rownames(geneexp) = toupper(rownames(geneexp))
#   gl = trimws(toupper(varName))
#   if (optsurv=="xsq") gl = unlist(strsplit(gl,"\\s+"))
#   
#   if (chkLasso & optsurv=="xsq" ) {
#     shiny::validate(need(nbpred>=1 & nbpred <=20, "Error: The number of genes should be between 1 and 20"))
#     shiny::validate(need(length(gsets) > 0, "Please select one or more gene sets."))
#     ## cat(gsets, "...gsets..\n")
#     # gl = union(gl, geneSetPathwayAnalysis::geneSets[[gsets]])
#     gl = union(gl, unique(c(geneSetPathwayAnalysis::geneSets[gsets], recursive = TRUE)))
#     comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp)) ## adding prefix
#     ## cat(length(gl), ".. \n")
#     ## ge <- geneexp[comgenes, which(clin$OncoTree3==cohort), drop=F]
#     # ge <- geneexp[comgenes, which(clin$OncoTree3 %in% cohort), drop=F]
#     ge <- geneexp[comgenes, which(clin$dataSource %in% cohort), drop=F]
#     rownames(ge) =  paste0(tolower(substr(rownames(ge),1,3)), substr(rownames(ge),4,nchar(rownames(ge))))
#     ## cat(dim(ge), "..\n")
#     
#     # check if samples have NA for all genes so remove them for Lasso --------------
#     vna = apply(ge,2, function(x) length(which(is.na(x))))
#     kna = which(vna==nrow(ge))
#     if (length(kna)!=0) ge = ge[,-kna]
#     # ------------------------------------------------------------------------------
#     shiny::validate(need(ncol(ge)>=10, "Error: Min number of samples is 10"))
#     
#     months_to_last_known_alive <- as.numeric(as.character(clin[colnames(ge),"osMonths"]))
#     vital_status <- clin[colnames(ge),"vitalStatus"]
#     vital_status <- ifelse(vital_status == "alive", 0, 1)
#     ge = data.frame(ge)
#     
#     # length(vital_status)
#     ## cat(dim(ge), " dim ge..................\n") 
#     # ww = cbind(t(ge), time = months_to_last_known_alive, status = vital_status)
#     # write.csv(ww,"tempcox.csv")
#     # set.seed(1)
#     ## cat(months_to_last_known_alive, " days alive .....\n")
#     
#     # checking NA in status or time
#     jna = which(is.na(months_to_last_known_alive) | is.na(vital_status))
#     if (length(jna)!=0) {
#       ge = ge[,-jna]
#       months_to_last_known_alive = months_to_last_known_alive[-jna]
#       vital_status = vital_status[-jna]
#       shiny::validate(need(ncol(ge)>=10, "Error: Min number of samples is 10"))
#     }
#     
#     y = cbind(time = months_to_last_known_alive, status = vital_status)
#     #test
#     # write.csv(y,"y.csv")
#     #
#     set.seed(1) ## right one 
#     fit1_cv = cv.glmnet(t(ge), y, family = "cox", alpha = 1) # default alpha=1 >> Lasso
#     
#     # summary(fit1_cv)
#     # plot(fit1_cv)
#     # title("TCGA ACC Lasso, Cox Family", line = 2.5)
#     lassoPredictorWts0 <- coef(fit1_cv, s = "lambda.min")[,1] # deviance
#     lassoPredictorWts0 <- lassoPredictorWts0[lassoPredictorWts0 != 0]
#     lassoPredictorWts0 <- lassoPredictorWts0[order(abs(lassoPredictorWts0), decreasing = TRUE)]
#     
#     ## nb max of genes
#     if (length(lassoPredictorWts0) > nbpred) lassoPredictorWts0 = lassoPredictorWts0[1:nbpred]
#     shiny::validate(need(length(lassoPredictorWts0)>0, "No Lasso predictors found"))
#     ## cat(lassoPredictorWts0, "selected beta ......")
#     # based on lambda.min
#     res = data.frame(beta = lassoPredictorWts0, HR = exp(lassoPredictorWts0), stringsAsFactors = F)
#     beta = lassoPredictorWts0
#     selgenes = names(lassoPredictorWts0)
#     ## cat(selgenes, "selected genes .....\n")
#     ## dat = data.frame(t(ge[selgenes,]))
#     dat = data.frame(t(ge[selgenes,]))
#     ## cat(dim(dat), "  dat  ..\n")
#     rscore = t(ge[selgenes,]) %*% beta ## risk score
#     
#     dat$variable = rowMeans(t(ge[selgenes,]))
#     ## dat$variable = rscore
#     dat$rscore = rscore
#     
#     ## cat(rscore, " ..rscore  ..\n")
#     ## dat$variable = rscore
#     dat$months_to_last_known_alive = months_to_last_known_alive
#     dat$vital_status = vital_status
#     ## cat(dim(dat), "  ..dat\n")
#     ## cat(dim(res), "  ..cox multi\n")
#     return(list(data=dat,comgenes=selgenes,coxmulti = res))
#   }
#   
#   else
#   {  
#     
#     # mutation1
#     if (optsurv=="mut") geneexp[geneexp>0] = 1
#     
#     nb=length(gl)
#     shiny::validate(need(nb<21, 
#                          "Error: Max number of genes is 20"))
#     comgenes = intersect(paste0(toupper(optsurv),gl),rownames(geneexp))
#     if (length(comgenes)> 0) 
#       ## 
#       ## if(!is.na(match(toupper(varName), rownames(geneexp))))
#     { 
#       ##   ge <- geneexp[toupper(varName), which(clin$OncoTree3==cohort)]
#       ## new
#       ## ge <- geneexp[comgenes, which(clin$OncoTree3==cohort), drop=F]
#       # ge <- geneexp[comgenes, which(clin$OncoTree3 %in% cohort), drop=F]
#       ge <- geneexp[comgenes, which(clin$dataSource %in% cohort), drop=F]
#       # dat = data.frame(patient = clin[colnames(ge),"Name"], stringsAsFactors = F)
#       
#       dat = data.frame(t(ge))
#       colnames(dat) = paste0(tolower(substr(colnames(dat),1,3)), substr(colnames(dat),4,nchar(colnames(dat))))
#       comgenes = colnames(dat)
#       ## dat$normalized_count <- colMeans(ge)
#       # then we need to do a little data-cleaning to get ready for our survival model
#       dat$months_to_last_known_alive <- as.numeric(as.character(clin[colnames(ge),"osMonths"]))
#       dat$vital_status <- clin[colnames(ge),"vitalStatus"]
#       dat$vital_status <- ifelse(dat$vital_status == "alive", 0, 1)
#       
#       ## cat("4 dat", dim(dat), class(dat), "\n") 
#     }
#     # testing
#     # print(nrow(dat))
#     # write.csv(dat,"dat_test.csv")
#     # 
#     shiny::validate(need(nrow(na.omit(dat))>1, 
#                          "data not found or only 1 sample"))
#     
#     shiny::validate(need(length(na.omit(dat[,-match(c("months_to_last_known_alive","vital_status"),colnames(dat))]))>0, 
#                          "all values are NA"))
#     #
#     #print(dat$normalized_count)
#     # first compute cox
#     
#     nbc = length(comgenes)
#     
#     if (nbc==1) {
#       dat$variable = dat[,comgenes]
#       cesp = NULL
#     } else {
#       cesp=coxph(Surv(months_to_last_known_alive, vital_status)~.,data=dat)
#       # print(summary(cesp))
#       beta = unlist(summary(cesp)[7])[1:nbc]
#       rscore = as.matrix(dat[,comgenes]) %*% beta ## risk score
#       
#       dat$variable = rowMeans(as.matrix(dat[,comgenes]))
#       ## dat$variable = rscore
#       dat$rscore = rscore
#     }
#     return(list(data=dat,comgenes=comgenes,coxmulti = cesp))
#     # new for Lasso
#   }
#   
#   #
#   
# } 
# 
# drawPlotExp6 <- function(cohort, varname, geneexp, clin, sampletype,priortrt,optsurv,chkLasso,gsets,nbpred, gender) {
#   #
#   # first make a call to BigQuery, and build the data frame
#   ## dat <- findInfoExp(toupper(varname), cohort, geneexp,clin,sampleinfo)
#   resu <- findInfoExp6(toupper(varname), cohort, geneexp,clin,sampletype,priortrt,optsurv,chkLasso,gsets,nbpred,gender)
#   dat = resu$data
#   # covariates = resu$comgenes
#   comgenes = substr(resu$comgenes,4,nchar(resu$comgenes)) ## remove xsq
#   
#   # write.csv(dat,"temp.csv")
#   # cat(dim(dat), "\n")
#   # shiny::validate(need(nrow(dat)>0, 
#   #                      "data not found"))
#   # 
#   # shiny::validate(need(length(na.omit(dat[,-c("months_to_last_known_alive","vital_status")]))>0, 
#   #                      "all values are NA"))
#   # #
#   # #print(dat$normalized_count)
#   # # first compute cox
#   # 
#   # nbc = length(covariates)
#   # 
#   # if (nbc==1) {
#   #   dat$variable = dat[,covariates]
#   # } else {
#   #   cesp=coxph(Surv(months_to_last_known_alive, vital_status)~.,data=dat)
#   #   beta = unlist(summary(cesp)[7])[1:nbc]
#   #   rscore = t(dat[,covariates]) %*% beta ## risk score
#   #   dat$variable = rscore
#   # }
#   # 
#   
#   
#   tlabchange=F
#   response.M25 = dat$variable
#   # mutation2
#   if (!((comgenes %in% c("CGA.STAIN","HORMONE.PROD","SYP.POSITIVE", "INSM1","MALEGENDER"))[1]) & (optsurv!="mut") ) {
#     response.M25=ifelse( dat$variable <=quantile(dat$variable,na.rm=T)[2],0,NA)
#     response.M25[which(dat$variable>quantile(dat$variable,na.rm=T)[4])]=1
#     tlabchange = T
#   } 
#   
#   
#   
#   # nb1 = length(which(response.M25==0))
#   # nb2 = length(which(response.M25==1))
#   
#   # new counts
#   inb1 = which(response.M25==0)
#   inb2 = which(response.M25==1)
#   nb1 = length(which(!is.na(dat$months_to_last_known_alive[inb1])))
#   nb2 = length(which(!is.na(dat$months_to_last_known_alive[inb2])))
#   # 
#   ## now using cox model
#   c0=coxph(Surv(months_to_last_known_alive, vital_status)~variable,data=dat)
#   nbobs= unlist(summary(c0)[4])
#   c=coxph(Surv(months_to_last_known_alive, vital_status)~response.M25,data=dat)
#   mfit.subtype1=survfit(c,newdata=data.frame(response.M25=c(0,1)))
#   
#   if (length(comgenes)<=10) genlab= paste(comgenes,collapse=",") else genlab = paste(paste(comgenes[1:10],collapse=","),paste(comgenes[11:length(comgenes)],collapse=","), sep="\n")
#   labcohort = paste(cohort,collapse=",")
#   ## improve cohort display
#   ## integer div  % / % ,  mod %% >> rest 
#   if (length(cohort)>10) {
#     labcohort = ""
#     nbl = length(cohort) %/% 10
#     reste = length(cohort) %% 10
#     for (k in seq(1,(nbl*10), by=10)) {
#       labcohort = paste0(labcohort,paste(cohort[k:(k+9)],collapse=","),"\n")
#     }
#     if (reste!=0) labcohort = paste0(labcohort,paste(cohort[(k+10):(k+9+reste)],collapse=","))
#   }
#   
#   ## end improvments
#   
#   if (optsurv=="xsq") {
#     if (chkLasso) { stitle = " Average Expression of Lasso selected genes: " } else {
#       if (length(comgenes)==1) stitle = " gene expression of "  else stitle = " Average Expression of entered genes: " }
#   }
#   else
#   {
#     if (optsurv=="mda") 
#     {stitle = " metadata using " } else 
#     {stitle = " mutation of " }
#   }
#   
#   if (gender=="Both") mesgender = "All" else mesgender = gender
#   
#   ##plot(mfit.subtype1,main=paste("Survival prediction in ",cohort," based on \n",stitle, paste(comgenes,collapse=",")," with ",nbobs, " samples.",sep=""),lwd=2,xlab="Time (days)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
#   par(mar = c(5.1,4.1,7.1,2.1))
#   plot(mfit.subtype1,main=paste("Survival prediction in ",labcohort," based on \n",stitle, genlab," with ",mesgender," ",nbobs, " patients.",sep=""),lwd=2,xlab="Time (months)",ylab="Survival probability", col=c("blue","red"),lty=1,cex.main=1.5,col.main="blue", cex.axis=1, cex.lab=1.2)
#   
#   ## mtext(paste("Log-rank P = ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1)
#   
#   # mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
#   # mtext(paste("HR 2groups: ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep=""),side=1,line=-2,adj=0,cex=1.2)
#   # mtext(paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep=""),side=1,line=-3,adj=0,cex=1.2)
#   
#   #  mtext(paste("Log-rank test 2groups p-val: ",round(summary(c)$logtest["pvalue"],4)," ",sep=""), side=1,line=-1,adj=0,cex=1.2)
#   mtext((paste("HR continuous: ",round(summary(c0)$coefficients[2],4)," p-val: ",round(summary(c0)$coefficients[5],4),sep="")),side=1,line=-1,adj=0,cex=1.3,col="red")
#   mtext((paste("HR 2 groups:    ",round(summary(c)$coefficients[2],4)," p-val: ",round(summary(c)$coefficients[5],4),sep="")),side=1,line=-2,adj=0,cex=1.3,col="red")
#   mtext((paste("Hazard Ratio (HR): ")),side=1,line=-3,adj=0,cex=1.3,col="red")
#   
#   
#   # legend('topright', c("Q1","Q4"),lty=1, col=c("blue","red"), bty='n', cex=.75,title=paste(comgenes,collapse=","))
#   if (tlabchange==F) 
#   {  if(optsurv!="mut") 
#   {
#     stitle21 = "No   : "
#     stitle22 = "Yes  : " } 
#     else 
#     {
#       stitle21 = "Wt   : "
#       stitle22 = "Mut  : " }
#   } else { 
#     stitle21 = "Exp <= 25% percentile: "
#     stitle22 = "Exp  > 75%  percentile: "
#   }
#   ## legend('topright', c(paste0(stitle21,nb1," patients"),paste0(stitle22,nb2," patients")),lty=1, col=c("blue","red"), bty='n', cex=1.2,title=paste(comgenes,collapse=","))
#   legend('topright', c(paste0(stitle21,nb1," patients"),paste0(stitle22,nb2," patients")),lty=1, col=c("blue","red"), bty='n', cex=1.3,title=genlab , inset=0.03)
#   
# }

scaleDataForHeatmap <- function(dat, scaleByRow = FALSE){
  if (is.vector(dat)){
    vecNames <- names(dat)
    dat <- matrix(dat, nrow = 1, ncol = length(dat))
    rownames(dat) <- "tmp1"
    colnames(dat) <- vecNames
  }
  
  scaledDat <- matrix(NA, nrow = nrow(dat), ncol = ncol(dat))
  rownames(scaledDat) <- rownames(dat)
  colnames(scaledDat) <- colnames(dat)
  rowDataTypes <- unname(rcellminer::getMolDataType(rownames(dat)))
  dataTypes <- unique(rowDataTypes)
  
  getQuantiles <- function(dataType){
    if (dataType == "mut"){
      loQtl <- 0
      hiQtl <- 1
    } else{
      loQtl <- 0.05
      hiQtl <- 0.95
    }
    return(c(loQtl, hiQtl))
  }
  
  if (scaleByRow){
    for (i in seq_len(nrow(scaledDat))){
      valQtls <- quantile(x = as.numeric(dat[i, ]), 
                          probs = getQuantiles(rowDataTypes[i]), na.rm = TRUE)
      dTypeMin <- valQtls[1]
      dTypeMax <- valQtls[2]
      dTypeRange <- dTypeMax - dTypeMin
      if (dTypeRange != 0){
        tmp <- (dat[i, ] - dTypeMin) / dTypeRange
        tmp[which(tmp < 0)] <- 0
        tmp[which(tmp > 1)] <- 1
        scaledDat[i, ] <- tmp
      } else{
        scaledDat[i, ] <- 0.5
      }
    }
  } else{
    for (dType in dataTypes){
      indexSet <- which(rowDataTypes == dType)
      valQtls <- quantile(x = as.numeric(dat[indexSet, ]), 
                          probs = getQuantiles(dType), na.rm = TRUE)
      dTypeMin <- valQtls[1]
      dTypeMax <- valQtls[2]
      dTypeRange <- dTypeMax - dTypeMin
      for (i in indexSet){
        if (dTypeRange != 0){
          tmp <- (dat[i, ] - dTypeMin) / dTypeRange
          tmp[which(tmp < 0)] <- 0
          tmp[which(tmp > 1)] <- 1
          scaledDat[i, ] <- tmp
        } else{
          scaledDat[i, ] <- 0.5
        }
      }
    }
  }
  
  return(scaledDat)			 
}

## plot2D edited

plotCellMiner2Dv2 <- function(df, xCol="x", yCol="y", xLabel=xCol, yLabel=yCol, 
                            title=NULL, colorPalette=NULL, classCol=NULL, tooltipCol=NULL, 
                            showLegend=FALSE, showTrendLine=TRUE, showTitle=TRUE, singleColor="#0000FF",
                            alpha=1, numberColPrefix="X", xLimVal=NULL, yLimVal=NULL, pointSize=3) {
  
  # nci60DrugActZ <- exprs(getAct(rcellminerData::drugData))
  # nci60GeneExpZ <- getAllFeatureData(rcellminerData::molData)[["exp"]]
  # # Load colors
  # colorTab <- loadNciColorSet(returnDf=TRUE)
  # tissueColorTab <- unique(colorTab[, c("tissues", "colors")])
  # # Merge data
  # xCol <- "SLFN11"
  # yCol <- "94600"
  # classCol <- "tissues"
  # xLabel <- xCol
  # yLabel <- yCol
  # title <- NULL
  # df <- data.frame(y=nci60DrugActZ[yCol,], x=nci60GeneExpZ[xCol,])
  # yCol <- "X94600" # MUST NOT BE A NUMBER
  # colnames(df) <- c(yCol, xCol)
  # df <- cbind(df, colorTab)
  # df[, classCol] <- as.factor(df[, classCol])
  # colorPalette <- tissueColorTab[, "colors"]
  # names(colorPalette) <- tissueColorTab[, classCol]
  # colors <- rep("blue", nrow(df))
  # colors[1:10, "colors"] <- "red"
  # showLegend <- FALSE
  # showTrendLine <- TRUE
  # showTitle <- TRUE
  # alpha <- 1
  
  # Fix column names if they start with a number 
  if (grepl("^[0-9]", xCol)) {
    tmpXCol <- paste0(numberColPrefix, xCol)
    df[, tmpXCol] <- df[, xCol]
    xCol <- tmpXCol
  }
  
  if (grepl("^[0-9]", yCol)) {
    tmpYCol <- paste0(numberColPrefix, yCol)
    df[, tmpYCol] <- df[, yCol]
    yCol <- tmpYCol
  }
  
  # Create title 
  if(is.null(title)) {
    corResults <- cor.test(df[,xCol], df[,yCol], use="pairwise.complete.obs")
    title <- paste0(paste(yLabel, 'vs.', xLabel),
                    '\nPearson correlation (r)=', round(corResults$estimate, 2),
                    ', p-value=', formatC(signif(corResults$p.value, 2)),'\nNb.samples= ',nrow(df))
  }
  
  # Plot image
  p1 <- ggplot(data=df, aes_string(x=xCol, y=yCol))
  p1 <- p1 + theme_bw()
  
  if(!is.null(colorPalette) && !is.null(classCol)) {
    # p1 <- p1 + suppressWarnings(geom_point(aes_string(color=classCol, text=tooltipCol), alpha=alpha, size = pointSize))
    # p1 <- p1 + scale_colour_manual(name="", values=colorPalette)
    p1 <- p1 + suppressWarnings(geom_point(aes_string(fill = classCol, text=tooltipCol), alpha = alpha, size = pointSize, shape=21, stroke=0.15))
    p1 <- p1 + scale_fill_manual(name = "", values = colorPalette)
  } else {
    p1 <- p1 + suppressWarnings(geom_point(aes_string(text=tooltipCol), color=singleColor, alpha=alpha, size = pointSize))
  }
  
  if(!is.null(xLabel)) {
    p1 <- p1 + xlab(xLabel)		
  } else {
    p1 <- p1 + xlab(xCol)				
  }
  
  if(!is.null(yLabel)) {
    p1 <- p1 + ylab(yLabel)		
  } else {
    p1 <- p1 + ylab(yCol)				
  }
  
  if(showTrendLine) {
    p1 <- p1 + geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.25)
  }
  
  if(showTitle) {
    p1 <- p1 + ggtitle(title)
  }
  
  if(!showLegend) {
    p1 <- p1 + theme(legend.position="none", plot.title = element_text(size=14),
                     axis.title.x = element_text(size=20),
                     axis.title.y = element_text(size=20))
  }
  
  if(!is.null(xLimVal)) {
    p1 <- p1 + xlim(xLimVal[1], xLimVal[2])
  }
  
  if(!is.null(yLimVal)) {
    p1 <- p1 + ylim(yLimVal[1], yLimVal[2])
  }
  
  return(p1)
}


