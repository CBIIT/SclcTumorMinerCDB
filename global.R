# NOTE: Size is not automatically set for rChartsAlternative output
plotHeight <- 1100
plotWidth <- 1300

OutOfSurvival <- c("CSCLC_NCI","LNET_MyPart","NSCLC_NCI","SCCE_NCI","TNET_NCI","TNET_MyPart","TNET_NIDDK")

 srcContent_cell <- readRDS("srcContent_cell.rds")
# ccleinfo <- read.csv("ccle_cells_depmap_infosimilar.csv", row.names = 1, stringsAsFactors = F)

 ccleinfo <- read.csv("ccle_cells_depmap_infosimilar_drugs.csv", row.names = 1, stringsAsFactors = F)

## ADC target genes
## adcgenes = read.delim("Genes44_targetADC_v2.txt",row.names = 1)

# adcgenes0 = read.delim("Genes52_targetADC_drug_status.txt")
# adcgenes0 = adcgenes[, c(1,3,5,2,4,6)]

## adcgenes = read.delim("Genes53_targetADC_drug_041224.txt")
## adcgenes = read.delim("Genes60_targetADC_drug_112524.txt")
## adcgenes = read.delim("Genes72_targetADC_drug_06022025.txt")
adcgenes = read.delim("Genes63_targetADC_drug_06102025.txt")

adcgenes = adcgenes[, c(1,3,5,2,4,6)]


tooltipCol <- "tooltip"
## enableBookmarking("url")

isDrugActivityDataType <- function(prefix){
	# TO DO: Make configurable.
	drugActTypePrefixes <- "act"
	if (prefix %in% drugActTypePrefixes){
		return(TRUE)
	} else{
		return(FALSE)
	}
}

isGeneProtDataType <- function(prefix){
	# TO DO: Make configurable. this is for regression models
	geneProtDataTypePrefixes <- c("cop", "mut", "met", "exp", "xai","pro","swa", "xsq","mth","his","cri","rrb","tmm","tpm","muc")
	if (prefix %in% geneProtDataTypePrefixes){
		return(TRUE)
	} else{
		return(FALSE)
	}
}

isGeneID <- function(prefix){
  # TO DO: Make configurable.
  geneProtDataTypePrefixes <- c("cop", "mut", "met", "exp", "xai", "xsq","mth","his","cri","rrb","tmm","tpm","muc")
  if (prefix %in% geneProtDataTypePrefixes){
    return(TRUE)
  } else{
    return(FALSE)
  }
}