# SclcTumorMinerCDB

SclcTumorMinerCDB is a Shiny/R app that simplifies access and exploration of Small Cell Lung Cancer (SCLC) patient clinical and genomics data across different sources.

# Quickstart 
To run the app, run the following commands in R inside the project folder:
```
shiny::runApp()
```

# Dependencies and Source Repositories
## Dependencies Created by the NCI DTB Team 
* rcellminer: https://github.com/CBIIT/rcellminer
* rcellminerData: https://github.com/CBIIT/rcellminerData
* rcellminerUtilsCDB: https://github.com/CBIIT/rcellminerUtilsCDB
* rcellminerElasticNet: https://github.com/CBIIT/rcellminerElasticNet
* geneSetPathwayAnalysis: https://github.com/CBIIT/geneSetPathwayAnalysis
* tumorcomparer: https://github.com/sanderlab/tumorcomparer

## Other Dependencies 
* shiny
* shinyalert
* shinyjs
* shinycssloaders
* shinyTree
* base64enc
* bslib
* clusterProfiler
* DT
* htmltools
* maftools
* memoise
* markdown
* reshape2
* tidyr
* survival
* survminer
* heatmaply
* dplyr
* jsonlite
* stringr
* glmnet
* ggplot2
* gplots
* plotly
* foreach
* doSNOW
* parallel


## Zenodo data
The following pre-built data content (rds files) are available on Zenodo (https://zenodo.org/records/17793701) for use to run a local SclcTumorMinerCDB version. 

* srcContent.rds: SCLC patient data
* srcContent_cell.rds: SCLC cell line data from the Broad Institute (Cancer Cell Line Encycopedia) 

## Cache directory
You need to create in your local machine a cache directory and specify the path to this directory in the appConfig.json file (cacheDir field)
 
