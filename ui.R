library(shiny)
library(rcellminer)
library(shinycssloaders)
library(plotly)
library(shinyTree)
library(shinyjs)

# options(repos = BiocInstaller::biocinstallRepos()) ## bioconductor 3.7
# getOption("repos")
#--------------------------------------------------------------------------------------------------
# LOAD CONFIGURATION AND REQUIRED DATA SOURCE PACKAGES.
#--------------------------------------------------------------------------------------------------
config <- jsonlite::fromJSON("config.json")
appConfig <- jsonlite::fromJSON("appConfig.json")
metaConfig <- jsonlite::fromJSON("configMeta.json")

toplinks <- appConfig$TopLinks
category <- appConfig$category
banner <- appConfig$banner

source("modal1.R")
source("appUtils.R")

if (!is.null(appConfig$appName)){
	appTitle <- appConfig$appName
} else{
	appTitle <- "CellMiner"
}

dataSourceChoices <- setNames(names(config),
															vapply(config, function(x) { x[["displayName"]] }, 
																		 character(1)))
options = "";
for(y in 1:length(dataSourceChoices)){
  if (dataSourceChoices[y]=="nci60")
  {
    options =  paste0(options,"<option value=",dataSourceChoices[y]," selected>",names(dataSourceChoices)[y],"</option>;")
  }
  else
   {
   options =  paste0(options,"<option value=",dataSourceChoices[y],">",names(dataSourceChoices)[y],"</option>;");
   }
  }

#print(options)
metaChoices <- setNames(names(metaConfig),
												vapply(metaConfig, function(x) { x[["displayName"]] }, 
															 character(1)))

metaoptions = "";
for(y in 1:length(metaChoices)){
  if (metaChoices[y]=="nci60")
  {
    metaoptions =  paste0(metaoptions,"<option value=",metaChoices[y]," selected>",names(metaChoices)[y],"</option>;")
  }
  else
  {
    metaoptions =  paste0(metaoptions,"<option value=",metaChoices[y],">",names(metaChoices)[y],"</option>;");
  }
}

listlinks = ''
for (k in 1:nrow(toplinks)) {
  listlinks=paste0(listlinks,tags$a(href=toplinks$url[k],toplinks$label[k],style="font-size: 18px;float: right;background-color: steelblue;color: white;display: inline-block;margin: 5px 5px;padding: 10px 10px;",target="_blank"),"\n")
}
# cat(listlinks)
if (category == "internal") mytitle="<p style='text-align: center; font-size: 20px; color:blue;' >~ Internal version ~</p>" else  
     if (category == "private") mytitle="<p style='text-align: center; font-size: 20px; color:red;' >~ Private version ~</p>" else 
          mytitle=""

#if("rCharts" %in% installed.packages()) {
#	options(RCHART_LIB='highcharts')	
#	library(rCharts)
#	hasRCharts <- TRUE
#} else {
#	hasRCharts <- FALSE
#}
## ---

shinyUI(
  function(request) {
  fluidPage(
   title=appTitle,
   ## title="SCLC TumorMinerCDB",
  tags$html(lang="en"), 
  ## tags$p(title="TumorMiner"), 
  ## tags$title("TumorMiner"),
  #tags$head(tags$style(type="text/css", ".body {color: blue;}",".clear {clear:both}")),
## rm1  tags$a(href="#skiplink","Skip over navigation",style="font-size: 10px; float: left"),
  
  # HTML("<p style='text-align: center; font-size: 20px; color:blue;' >~ Internal version ~</p>"),
## rm1  HTML(mytitle),
  #tags$h4("~Internal version~",style="color: blue"),
  # br(),
  # tags$html("~Internal version~",style="text-align: center; font-size: 20px"),
 
  # tags$a(href="https://discover.nci.nih.gov/cellminer/"," CellMiner NCI-60 ",style="font-size: 14px;float: right;background-color: steelblue;color: white;display: inline-block;margin: 5px 5px;padding: 10px 10px;",target="_blank"),
  # tags$a(href="https://dtp.cancer.gov"," NCI/DCTD/DTP ",style="font-size: 14px;float: right;background-color: steelblue;color: white;display: inline-block;margin: 5px 5px;padding: 10px 10px;",target="_blank"),
## rm1  HTML(listlinks),
  
  ###tags$p("CellMinerCDB",style="font-size: 24px;color: white;background-color: dodgerblue;text-align:center;height:50px;"),
  ### tags$img(src = "files/banner.jpg",height="110px",width="1650px"),
  # tags$img(src = "files/banner.png",alt= "banner",height="100%",width="100%", border="0"),
## rm1  tags$img(src = banner,alt= "banner",height="100%",width="100%", border="0"),
  
  #tags$img(src = "files/banner.png",alt= "banner",height="100%",width="100%", border="0", style="padding: 0px; display: block; line-height: 0; font-size: 0px; border: 0px; clear: both; vertical-align: top; margin: 0px 0px 0px 0px;"),
   #navbarPage(h6(style="vertical-align:top;font-size: 24px;color: dodgerblue;",appTitle), 
   # navbarPage(HTML("<p style='font-size: 24px;color: dodgerblue;'>", appTitle,"</p>"), 

## new 1 
## JMR new GUI ----------------------
## part1
includeHTML("www/uswds/ui/header.html"), ## new update JW***
## big div
tags$div(class="usa-section",
         tags$div(class="grid-container margin-bottom-10",
                  tags$div(class="gpf-content usa-prose site-prose",
                           tags$header(
                             ##  tags$p("Tools", class="site-subheading"),  # JW***
                             tags$h1(class="site-page-title tablet:margin-bottom-0",
                                     HTML(appConfig$appName))),
                           tags$div(class="margin-bottom-5",
                                    tags$span("Version:", class="post-date site-subheading",
                                              HTML(appConfig$appVersion)),
                                    tags$span("- Release:", class="post-date site-subheading",
                                              HTML(appConfig$appRelease))),
                           
## end part1        -----------------                 
                           
  	navbarPage(title="", id="nv",
						 inverse=FALSE,
						 header = list(tags$head(includeCSS("www/css/hacks.css")),
						 							 #tags$head(includeCSS("www/css/tooltip.css")),
						 							 # Add/run startup Javascript
						 							 
						 							 ## JMR part2 
						 							 tags$head(tags$link(rel="icon", type="image/png", href="uswds/img/favicons/favicon-32x32.png")),
						 							 tags$head(includeCSS("www/uswds/css/styles.css")),
						 							 tags$head(includeCSS("www/uswds/css/gpf_theme.css")),
						 							 ## end part2
						 							 
						 							 tags$head(tags$script(onloadJs)),
						 							 # Use JQuery (built into Shiny) calls to show/hide modal based on message
						 							 tags$head(includeScript("www/js/showLoading.js")),
						 							 tags$head(includeScript("www/js/showSkip.js")),
						 							 tags$head(includeScript("www/js/leaving.js")),
						 							 # load Javascript snippet to parse the query string.
						 							 #tags$script(includeScript("www/js/parse_input.js")),
						 							 tags$head(includeScript("www/js/google-analytics.js")),
						 							 
						 							 ## part3
						 							 tags$head(tags$script(src="uswds/js/uswds-init.min.js")), 
						 							 ## end part3
						 							 
						 							 tags$head(HTML("<script async type='text/javascript' src='https://dap.digitalgov.gov/Universal-Federated-Analytics-Min.js?agency=HHS&subagency=NCI' id='_fed_an_ua_tag'> </script>")),
						 							 tags$head(
						 							   tags$style(type="text/css", ".irs-grid-text { font-size: 8pt;color: black; }",
						 							              ".irs-min { font-size: 8pt; background: white; }", ".irs-max { font-size: 8pt; background: white;}",
						 							              ".irs-from { font-size: 8pt; color: black;background: white;}", ".irs-to { font-size: 8pt;  color: black;background: white;}"
						 							              , "body {font-size: 14pt;}", "img {display: block;}", ".clear {clear: both}"
						 							   )
						 							 ),
						 							 tags$head(
						 							   tags$style(HTML(
						 							     paste0(".navbar-nav { font-size: 24px; color: black; }"),
						 							     paste0(".navbar-default .navbar-brand { font-size: 24px; color: dodgerblue; }")
						 							   )
						 							   )
						 							 ),
						 							 # tags$head(tags$title("TumorMiner")),
						 							 tags$head(HTML("<title>TumorMiner</title>")),
						 							 tags$head(tags$meta(name="description",content="CellMiner Cross Database (CellMinerCDB) is the first web application to allow translational researchers to conduct analyses across all major cancer cell line pharmacogenomic data sources from NCI-DTP NCI-60, Sanger GDSC, and Broad CCLE/CTRP")),
						 							 tags$head(HTML("<script type=\"application/ld+json\">
  {
  \"@context\":\"https://schema.org/\",
  \"@type\":\"Dataset\",
  \"name\":\"NCI60 and other cancer cell line datasets\",
  \"description\":\" CellMinerCDB is a resource that simplifies access and exploration of cancer cell line pharmacogenomic data across different sources\",
  \"url\":\" https://discover.nci.nih.gov/cellminercdb/\",
  \"keywords\":[
  \"NCI60\",
  \"GDSC\",
  \"CCLE\",
  \"CTRP\"
  ],
  \"creator\":{
  \"@type\":\"Organization\",
  \"url\": \" https://discover.nci.nih.gov/ \",
  \"name\":\"GPF/DTB/CCR/NCI/NIH\",
  \"contactPoint\":{
  \"@type\":\"ContactPoint\",
  \"contactType\": \"customer service\",
  \"email\":\"Webadmin@discover.nih.gov\"
  }
  }
  }
  </script>"))
						              ),
		#background-color: blue; font-color: white;
		
		## new Mypatient Module -----------------
		## myPatientModuleUI("myPatient"),


### New navbarMenu ********************************************************************
### navbarMenu("Cross-Patients Analysis",
           
		#------[NavBar Tab: Univariate Analyses]---------------------------------------------------------
		## tabPanel("Univariate Analyses",
		tabPanel("2 variable Analyses",
			fluidPage(
    		loadingModal(),
	    	sidebarLayout(
	        sidebarPanel(
	          style = "height: 120vh; overflow-y: auto;", 
	        	width=3, 
	        	tags$div(
	        	  id="input_container",
	        	  tags$a(id="skiplink"),
	            #selectInput("xDataset", "x-Axis patient Set", choices=dataSourceChoices, selected = "nci60"),
	        	  HTML(
	        	    paste("<label class='control-label' for='xDataset'>Horizontal Axis</label>","<select hidden id='xDataset'>",options,"</select>")
	        	  ),  ## hidden for now 
	        	  uiOutput("xPrefixUi"),
	            ## textInput("xId", "Identifier: (e.g. TOP2A)", "TOP2A"),
	        	  ## 11/07 adc
	   ##     	  uiOutput("xAdcUi"),
	        	  ## end adc
	        	  uiOutput("xIdUi"),
	        	  conditionalPanel(condition="input.ts==1 || input.ts==2 || input.ts==4",
	        	  uiOutput("xAxisRangeUi") ),
	        	  br(),br(),
	            #selectInput("yDataset", "y-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
	        	  HTML(
	        	    paste("<label class='control-label' for='yDataset' id='lyd'>Vertical Axis</label>","<select hidden id='yDataset'>",options,"</select>")
	        	  ),
	        	  conditionalPanel(condition="input.ts==1 || input.ts==2 || input.ts==4",
	        	  uiOutput("yPrefixUi"),
	          	## textInput("yId", "Identifier: (e.g. MKI67)", "MKI67"),
	        	  ## 11/07 adc
	   ##     	  uiOutput("yAdcUi"),
	        	  ## end adc
	        	  uiOutput("yIdUi"),
	          	uiOutput("yAxisRangeUi"),
	          	
	            # checkboxInput("showColor", "Show Color?", value=TRUE),
              br()
	        	  ) # end conditional panel
	        	  , 
	          	## radioButtons("tissueSelectionMode", "Select Tissues", c("To include", "To exclude")),
	        	  radioButtons("tissueSelectionMode", "Select Subsets", c("To include", "To exclude")),
	        	  ## HTML("<p><b>Select Tissue/s of Origin</b></p>"),
	        	  HTML("<p><b>Select Subset(s)</b></p>"),
	        	  uiOutput("selectTissuesUi"),
	        	  
	        	  conditionalPanel(condition="input.ts==1 || input.ts==2 || input.ts==4",
	        	  checkboxInput("showColor", "Show Color?", value=TRUE),
	        	  ## HTML("<p><b>Select Tissue/s to color</b></p>"),
	        	  HTML("<p><b>Select Subset(s) to color</b></p>"),
	            uiOutput("showColorTissuesUi")
	        	  ) # end conditional panel
	        	  # br(),
	        	  # uiOutput("showCellsUi")
	        	)
	       
	        ),
        mainPanel(
         ## div(style="font-size: 16px", align="center", "CellMinerCDB enables exploration and analysis of cancer cell line pharmacogenomic data across different sources. If publishing results based on this site, please cite: ", a("Rajapakse.VN, Luna.A, Yamade.M et al. iScience, Cell Press. 2018 Dec 12.", href="https://www.cell.com/iscience/fulltext/S2589-0042(18)30219-0", target = "_blank", style="font-size: 16px;", class = "dm")),
         ## tags$head(tags$style(type='text/css', ".nav-tabs {font-size: 16px} ")),
         style = "font-size: 20px;",
         uiOutput('tabsetPanel')
        )
    	 )
			)
		),
		#-----[NavBar Tab: Regression Models]------------------------------------------------------------
		regressionModelsInput("rm", dataSourceChoices),
#-----[NavBar Tab: Expression Distributions]---------------------------------------------------------------------
		tabPanel("Expression distribution",
		         fluidPage(
		           sidebarLayout(
		             sidebarPanel(
		               width=3,
		               tags$div(
		                 id="input_container",
		                 tags$a(id="skiplink"),
		                 HTML(
		                   paste("<label class='control-label' for='dataDist'>Patient Set</label>","<select hidden id='dataDist'>",options,"</select>")
		                 ),
		                 br(),br(),br(),br(),
		                 textInput("geneDist", "Gene symbol: (e.g. SLFN11)", "SLFN11"),
		                 br(),
		                 # radioButtons("kmtype",label="Data type:", choices=list("Somatic mutation"=0,"Gene expression"=1)),
		                 checkboxInput("sampletype2", label="With only primary", value = F),
		                 checkboxInput("therapy2", label="Without prior therapy", value = F),
		                 br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br()
		               )
		             ), #end sidebarPanel
		             mainPanel(
		               
	#	               HTML("The following table displays the selected gene mutation variants for all patients with details such as variant allele frequency (VAF), mutation type (Exonic function) and Amino Acid changes (AAChange). The AAChange specfies the cDNA (c.) and the protein (p.) changes for all gene transcripts. All other details can be found by scrolling to the right of the screen."),
		               br(),br(),
		               withSpinner(plotlyOutput('genePlot', width=1400, height= 700)),
		               
		               p("The number on top of each boxplot represent its number of patients."),br(),
		               downloadButton("downloadPlotData", "Download Plot Data")
		               
		             )
		           )
		         ) #end fluidPage
		), #end tabPanel
		
#-----[NavBar Tab: Mutation variants]---------------------------------------------------------------------
tabPanel("Mutation variants",
         fluidPage(
           sidebarLayout(
             sidebarPanel(
               width=3,
               tags$div(
                 id="input_container",
                 tags$a(id="skiplink"),
                 #selectInput("mdataSource", "Data Source", choices=metaChoices, selected = "nci60")
                 HTML(
                   paste("<label class='control-label' for='dataMut'>Patient Set</label>","<select hidden id='dataMut'>",options,"</select>")
                 ),
                 br(),br(),br(),
                 textInput("geneId", "Gene Identifiers: (max 5, separated by a space e.g. TP53 RB1 MUC17 FLG2)", "TP53 RB1 MUC17 FLG2"),
                
               )
             ), #end sidebarPanel
             mainPanel(
               
               # HTML("The following table displays the selected gene mutation variants for all patients with details such as variant allele frequency (VAF), mutation type (Exonic function) and Amino Acid changes (AAChange). The AAChange specfies the cDNA (c.) and the protein (p.) changes for all gene transcripts. All other details can be found by scrolling to the right of the screen."),
               # br(),br(),
               ## uiOutput('mutationPanel')
               style = "font-size: 20px;",
               uiOutput('mutationPanelV2')
             )
           )
         ) #end fluidPage
) ,#end tabPanel

#-----[NavBar Tab: Survival]---------------------------------------------------------------------
## tabPanel("Survival Analyses", 
## tabPanel("test"),
tabPanel("Prognosis Biomarkers", 
         fluidPage(	
           sidebarLayout(
             sidebarPanel(
               width=3, 
               tags$div(
                 id="input_container", 
                 tags$a(id="skiplink"),
                 #selectInput("mdataSource", "Data Source", choices=metaChoices, selected = "nci60")
                 # HTML(
                 #   paste("<label class='control-label' for='survSource'>Patient Set</label>","<select id='survSource'>",metaoptions,"</select>")
                 # ),
                 HTML(
                   paste("<label class='control-label' for='survSource'>Patient Survival</label>","<select hidden id='survSource'>",metaoptions,"</select>")
                 ),
                 
                 br(),br(),
                 ## uiOutput("dataTypeUi"),
                 uiOutput("onco3Ui"),
                 br(),
                 
##                 radioButtons("optsurv","Data type", choices = c("Gene expression" = "xsq", "Gene mutation" = "mut","Metadata/signatures" = "mda" )),
                 radioButtons("optsurv","Data type", choices = c("Gene expression" = "xsq", "Gene mutation" = "mut","Metadata" = "mda" ,"Signatures" = "sig", "NMFs" = "nmf")),
                 
                 ### textInput("varname",  "Gene(s) Symbol(s) separated by space(s) or any metadata feature", value = "PINK1 BUB1B"),
                 # textInput("varname",  "Gene(s) Symbol(s) separated by space(s) or any metadata feature", value = "MKI67"),
                 uiOutput("varUi"),
                 # br(),
                 # radioButtons("kmtype",label="Data type:", choices=list("Somatic mutation"=0,"Gene expression"=1)),
                 checkboxInput("sampletype", label="With only primary", value = F),
                 checkboxInput("therapy", label="Without prior therapy", value = F),
                 ## br(),br(),br(),br(),br(),br(),br(),br(),br(),
                 ## radioButtons("optgroup","Groups selection by", choices = c("Q1-Q4 (lower quartile vs higher quartile)" = "25", "T1-T3 (lower tercile vs higher tercile)" = "33","Median" = "50" )),
                 radioButtons("optgroup","Groups selection by", choices = c("Median" = "50", "T1-T3 (lower tercile vs higher tercile)" = "33", "Q1-Q4 (lower quartile vs higher quartile)" = "25" )),
                 radioButtons("optgender","Gender", choices = c("Both" = "Both", "Female" = "Female","Male" = "Male" )),
                 ## checkboxInput("patientuniq", label="With unique patient sample", value = TRUE),
                 
                 ## new stuff --------------
                 ## br(),
                 conditionalPanel(condition="input['optsurv'] == 'xsq'",
                                  checkboxInput("Lasso","Run Lasso", value=FALSE),
                                  conditionalPanel(condition="input['Lasso'] == 1",
                                                   uiOutput("selectInputSurvLassoUI"),
                                                   numericInput("maxSurvPredictorsLasso", "Maximum Number of Genes", value = 4, 
                                                                min = 1, max = 10 , step = 1)
                                  ),
                                  conditionalPanel(condition="input['Lasso'] == 0",
                                                   br(),br(),br()     
                                  )
                                  
                 ),
                 conditionalPanel(condition="input['optsurv'] != 'xsq'",  br(),br(),br())
                 ## br(),br()
                 # conditionalPanel(condition="input['Lasso'] == 1",
                 #                  uiOutput("selectInputSurvLassoUI"),
                 #                  numericInput("maxSurvPredictorsLasso", "Maximum Number of Genes", value = 4, 
                 #                               min = 1, max = 30 , step = 1)
                 # )
                 
                 
                 # conditionalPanel(
                 #   # condition must be a Javascript expression.
                 #   condition = paste0("input['", ns("algorithm"), "'] == 'Lasso'"),
                 #   uiOutput(ns("selectInputGeneSetsUi"))),
                 # conditionalPanel(
                 #   condition = paste0("input['", ns("algorithm"), "'] == 'Lasso'"),
                 #   numericInput(ns("maxNumPredictors"), 
                 #                "Maximum Number of Predictors", value = 4, 
                 #                min = 1, max = 100, step = 1))
                 ## end new stuff
                 ##)
                 
               )
             ), #end sidebarPanel
             mainPanel(
               style = "font-size: 20px;",
               uiOutput('tabsetSurvival')
               # HTML("The survival plots are based on Cox proportional-hazards model. For now we can consider only one gene or one metadata. If you enter a set of genes, we compute the survival based on their average across all patient samples."),
               # br(),br(),
               # withSpinner(plotOutput("Survplot"))
               
             )
           )
         ) #end fluidPage
), #end tabPanel

### End navbarMenu ********************************************************************
### ),

##
		#-----[NavBar Tab: Metadata]---------------------------------------------------------------------
		## tabPanel("Metadata", 
		# tabPanel("Download data", 
		# 				 fluidPage(	
		# 				 	sidebarLayout(
		# 				 		sidebarPanel(
		# 				 			width=3, 
		# 				 			tags$div(
		# 				 				id="input_container", 
		# 				 				tags$a(id="skiplink"),
		# 				 				#selectInput("mdataSource", "Data Source", choices=metaChoices, selected = "nci60")
		# 				 				HTML(
		# 				 				  paste("<label class='control-label' for='mdataSource'>Patient Set</label>","<select hidden id='mdataSource'>",metaoptions,"</select>")
		# 				 				),
		# 				 			###	br(),br(),br(),br(),br(),br(),
		# 				 				br(),br(),
		# 				 				uiOutput("dataTypeUi"),
		# 				 				br(),
		# 				 				downloadButton('downloadExp', 'Download Data'),
		# 				 				br(),br(),
		# 				 				downloadButton('downloadFoot', 'Download Footnotes'),
		# 				 				br(),br(),br(),br(),br(),br(),
		# 				 				HTML("<b>Download current patient set information</b>"),
		# 				 				downloadButton('downloadCell', 'Download patients annotation'),
		# 				 				br(),br()
		# 				 				# HTML("<b>Download drug synonyms table with matching IDs for all patient sets</b>"),
		# 				 				# downloadButton('downloadSyn', 'Download Table'),
		# 				 				# br(),br()
		# 				 				#uiOutput(""),
		# 				 			)
		# 				 		), #end sidebarPanel
		# 				 		mainPanel(
		# 				 		  # htmlOutput('sourceLink'),
		# 				 		 ##  uiOutput('sourceLink'),
		# 				 			uiOutput('metadataPanel')
		# 				 			#h4(htmlOutput('sourceLink'))
		# 				 			# htmlOutput('sourceLink')
		# 				 		)
		# 				 	)
		# 				 ) #end fluidPage
		# ), #end tabPane 

### new predictive biomarkers , coming soon--------------------------------------------------------------------

tabPanel("Predictive biomarkers",
         fluidPage(
           loadingModal(),
           sidebarLayout(
             sidebarPanel(
               style = "height: 120vh; overflow-y: auto;", 
               width=3, 
               tags$div(
                 id="input_containerpb",
                 tags$a(id="skiplink"),
                 #selectInput("xDataset", "x-Axis patient Set", choices=dataSourceChoices, selected = "nci60"),
                 # HTML(
                 #   paste("<label class='control-label' for='xDatasetpb'>Horizontal Axis</label>","<select hidden id='xDatasetpb'>",options,"</select>")
                 # ),  ## hidden for now 
                 HTML(
                   paste("<label class='control-label' for='xDatasetpb'>Horizontal Axis</label>","<select hidden id='xDatasetpb'>",options,"</select>")
                 ),  ## hidden for now 
                 uiOutput("xPrefixUipb"),
                 ## textInput("xId", "Identifier: (e.g. TOP2A)", "TOP2A"),
                 
                 uiOutput("xIdUipb"),
                 conditionalPanel(condition="input.tspb==1 || input.tspb==2 || input.tspb==4",
                                  uiOutput("xAxisRangeUipb") ),
                 br(),br(),
                 #selectInput("yDataset", "y-Axis Dataset", choices=dataSourceChoices, selected = "nci60"),
                 HTML(
                   paste("<label class='control-label' for='yDatasetpb' id='lyd'>Vertical Axis</label>","<select hidden id='yDatasetpb'>",options,"</select>")
                 ),
                 conditionalPanel(condition="input.tspb==1 || input.tspb==2 || input.tspb==4", br(),
                                  uiOutput("yPrefixUipb"),
                                  ## textInput("yId", "Identifier: (e.g. MKI67)", "MKI67"),
                                  ## 11/07 adc
                                  # uiOutput("yAdcUipb"),
                                  ## end adc
                                  uiOutput("yIdUipb"),
                                  uiOutput("yAxisRangeUipb"),
                                  
                                  # checkboxInput("showColor", "Show Color?", value=TRUE),
                                  br()
                 ) # end conditional panel
                 , 
                 ## radioButtons("tissueSelectionMode", "Select Tissues", c("To include", "To exclude")),
                 radioButtons("tissueSelectionModepb", "Select Subsets", c("To include", "To exclude")),
                 ## HTML("<p><b>Select Tissue/s of Origin</b></p>"),
                 HTML("<p><b>Select Subset(s)</b></p>"),
                 uiOutput("selectTissuesUipb"),
                 
                 conditionalPanel(condition="input.tspb==1 || input.tspb==2 || input.tspb==4",
                                  checkboxInput("showColorpb", "Show Color?", value=TRUE),
                                  ## HTML("<p><b>Select Tissue/s to color</b></p>"),
                                  HTML("<p><b>Select Subset(s) to color</b></p>"),
                                  uiOutput("showColorTissuesUipb")
                 ) # end conditional panel
                 # br(),
                 # uiOutput("showCellsUi")
               )
               
             ),
             mainPanel(
               ## div(style="font-size: 16px", align="center", "CellMinerCDB enables exploration and analysis of cancer cell line pharmacogenomic data across different sources. If publishing results based on this site, please cite: ", a("Rajapakse.VN, Luna.A, Yamade.M et al. iScience, Cell Press. 2018 Dec 12.", href="https://www.cell.com/iscience/fulltext/S2589-0042(18)30219-0", target = "_blank", style="font-size: 16px;", class = "dm")),
               style = "font-size: 20px;",
               uiOutput('tabsetPanelpb')
             )
           )
         )
),


## new Mypatient Module -----------------

myPatientModuleUI("myPatient"),

#-----[NavBar Tab: Download]---------------------------------------------------------------------
tabPanel("Download data", 
         fluidPage(	
           HTML(paste("<label class='control-label' for='mdataSource'>Patient data download page</label>","<select hidden id='mdataSource'>",metaoptions,"</select>")
             				 				),
                 br(),br(),
                 uiOutput("dataTypeUi"),
                 br(),
                 downloadButton('downloadExp', 'Download Data'),
                 br(),br(),
                 HTML("<b>Download Footnotes</b>"),br(),
                 downloadButton('downloadFoot', 'Download Footnotes'),
                 br(),br(),
                 HTML("<b>Download RNASeq data without Batch correction</b>"),br(),
                 downloadButton('downloadExpNOBER', 'Download RNAseq Data'),
                 br(),br(),br(),br(),
                 HTML("<b>Download current patient set information</b>"),br(),
                 downloadButton('downloadCell', 'Download patients annotation'),
                 br(),br()
                 
               
         ) #end fluidPage
), #end tabPane 


		#-----[NavBar Tab: Search]---------------------------------------------------------------------
		tabPanel("Search IDs",
		         fluidPage(
		           sidebarLayout(
		             sidebarPanel(
		               width=3,
		               tags$div(
		                 id="input_container",
		                 tags$a(id="skiplink"),
		                 #selectInput("mdataSource", "Data Source", choices=metaChoices, selected = "nci60")
		                 HTML(
		                   paste("<label class='control-label' for='dataSrc'>patient Set</label>","<select hidden id='dataSrc'>",options,"</select>")
		                 ),
		                 br(),br(),
		                 uiOutput("dataTypeUi_s")
		                 # uiOutput("dataTypeUi"),
		                 # br(),
		                 # downloadButton('downloadExp', 'Download data for selected type')
		                 #uiOutput(""),
		               )
		             ), #end sidebarPanel
		             mainPanel(
		               #includeMarkdown("www/files/help.md"),
		             ##  DT::dataTableOutput("ids2")
		        ###       DT::dataTableOutput("ids_s")
	               uiOutput('searchPanel')
		               #h4(htmlOutput('sourceLink'))
	 #              htmlOutput('sourceLink')
		             )
		           )
		         ) #end fluidPage
		), #end tabPane
		#-----[NavBar Tab: About]------------------------------------------------------------------------
# 		tabPanel("TCGA distribution",
# 		            fluidPage(	
# 		              sidebarLayout(
# 		                sidebarPanel(
# 		                  width=3, 
# 		                  tags$div(
# 		                    id="input_container", 
# 		                    tags$a(id="skiplink"),
# 		                    #selectInput("mdataSource", "Data Source", choices=metaChoices, selected = "nci60")
# 		                    HTML(
# 		                      paste("<label class='control-label' for='mdataSource'>patient Set</label>","<select id='cmpSource'>",metaoptions,"</select>")
# 		                    ),
# 		                    br(),br(),
# 		                    textInput("vgene", "Gene symbol: ", "SLFN11"),
# 		                    #uiOutput("cmpTypeUi"),
# 		                    br(),br(),
# 		                    checkboxInput("zsid","z-score ",value=T),
# 		                    br(),br(),
# 		                    br(),br(),
# 		                    HTML("<b>Compare current patient set expression to TCGA</b>"),
# 		                    br(),
# 		                    actionButton('subtcga', 'Submit'),
# 		                    br(),br(),br(),br(),
# 		                    br(),br(),br(),br()
# 		                  )
# 		                ), #end sidebarPanel
# 		                mainPanel(
# 		                  #htmlOutput('sourceLink'),
# 		                  #uiOutput('sourceLink'),
# 		                  ## includeMarkdown("www/files/tcga.md"),
# 		                  withSpinner(plotlyOutput('tcgaPlot'))
# 		                  
# 		                )
# 		              )
# 		            ) #end fluidPage
#            ), 
		tabPanel("Help",
		         tags$a(id="skiplink"),
		         includeMarkdown("www/files/guide_sclcpat.md")
		         ## includeMarkdown("www/files/guide_pat.md")
		         ## includeHTML("www/files/guide2.html")
		         #h1("For testing"),
		         #textOutput("ipAddress")
		)
    #  tabPanel("Video tutorial",
    #       tags$a(id="skiplink"),
    #     includeMarkdown("www/files/video.md")
    # #      #h1("For testing"),
    # #      #textOutput("ipAddress")
    #  )
	),
br(),br(),hr(),
## rm 2
# tags$div(style="font-size: 12px",
# tags$p("CellMinerCDB is a development of the ",
# tags$a("Genomics and Pharmacology Facility,", href="https://discover.nci.nih.gov/", target = "_blank",style="font-size: 12px;"),
# tags$a(" Developmental Therapeutics Branch (DTB), ",href='https://ccr.cancer.gov/Developmental-Therapeutics-Branch', target='_blank',style="font-size: 12px;"),
# tags$a("Center for Cancer Research (CCR), ", href="https://ccr.cancer.gov/", target = "_blank",style="font-size: 12px;"),
# tags$a("National Cancer Institute (NCI) ", href="https://www.cancer.gov/", target = "_blank",style="font-size: 12px;"),
# "prepared in collaboration with the ",
# tags$a("cBio Center", href="http://www.sanderlab.org/", target = "_blank",style="font-size: 12px;", class = "dm"),
# " at the Dana-Farber Cancer Institute.",
# br(),br(),
# # tags$html("Please email 'Webadmin@discover.nci.nih.gov' with any problems, questions or feedback on the tool",style="font-size: 12px; float: left"),
# "Please ", 
# tags$a("email us", href="mailto:Webadmin@discover.nci.nih.gov&subject=CellMinerCDB",style="font-size: 12px;"),
# " with any problems, questions or feedback on the tool",
# br(),br(),
# tags$a("Notice and Disclaimer. ", href="files/disclaimer.html", target = "_blank"), 
# tags$a(" HHS Vulnerability Disclosure.",href='https://www.hhs.gov/vulnerability-disclosure-policy/index.html', target='_blank',style="font-size: 12px;"),
# 
# )
# )
                  ))), # end big div   
#JMR part 4
  includeHTML("www/uswds/ui/footer.html"),
  tags$div(includeScript("www/uswds/js/uswds.min.js"))

  )

  } # end function 
)
