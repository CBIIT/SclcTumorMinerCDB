# ORIGINAL CODE 
# gene_expression <- list(
#   "SLFN11",
#   "ASCL1",
#   "POU2F3",
#   "EZH2"
# )
# library(jsonlite)
# json_expression <- toJSON(gene_expression)
# 
# write(json_expression, "expression_biomarkers.json")

library(jsonlite)

#ACC Biomarkes (gene expression, mutation, methylation)

# acc_gene_expression <- list(
#   "TERT", "IGF2", "MEK", "BRAF", "CYP21A2", "CYP11B1", "HSD3B2", "MGARP", "STAR", "CYP17A1", "CCN3", "SULT2A1", 
#   "ATP4A", "CYP11A1", "FDXR", "SOAT1", "ABCB1", "SCARB1", "CYP11B2", "INHA", "HOXA5", "FDX1", "TM7SF2", "TBX3", 
#   "AOX1", "SLC16A9", "DLK1", "MC2R", "NR5A1", "MK167", "TOP2A", "MK167", "TOP2A", "TERF2")
# 
# acc_mutation <- list(
#   "TP53", "CTNNB1", "RPL22", "MEN1", "NF1", "RB1", "ZNRF3", "NF2")
# 
# acc_methylation <- list("ACCN4", "ACE", "ALX4", "CD14", "CDKN2A", "COL61A", "CRIP3", "DAPK1", "DMRT2", "FARP1", "G0S2", "GFRA3", "TAL1", "TNK2", "TRH")
# 
# #SCLC Biomarkers (gene expression, mutation, methylation)
# 
# sclc_gene_expression <- list (
#   "SRSF1", "SUV39H1", "GINS1", "PRPS1", "KPNA2", "AURKB", "TNPO2", "ORC6", "CCNA2", "LIG3", "MTF2", 
#   "GADD45G", "POLA1", "POLD4", "POLE4", "RFC5", "RMI1", "RRM1", "ABCC3", "AHNAK", "ANXA1", "ARHGDIB", 
#   "ASCL1", "BEX1", "BSN", "CAV1", "CAV2", "CCN1", "CCND1", "CELF3", "CHGA", "CHGB", "CRMP1", "EMP1", 
#   "EPHA2", "IFITM2", "IFITM3", "INSM1", "ITGB4", "KIF1A", "KIF5C", "LGALS3", "MYOF", "MYT1", "PLAU", 
#   "PTGES", "RAB27B", "RTN1", "RUNDC3A", "S100A10", "S100A16", "SCG3", "SEZ6", "SH3GL2", "SLC16A5", 
#   "SYN1", "SYP", "SYT11", "SYT4", "TACSTD2", "TAGLN3", "TFF3", "TGFBI", "TGFBR2", "TLCD3B", "TMSB15A", 
#   "TMSB15B:ENSG00000158427.15", "YAP1", "NEUROD1", "MYC", "MYCL", "MCYN", "KRAS", "FGFR1", "TP53", "RB1", "FHIT", 
#   "SLFN11", "POU2F3", "EZH2", "PD-L1", "YAP1", "NOTCH1", "REST", "HLA-A", "HLA-B", "HLA-C")
# 
# sclc_mutation <- list (
#   "TP53", "RB1", "RAD51D", "CHEK1", "BRAC2", "MLH1", "SMARCA4")
# 
# sclc_methylation <- list(
#   "SLFN11", "MGMT", "SMARCA1", "DNMT1", "DNMT3A", "DNMT3B"
# )
# 
# biomarkers <- list(
#   ACC = list(
#     gene_expression = acc_gene_expression,
#     mutation = acc_mutation,
#     methylation = acc_methylation
#   ),
#   SCLC = list(
#     gene_expression = sclc_gene_expression,
#     mutation = sclc_mutation,
#     methylation = sclc_methylation
#   )
# )
# 
# # Convert the nested structure to JSON
# json_biomarkers <- toJSON(biomarkers)
# 
# # Write the JSON to a file
# write(json_biomarkers, "biomarkers.json")





# Read the data from the CSV file
data <- data.frame(read.csv("TumorMinerBiomarkersListUpdated.csv"))

# Group the data by disease and data type and create the desired JSON structure
grouped_data <- split(data$Symbol, list(data$Disease, data$Data.Type))
disease_names <- unique(data$Disease)
 
# Initialize an empty list to store the final JSON data
json_data <- list()

# Loop through each disease
for (disease in disease_names) {
  json_data[[disease]] <- list()
  
  # Loop through each data type (gene_expression, mutation, methylation)
  for (data_type in c("gene_expression", "mutation", "methylation")) {
    symbols <- grouped_data[[paste(disease, data_type, sep = ".")]]
    json_data[[disease]][[data_type]] <- list(unlist(symbols))
  }
}

# Convert the nested list to JSON format
json_string <- toJSON(json_data, auto_unbox = TRUE, pretty = TRUE)

# Save the JSON to a file
write(json_string, file = "biomarkers.json")



#Convert json to excell to check

# library(jsonlite)
# library(writexl)
#
# json_data <- fromJSON("output.json")
# 
# data <- data.frame(Disease = character(),
#                    Data.Type = character(),
#                    Symbol = character(),
#                    stringsAsFactors = FALSE)
# 
# # Loop through the JSON data and append to the data frame
# for (disease in names(json_data)) {
#   for (data_type in names(json_data[[disease]])) {
#     symbols <- json_data[[disease]][[data_type]]
#     data <- rbind(data, data.frame(Disease = disease,
#                                    Data.Type = data_type,
#                                    Symbol = unlist(symbols),
#                                    stringsAsFactors = FALSE))
#   }
# }
# 
#
# write_xlsx(data, path = "output.xlsx")
