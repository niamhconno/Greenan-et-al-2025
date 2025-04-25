######################################################################
### Pipeline to identify miRNA targets, filter, visualise, perform 
###  pathway enrichment analysis
### Performed in R version 4.3.3
### Code written by Niamh Connolly, Royal College of Surgeons in Ireland (RCSI)
### If using this code, please cite Greenan et al 2025:
### "Changes in Molecular Signature and Cytokine Expression on the Ocular 
### Surface in Response to Optimising Therapy in Aqueous-Deficient Dry Eye 
### Disease"
######################################################################

# Load required libraries
library(miRBaseConverter)
library(sqldf)
library(ggplot2)
library(httr)
library(tidyr)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)

# Set working directory 
setwd("C:/Users/niamhmconnolly/OneDrive - Royal College of Surgeons in Ireland/mirna-projects/emily-ocular/run2-92a-et-al") 

## Define miRNAs of interest and check their sequence and species conservation
mir_list <- as.data.frame(c("hsa-miR-4484","hsa-miR-6126","hsa-miR-7107-5p"))
colnames(mir_list) <- 'miRNA_ID'
head(mir_list)

#import mirbase file
mirbase22 <- getAllMiRNAs(version="v22", type="mature", species = )
colnames(mirbase22) <- c("Mature_Acc_22" , "Mature_ID_22" , "Mature_Seq_22")

# extract hsa+mmu+rno from mirbase
mirbase22 <- sqldf("SELECT * FROM mirbase22 WHERE 
                   mirbase22.Mature_ID_22 LIKE '%hsa%' OR
                   mirbase22.Mature_ID_22 LIKE '%mmu%' OR 
                   mirbase22.Mature_ID_22 LIKE '%rno%'")

# Extract miRNA ID and sequence for the hsa miRNAs of interest
hsa_mirna <- sqldf("SELECT mirbase22.Mature_ID_22, mirbase22.Mature_Seq_22
                   FROM mir_list, mirbase22
                   WHERE mirbase22.Mature_ID_22 LIKE '%hsa%'
                   AND mir_list.miRNA_ID = mirbase22.Mature_ID_22")

#combine tables and set column names
miRNAs_r_m = hsa_mirna
colnames(miRNAs_r_m) = c("miRNA_ID","miR_seq")

### Add miRNA family info, for information
version <- checkMiRNAVersion(miRNAs_r_m$miRNA_ID,verbose=FALSE)
# Get MIMAT ID
result <- miRNA_NameToAccession(miRNAs_r_m$miRNA_ID,version=version)
family <- checkMiRNAFamily(result$Accession)
miRNAs_r_m <- merge(miRNAs_r_m,family,by.x="miRNA_ID", by.y="miRNAName_v21")
miRNAs_r_m <- miRNAs_r_m[,c(1,2,5)]

head(miRNAs_r_m)

# remove unnecessary tables
rm(hsa_mirna, version, result, family)

######################################################################
## Identify predicted targets of our miRNA from miRDIP
### Code adapted from: http://ophid.utoronto.ca/mirDIP/api_R.jsp
### Citation: Tokar T, Pastrello C, Rossos AEM, Abovsky M, Hauschild AC, 
### Tsay M, Lu R, Jurisica I. mirDIP 4.1-integrative database of human 
### microRNA target predictions. Nucleic Acids Res. 2018 Jan 4;46(D1):D360-D370. 
### doi: 10.1093/nar/gkx1144. PubMed PMID: 29194489; PubMed Central 
######################################################################

url <- "http://ophid.utoronto.ca/mirDIP"
mapScore <- list("0", "1", "2", "3");
names(mapScore) <- c("Very High", "High", "Medium", "Low")

unidirectionalSearchOnMicroRNAs <- function(microRNAs, minimumScore) {
  parameters <- list(
    genesymbol = "",
    microrna = microRNAs,
    scoreClass = mapScore[minimumScore]
  )
  # ... send http POST
  res <- POST(paste(url, "/Http_U", sep = ""), body = parameters, encode = "form", verbose())
}

# make results-map as keyword - value
makeMap <- function(res) {
  
  ENTRY_DEL = "\001"
  KEY_DEL = "\002"
  response = httr::content(res, "text")
  arr = unlist(strsplit(response, ENTRY_DEL, fixed = TRUE))
  list_map <- list("")
  vec_map_names <- c("");
  
  for (str in arr) {
    arrKeyValue = unlist(strsplit(str, KEY_DEL, fixed = TRUE));
    
    if (length(arrKeyValue) > 1) {
      list_map[length(list_map) + 1] <- arrKeyValue[2]
      vec_map_names[length(vec_map_names) + 1] <- arrKeyValue[1]
    }
  }
  names(list_map) <- vec_map_names
  list_map
}

# Unidirectional search on MicroRNAs
if (TRUE) { 
  
  # Extract hsa from miRNAS_r_m file (miRDIP will only search hsa miRNAs)
  temp <- miRNAs_r_m$miRNA_ID[grep('hsa',miRNAs_r_m$miRNA_ID)]
  microRNAs = toString(temp)
  
  # Minimum Score - Use one of those:'Very High', 'High', 'Medium', 'Low' .
  minimumScore = "Medium"
  
  res <- unidirectionalSearchOnMicroRNAs(microRNAs, minimumScore)
  
  responseCode = status_code(res)
  if (responseCode != 200) {
    
    cat("Error: Response Code : ", responseCode, "\r\n")
  } else {
    
    list_map <- makeMap(res)
    
    # print results
    cat("\r\n", "Unidirectional Search on MicroRNA(s):", "\r\n")
    cat("Generated at: ", unlist(list_map["generated_at"]), "\r\n")
    cat("Micro RNAs: ", unlist(list_map["micro_rnas"]), "\r\n")
    cat("Minimum Score: ", unlist(list_map["minimum_score"]), "\r\n")
    cat("\r\n", "Results Size: ", unlist(list_map["results_size"]), "\r\n")
    # Suppress output here
    #cat("Results: \r\n", unlist(list_map["results"]), "\r\n") # formatted as tab - delimited spreadsheet
    
    #convert results into data table 
    temp <- list_map[["results"]] 
    mirdip <- read.table(text = temp, sep = "\t", colClasses = "character", header = TRUE)
    
    #Remove Pseudogene column 
    mirdip = mirdip[!names(mirdip) %in% c("Pseudogene")]
    rm(list_map, mapScore, res)
    
    #rename columns
    colnames(mirdip) = c("Gene_Symbol" , "Uniprot" , "MicroRNA" , "Integrated_Score" , "Number_of_Sources" ,	"Score_Class" ,	"Sources")
    
    #Combine two columns to create miRtarget column
    mirdip$miR_target <- paste(mirdip$MicroRNA, mirdip$Gene_Symbol) 
  }
}
head(mirdip)

######################################################################
## Identify validated targets (MTIs with experimental evidence)
# Load validated database (processed 2020-May-27 to include conversion to miRBase V22)
# Includes miRTarBase, TarBase V8 and miRecords
mtis_validated <- read.delim("C:/Users/niamhmconnolly/OneDrive - Royal College of Surgeons in Ireland/EpimiRNA/Network Map & Target ID/R/Database downloads/Validated dbs/f_validation_db_200527.txt",
                            stringsAsFactors=FALSE)

# Select miRNAs of interest from mtis_validated 
temp_valid <- sqldf("SELECT * FROM mtis_validated
                    WHERE EXISTS
                    (SELECT miRNAs_r_m.miRNA_ID FROM miRNAs_r_m WHERE 
                    mtis_validated.mirna = miRNAs_r_m.miRNA_ID)")

############################################
############################################
# Merge validated interactions with mirdip predicted interactions
merged_mirdip_valid <- merge(mirdip, temp_valid, by.x = 'miR_target',
                             by.y = 'miR_target', suffixes = c('.p','.v'), 
                             all.x = TRUE, all.y = TRUE)
rm(temp_valid,mtis_validated)

#drop columns not reqd
merged_mirdip_valid <- merged_mirdip_valid[ , !names(merged_mirdip_valid) %in% 
                                              c("mirna","gene", "miRNA", "target_gene")]

#split the miR_target column to get mirnas+genes for all rows
# Gives a warning if there is an entry with > 1 space (e.g. mmu-miR-124-3b CDK6 (HSA))
merged_mirdip_valid <- separate(data = merged_mirdip_valid, col = miR_target, 
                                into = c("miRNA","target_gene"), sep = " ")

#Recreate column with miRNA-gene target combined
merged_mirdip_valid$miR_target <- paste(merged_mirdip_valid$miRNA,
                                        merged_mirdip_valid$target_gene)
merged_mirdip_valid[, "miR_target"] <- merged_mirdip_valid$miR_target

head(merged_mirdip_valid)

# Number of miRNA-target interactions that are predicted and/or validated
nrow(merged_mirdip_valid)

######################################################################
##Calculate validated_score for validated MTIs (based on number of references) 

# Calculate validated score:
# If #pubs = 1, score = 0.5
# If #pubs = 2, score = 1.0
# If #pubs >=3, score = 1.5

validity_score = function(cnt){
  if(is.na(cnt)){score = 0}
  else if(as.numeric(cnt)==0){score = 0}
  else if(as.numeric(cnt)==1){score = 0.5}
  else if(as.numeric(cnt)==2){score = 1.0}
  else if(as.numeric(cnt)>=3){score = 1.5}
  return(score)
}

# x[16] here should be cnt column
merged_mirdip_valid$Validity_score = apply(merged_mirdip_valid,1, function(x){
  return(validity_score(x[16]))
})

######################################################################
# Sum predicted and validated scores to get Final_score

# Create numeric matrix from the two columns
tempMat <- cbind(as.numeric(merged_mirdip_valid$Integrated_Score), merged_mirdip_valid$Validity_score)
# Then sum these columns (na.rm = TRUE will ignore NA values)
merged_mirdip_valid$Final_score = rowSums(tempMat,na.rm = TRUE)

rm(mirdip, tempMat)

######################################################################
## Combine MTIs from different species 
# Code adapted from that originally written by Diana Smirnovova

# Create new dataframe for clean start ;)
mti_result_file <- merged_mirdip_valid

# Change all NAs to 0 (will use this later)
mti_result_file[is.na(mti_result_file)] <- 0

######################################################################
### Evidence classification, Identify MTIs with 'strong' evidence
# Labeling MTIs as STRONG or WEAK based on the Evidence and Experiment columns
# An MTI is labelled as STRONG if it contains either 'PCR', 'Western', or 'Luciferase Reporter Assay' 
# An MTI is labelled as PRED_ONLY if it has a predicted score but no experimental evidence
# Otherwise, the MTI is labelled weak (experimental evidence exists but not STRONG)

# NB. Spelling inaccuracies in Experiments column will not be accounted for

evidence_classification2 <- function(Integrated_Score, Experiments, cnt){
  if(Integrated_Score!=0 & (grepl('PCR', Experiments, ignore.case = TRUE) || 
                            grepl('Western', Experiments, ignore.case = TRUE) ||
                            grepl('Luciferase Reporter Assay', Experiments, ignore.case = TRUE))) {class = 'STRONG+PRED'}
  else if(Integrated_Score==0 & (grepl('PCR', Experiments, ignore.case = TRUE) || 
                                 grepl('Western', Experiments, ignore.case = TRUE) ||
                                 grepl('Luciferase Reporter Assay', Experiments, ignore.case = TRUE))) {class = 'STRONG_NOPRED'}
  else if (Integrated_Score!=0 & cnt==0) {class = 'PRED_ONLY'}
  else if (Integrated_Score!=0 & cnt!=0) {class= 'WEAK*+PRED'}
  else {class = 'WEAK*_NOPRED'}
  return(class)
}

# Call function
mti_result_file$expt_strength <- apply(mti_result_file, 1, function(x) {
  return(evidence_classification2(x[6], x[14], x[16]))
})

######################################################################
## Combine species so species-independent MTIs are all on one row

# Specify species to be analysed
species = c('hsa')

# Extract miR_target column
all_MTIs <- mti_result_file$miR_target

# Remove species info from all_MTIs
for (x in seq(1, length(species))){
  all_MTIs <- sub(paste(species[x],'-', sep=''), '', all_MTIs)
}

# Retain only unique MTIs
unique_MTIs <- unique(all_MTIs)

# Define column names to be extracted from mti_result_file
extract_cols = c('Integrated_Score', 'Score_Class', 'References', 'Experiments','cnt', 'Validity_score', 'Final_score',
'expt_strength')

# Define column names of output file
combined_cols = c('Iscore', 'conf_class', 'References', 'Experiments', 'pub_cnt', 'Vscore', 'Fscore', 'EvType')

# Combine columns function
combine_species_columns <- function(MTI, species, extract_cols, combined_cols){
  #Pre-define empty variable
  MTI_info <- c()
  for (x in seq(1, length(species))){
    # for each species obtain desired columns from one miRNA 
    selected_cols <- mti_result_file[(mti_result_file$miR_target==paste(species[x],'-', MTI,sep='')),extract_cols]
    # Name them as specified in 'prefix'
    colnames(selected_cols) <- paste(species[x], combined_cols, sep='_')
    # concatenate current species information (for one miRNA) into MTI_info
    MTI_info <- c(MTI_info, selected_cols)
  }
  #Transform MTI_info (list) for each species into character with all column headings combined
  MTI_info <- sapply(MTI_info, function(x) {ifelse(length(x) == 0, 0, x)})
  return(MTI_info)
}

# Combine all MTIS function
combine_all_MTIs <- function(MTI, species, extract_cols, combined_cols){
  # Pre-define empty variable
  mti_all <- c()
  for (x in seq(1, length(MTI))){
    # call combine_columns function for each miRNA
    mti_current <- rbind(unlist(combine_species_columns(MTI[x], species, extract_cols, combined_cols)))
    # Save column labels
    col_headings <- colnames(mti_current)
    # Concatenate MTI name with species info
    mti_current_named <- c(MTI[x],mti_current)
    # Append new MTI onto mti file
    mti_all <- rbind(mti_all,mti_current_named)
  }
  # Reapply column names
  colnames(mti_all) <- c("miR_target", col_headings)
  rownames(mti_all) <- NULL
  return(as.data.frame(mti_all))
}

# Call function to loop through all MTIs and all species
# Takes a while to run...
all_MTIs_combined <- combine_all_MTIs(unique_MTIs, species, extract_cols, combined_cols)

# Recreate miRNA and target gene columns, and add to all_* dataframe
miR_target = all_MTIs_combined$miR_target
miRNA <- apply(data.frame(miR_target), 1, function(x) unlist(strsplit(x, ' '))[1])
target_gene <- apply(data.frame(miR_target), 1, function(x) unlist(strsplit(x, ' '))[2]) 

all_MTIs_combined <- cbind(miRNA, target_gene, all_MTIs_combined)

# The number of unique (species-independent) MTIs
nrow(all_MTIs_combined)

######################################################################
## Add gene names to dataset
# Extract genenames using bitr
temp <- bitr(all_MTIs_combined$target_gene, 
             fromType="SYMBOL", toType="GENENAME", 
             OrgDb="org.Hs.eg.db")
# Remove duplicates
temp <- temp[!duplicated(temp$SYMBOL), ]
colnames(temp) <- c("SYMBOL", "gene_name")

# Merge files to add genename column
all_MTIs_combined <- merge(all_MTIs_combined, temp, by.x="target_gene", by.y="SYMBOL", all.x=TRUE)

rm(temp)

# Move gene name to beginning of dataframe
all_MTIs_combined <- dplyr::relocate(all_MTIs_combined, gene_name, .before=hsa_Iscore)

######################################################################
######################################################################
## Further analysis of MTIs - calculate total FScore
all_MTIs_combined$total_Fscore <- as.numeric(all_MTIs_combined$hsa_Fscore)
  
######################################################################
## Further analysis of MTIs - count occurrence of genes
# Edit so it outputs the number of genes rather than the interations

# Count occurrence of each gene
count_genes <- function(genename){
  # Pre-define empty variable
  gene_count <- c()
  # For each gene, count the number of time it occurs in the whole list
  for (x in seq(1,length(genename))){
    gene_count[x] <- sum(genename == genename[x], na.rm=TRUE)
  }
  return(as.numeric(gene_count))
}

# Call function
all_MTIs_combined$gene_count <- count_genes(all_MTIs_combined$target_gene)

######################################################################
## Further analysis of MTIs - count number of publications per MTI

# I don't know why you have to convert into character first and then numeric... ;-/
all_MTIs_combined$hsa_pub_cnt <- as.numeric(as.character(all_MTIs_combined$hsa_pub_cnt))
all_MTIs_combined$total_pub_cnt <- all_MTIs_combined$hsa_pub_cnt
  
# Convert scores to numeric
all_MTIs_combined$hsa_Iscore <- as.numeric(as.character(all_MTIs_combined$hsa_Iscore))
all_MTIs_combined$hsa_Vscore <- as.numeric(as.character(all_MTIs_combined$hsa_Vscore))

######################################################################
######################################################################
## Filter MTIs to retain only those of most interest

#Remove predicted only MTIs unless confidence is Very High
MTIs_filtered <- all_MTIs_combined[which(all_MTIs_combined$hsa_EvType != "PRED_ONLY" 
                                          | all_MTIs_combined$hsa_conf_class == "Very High"),]

# Remove Weak evidence with no prediction if only 1 pub
MTIs_filtered <- MTIs_filtered[which(MTIs_filtered$hsa_EvType != "WEAK*_NOPRED"
                                     | MTIs_filtered$hsa_pub_cnt > 1),]

######################################################################
# Redo gene count for filtered dataset
MTIs_filtered$gene_count_filtered <- count_genes(MTIs_filtered$target_gene)

# Write to file - SUpp Table 1
write.csv(MTIs_filtered, file="MTIs_filtered-run2b_v2.csv", row.names = FALSE)

######################################################################
######################################################################
## Create heatmap (numerical matrix) of miRNA-mRNA interactions 
# [Greenan et al, 2025; Figure 2E]
# First, extract relevant information
# Retain mRNA targets that are targeted by >X miRNA
heatmap_data <- MTIs_filtered[MTIs_filtered$gene_count_filtered>1, ]
#temp_mat <- heatmap_data[,c(1:2,27,28,38)]
temp_mat <- heatmap_data[,c("miRNA","target_gene","gene_count_filtered")]
#temp_mat$score <- heatmap_data$total_Fscore
temp_mat$score <- as.numeric(as.factor(heatmap_data$hsa_EvType))

# Create wide matrix
temp_mat_wide <- pivot_wider(temp_mat, names_from = miRNA, values_from = score)
#Replace 'NA's with 0 so they will be plotted
temp_mat_wide[is.na(temp_mat_wide)] <- 0
temp_mat_wide_mat <- as.data.frame(temp_mat_wide)
rownames(temp_mat_wide_mat) <- temp_mat_wide$target_gene
temp_mat_wide_mat <- temp_mat_wide_mat[,-1]

# Complex Heatmap
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
my_palette = colorRampPalette(c("white", "red"))(n=100)

Heatmap(name="EvType", as.matrix(temp_mat_wide_mat[,-c(1)]), col=my_palette, 
        row_names_gp = gpar(fontsize = 8)
        )

  
######################################################################
## Reactome Pathway enrichment analysis and visualisation
## If you use ReactomePA in published research, please cite:
## Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for 
## reactome pathway analysis and visualization. Molecular BioSystems 2016, 
## 12(2):477-479 
## See also: https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
## and https://bioconductor.org/packages/release/bioc/manuals/ReactomePA/man/ReactomePA.pdf
##
## clusterProfiler v3.14.3 
## For help: https://guangchuangyu.github.io/software/clusterProfiler
## If you use clusterProfiler in published research, please cite:
## Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: 
## an R package for comparing biological themes among gene clusters. OMICS: A
## Journal of Integrative Biology. 2012, 16(5):284-287.
  
  # Define RPEA function
  myRPEA <- function(MTIs_of_interest,pval_cutoff){
    # Create a vector containing only the target_gene column
    geneSymbol <- unique(MTIs_of_interest$target_gene) 
    
    # Convert geneSymbol into Entrez Gene ID for input to pathway enrichment
    # (to see Types that can be converted)
    keytypes(org.Hs.eg.db)
    entrez <- bitr(geneSymbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    
    # Number of genes as input to pathway enrichment analysis
    nrow(entrez)
    
    # Reactome pathway enrichment analysis
    # Only include pathways with between 10-500 members
    # *universe may need to be deleted on some systems*
    pathway <- enrichPathway(gene = entrez$ENTREZID, organism = "human", pvalueCutoff = pval_cutoff, 
                             pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10, 
                             maxGSSize = 500, readable = TRUE)
                            #universe)
    
    return(pathway)
  }
  
  # Function to generate dataframe from pathway output
  pathway_df <- function(pathway){
    
    rpea_output <- as.data.frame(pathway)
    
    # calculate ratio of genes targeted/#genes in pathway
    # number of genes in pathway
    rpea_output$pathwayGenes <- as.numeric(gsub("/.*$","",rpea_output$BgRatio))
    rpea_output$genePathwayRatio <- (rpea_output$Count/rpea_output$pathwayGenes)*100
    
    return(rpea_output)
  }
  
  # Define MTIs to be investigated!!!
  MTIs_of_interest <- MTIs_filtered
  #MTIs_of_interest <- MTIs_filtered[MTIs_filtered$gene_count_filtered>1, ]
  
  # Run pathway enrichment, convert output into dataframe
  pathwayRPEA <- myRPEA(MTIs_of_interest,pval_cutoff=0.05)
  pathwaydfRPEA <- pathway_df(pathwayRPEA)
  
  write.csv(pathwaydfRPEA,"RPEA-run2bV2_updated.csv", row.names = FALSE) # Supp Table 2
  
  # Visualise connections between pathways! (Enrichment map; formerly enrichMap)
  # if two categories have shared genes, they are connected; the edge width is scaled by the number of common genes between them.
  # pathway2 needs to be generated to avoid error in emapplot, then emapplot should be run on pathway2
  # Figure 2F, Greenan et al 2025
  dotplot(pathwayRPEA, x = "GeneRatio", showCategory = 30)
  
  #barplot(pathwayRPEA)
  #pathway2 <- pairwise_termsim(pathwayRPEA)
  #emapplot(pathway2, showCategory = nrow(pathway2))
  #emapplot(pathway2, showCategory = 30, cex_label_category=0.5)
  #treeplot(pathway2)
  #cnetplot(pathwayRPEA)
  