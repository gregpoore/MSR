#-----------------------------------------------------------------------------
# 03-kraken-db-filtering-and-ml-setup.R
# Copyright (c) 2023--, Greg Poore
# Purposes:
# - Filter Kraken db and related data based on Conterminator results
#-----------------------------------------------------------------------------

#----------------------------------------------------------#
# Load environments
#----------------------------------------------------------#

# Load dependencies
require(doMC)
require(plyr)
require(dplyr)
require(tidyr)
require(reshape2)
require(ggpubr)
require(ggsci)
require(ggrepel)
require(ggpmisc)

numCores <- detectCores()
registerDoMC(cores=numCores)

abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

#----------------------------------------------------------#
# Load data
#----------------------------------------------------------#
# Load raw Kraken data and metadata
load("Input_data/tcgaVbDataAndMetadataAndSNM_consolidated_Nov23.RData", verbose = TRUE)

# Load Kraken lineage file
krakenDbMetaJenya <- read.csv("Supporting_data/Kraken_db_filtering/repophlan_microbes_wscores_21Sept23.csv",
                            stringsAsFactors = FALSE, row.names = 1)
krakenDbMeta <- read.csv("Supporting_data/Kraken_db_filtering/repophlan_microbes_wscores.csv",
                         stringsAsFactors = FALSE, row.names = 1)
dim(krakenDbMeta) # 69286    35

# Extract missing genome IDs
krakenDbMetaJenya_New <- krakenDbMetaJenya[!(rownames(krakenDbMetaJenya) %in% rownames(krakenDbMeta)),]
dim(krakenDbMetaJenya_New) # 4245   35

# Find rows with taxa lineage --> only in "taxonomy" and "version_status" columns
# NOTE: Some taxonomy entries were shifted by one row and this code find them
sapply(krakenDbMetaJenya_New, function(x)
       length(which(grepl("\\|",x))))

# Create cleaned version of missing lineages
krakenDbMetaJenya_New_Form <- data.frame(taxonomy = c(
  krakenDbMetaJenya_New[which(grepl("\\|",krakenDbMetaJenya_New$taxonomy)),"taxonomy"],
  krakenDbMetaJenya_New[which(grepl("\\|",krakenDbMetaJenya_New$version_status)),"version_status"]
),
row.names = c(
  rownames( krakenDbMetaJenya_New[which(grepl("\\|",krakenDbMetaJenya_New$taxonomy)),] ),
  rownames( krakenDbMetaJenya_New[which(grepl("\\|",krakenDbMetaJenya_New$version_status)),] )
))

# Add cleaned missing lineages to full table
krakenDbMeta_Full <- rbind(krakenDbMeta[,"taxonomy",drop=FALSE],
                           krakenDbMetaJenya_New_Form)

krakenDbMeta_Full_Lineage <- krakenDbMeta_Full %>% #select(taxonomy) %>%
  mutate(taxForm = gsub("\\|t__.+","",taxonomy)) %>%
  mutate(taxForm = gsub("[k|p|c|o|f|g|s]__","",taxForm)) %>%
  separate(col = taxForm, sep = "\\|", 
           into = c("kingdom","phylum","class","order","family","genus","species")) %>%
  select(-taxonomy)

# Load WoLr1 OGUs with Conterminator hits
krakenDbConterminator <- read.csv("Supporting_data/Kraken_db_filtering/conterminator.problematic_taxids_resolved.result_conterm_prediction_unique_gIDs.csv",
                               stringsAsFactors = FALSE) %>% pull(gIDs)
length(krakenDbConterminator) # 583

# Subset Kraken db ranks by Conterminator hits
sum(krakenDbConterminator %in% rownames(krakenDbMeta_Full_Lineage)) # 583
sum(rownames(krakenDbMeta_Full_Lineage) %in% krakenDbConterminator) # 583

# Find which gIDs don't have matches <-- now fixed, so empty (Sept 21 2023)
missingGIDs <- krakenDbConterminator[which(! (krakenDbConterminator %in% rownames(krakenDbMeta_Full_Lineage)))]

krakenDbConterminatorRanks <- krakenDbMeta_Full_Lineage[intersect(krakenDbConterminator,
                                                             rownames(krakenDbMeta_Full_Lineage)),]
dim(krakenDbConterminatorRanks) # 583   7

# Subset to affected genera
krakenDbConterminatorRanksUniqueGenera <- sort(unique(krakenDbConterminatorRanks$genus[nzchar(krakenDbConterminatorRanks$genus)]))
length(krakenDbConterminatorRanksUniqueGenera) # 151

#----------------------------------------------------------#
# Find overlap with WIS genera
#----------------------------------------------------------#

# Find overlap with full WIS hit list
load("Interim_data/wisGeneraBactKnownUniqueRev_21Sept23.RData", verbose = TRUE)
# wisGeneraBactKnownUniqueRev
length(wisGeneraBactKnownUniqueRev) # 218
int_WIS_All_KrakenDbConterminator <- intersect(wisGeneraBactKnownUniqueRev, 
                                         krakenDbConterminatorRanksUniqueGenera)
length(int_WIS_All_KrakenDbConterminator) # 60 --> 27.5% (60 / 218)

# Find overlap within intersection of WIS and Kraken data
load("Interim_data/tcga-wis-overlapping-data-and-metadata-subset-24July22.RData", verbose=TRUE)
krakenFeatsWIS <- colnames(tcgaGenusKrakenAllFiltWIS)
length(krakenFeatsWIS) # 184

int_WIS_TCGA_KrakenDbConterminator <- intersect(krakenFeatsWIS, 
                                         krakenDbConterminatorRanksUniqueGenera)
length(int_WIS_TCGA_KrakenDbConterminator) # 60 --> 32.6% (60 / 184)

#----------------------------------------------------------#
# Identify and remove features in Kraken **genus** data
#----------------------------------------------------------#

# Subset to known genus level hits
vbDataBarnDFReconciledQC_Genus <- vbDataBarnDFReconciledQC
krakenFeats_Genus <- colnames(vbDataBarnDFReconciledQC_Genus)
krakenFeats_GenusForm <- gsub("^k__.+","",gsub("^k__.+\\.g__","",krakenFeats_Genus))

genusMatchesIdx <- unique(grep(paste(krakenDbConterminatorRanksUniqueGenera,collapse="|"), 
                               krakenFeats_GenusForm, value=FALSE))

vbDataBarnDFReconciledQC_Genus_NoHu <- vbDataBarnDFReconciledQC_Genus[,-c(genusMatchesIdx)]
dim(vbDataBarnDFReconciledQC_Genus) # 17625  1993
dim(vbDataBarnDFReconciledQC_Genus_NoHu) # 17625  1848 --> 145 dropped (7.3% of 1993)

# save(vbDataBarnDFReconciledQC_Genus,
#      vbDataBarnDFReconciledQC_Genus_NoHu,
#      metadataSamplesAllQC,
#      file = "Interim_data/kraken_data_genus_level_and_Conterminator_filtered_rerun_22Nov23.RData")

# Check if any zero-sum samples after subsetting
sum(rowSums(vbDataBarnDFReconciledQC_Genus_NoHu)==0) # =0

#----------------------------------------------------------#
# Subset to Illumina HiSeq
#----------------------------------------------------------#

metadataSamplesAllQC %>% count(data_submitting_center_label, experimental_strategy)
metadataSamplesAllQC %>% count(platform) # HiSeq accounts for 91.27% of samples
metadataSamplesAllQC %>% count(platform, data_submitting_center_label)

# Subset metadata to Illumina HiSeq samples and remove the only 2 Hopkins samples
metadataSamplesAllQC_HiSeq <- metadataSamplesAllQC %>%
  filter(platform == "Illumina HiSeq") %>% droplevels()

dim(metadataSamplesAllQC_HiSeq) # 16087    41

metadataSamplesAllQC_HiSeq %>% count(platform, data_submitting_center_label)

## Subset count data to Illumina HiSeq samples
# Genus
vbDataBarnDFReconciledQC_Genus_HiSeq <- vbDataBarnDFReconciledQC_Genus[rownames(metadataSamplesAllQC_HiSeq),]
vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq <- vbDataBarnDFReconciledQC_Genus_NoHu[rownames(metadataSamplesAllQC_HiSeq),]

dim(vbDataBarnDFReconciledQC_Genus_HiSeq) # 16087  1993
dim(vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq) # 16087  1848

#----------------------------------------------------------#
# Subset full data to individual seq centers
# metadataSamplesAllQC_HiSeq
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_HMS <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_BCM <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_MDA <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_WashU <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_Broad_WGS <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_UNC <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_CMS <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_Broad_RNA <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

#--------------------Subset **genus** raw count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
vbDataBarnDFReconciledQC_Genus_HiSeq_HMS <- vbDataBarnDFReconciledQC_Genus_HiSeq[rownames(metadataSamplesAllQC_HiSeq_HMS),]
vbDataBarnDFReconciledQC_Genus_HiSeq_BCM <- vbDataBarnDFReconciledQC_Genus_HiSeq[rownames(metadataSamplesAllQC_HiSeq_BCM),]
vbDataBarnDFReconciledQC_Genus_HiSeq_MDA <- vbDataBarnDFReconciledQC_Genus_HiSeq[rownames(metadataSamplesAllQC_HiSeq_MDA),]
vbDataBarnDFReconciledQC_Genus_HiSeq_WashU <- vbDataBarnDFReconciledQC_Genus_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WashU),]
vbDataBarnDFReconciledQC_Genus_HiSeq_Broad_WGS <- vbDataBarnDFReconciledQC_Genus_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
vbDataBarnDFReconciledQC_Genus_HiSeq_UNC <- vbDataBarnDFReconciledQC_Genus_HiSeq[rownames(metadataSamplesAllQC_HiSeq_UNC),]
vbDataBarnDFReconciledQC_Genus_HiSeq_CMS <- vbDataBarnDFReconciledQC_Genus_HiSeq[rownames(metadataSamplesAllQC_HiSeq_CMS),]
vbDataBarnDFReconciledQC_Genus_HiSeq_Broad_RNA <- vbDataBarnDFReconciledQC_Genus_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------Subset **genus** filtered count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_HMS <- vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_HMS),]
vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_BCM <- vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_BCM),]
vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_MDA <- vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_MDA),]
vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_WashU <- vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WashU),]
vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_Broad_WGS <- vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_UNC <- vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_UNC),]
vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_CMS <- vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_CMS),]
vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_Broad_RNA <- vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset genus raw count data
  vbDataBarnDFReconciledQC_Genus_HiSeq_HMS,
  vbDataBarnDFReconciledQC_Genus_HiSeq_BCM,
  vbDataBarnDFReconciledQC_Genus_HiSeq_MDA,
  vbDataBarnDFReconciledQC_Genus_HiSeq_WashU,
  vbDataBarnDFReconciledQC_Genus_HiSeq_Broad_WGS,
  vbDataBarnDFReconciledQC_Genus_HiSeq_UNC,
  vbDataBarnDFReconciledQC_Genus_HiSeq_CMS,
  vbDataBarnDFReconciledQC_Genus_HiSeq_Broad_RNA,
  
  # Subset genus filtered count data
  vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_HMS,
  vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_BCM,
  vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_MDA,
  vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_WashU,
  vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_Broad_WGS,
  vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_UNC,
  vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_CMS,
  vbDataBarnDFReconciledQC_Genus_NoHu_HiSeq_Broad_RNA,
  
  # Subset metadata
  metadataSamplesAllQC_HiSeq_HMS,
  metadataSamplesAllQC_HiSeq_BCM,
  metadataSamplesAllQC_HiSeq_MDA,
  metadataSamplesAllQC_HiSeq_WashU,
  metadataSamplesAllQC_HiSeq_Broad_WGS,
  metadataSamplesAllQC_HiSeq_UNC,
  metadataSamplesAllQC_HiSeq_CMS,
  metadataSamplesAllQC_HiSeq_Broad_RNA,
  
  file = "Interim_data/kraken_genus_raw_and_NoHu_data_for_ml_tcga_by_seq_center_rerun_22Nov23.RData")

# Scripts: SXXX

#----------------------------------------------------------#
# Plot results: Kraken Full vs Full_NoHu
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals
source("00-functions.R") # for barplotSummaryPerf() and barplotPerf() functions

## Load data
mlPerfAll10k_Allcancer_Kraken_NoHu <- read.csv("Supporting_scripts/S21-ML-10k-tcga-Kraken-raw-vs-NoHu-seqcenter-rerun/rep_perfML_10k_tcga_Kraken_raw_NoHu_rerun_ALL_22Nov23.csv", 
                                               stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

source("00-functions.R")
#----------------All----------------#
# Note: The results are *not* split by seqcenters, so only use summary plot FXN
barplotSummaryPerfNoHu(inputData = mlPerfAll10k_Allcancer_Kraken_NoHu,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Blood Derived Normal",
                   intFlag = FALSE,
                   plotWidthSingle = 8,
                   plotWidthCombined = 8,
                   fileNameString = "Kraken_huVSnohu_rerun")

barplotSummaryPerfNoHu(inputData = mlPerfAll10k_Allcancer_Kraken_NoHu,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor",
                   intFlag = FALSE,
                   plotWidthSingle = 12,
                   plotWidthCombined = 12,
                   outputData = TRUE,
                   fileNameString = "Kraken_huVSnohu_rerun")

barplotSummaryPerfNoHu(inputData = mlPerfAll10k_Allcancer_Kraken_NoHu,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                   intFlag = FALSE,
                   plotWidthSingle = 5,
                   plotWidthCombined = 5,
                   fileNameString = "Kraken_huVSnohu_rerun")

#----------------------------------------------------------#
# Identify and remove features in Kraken WIS **genus** data
#----------------------------------------------------------#

load("Interim_data/tcga-wis-overlapping-data-and-metadata-subset-24July22.RData", verbose=TRUE)
# tcgaGenusKrakenAllFiltWIS
# tcgaGenusKrakenQCFiltWIS
# metadataSamplesAll
# metadataSamplesAllQC

genusMatchesWISIdx <- unique(grep(paste(krakenDbConterminatorRanksUniqueGenera,collapse="|"), 
                                  colnames(tcgaGenusKrakenQCFiltWIS), value=FALSE))

tcgaGenusKrakenQCFiltWIS_NoHu <- tcgaGenusKrakenQCFiltWIS[,-c(genusMatchesWISIdx)]
dim(tcgaGenusKrakenQCFiltWIS) # 17625  184
dim(tcgaGenusKrakenQCFiltWIS_NoHu) # 17625  124 --> 60 dropped (32.6% of 184)

# save(tcgaGenusKrakenQCFiltWIS,
#      tcgaGenusKrakenQCFiltWIS_NoHu,
#      metadataSamplesAllQC,
#      file = "Interim_data/kraken_WIS_data_genus_level_and_Conterminator_filtered_rerun_22Nov23.RData")

# Check if any zero-sum samples after subsetting
sum(rowSums(tcgaGenusKrakenQCFiltWIS_NoHu)==0) # =0

## Subset count data to Illumina HiSeq samples (metadataSamplesAllQC_HiSeq object is from line ~141)
# Genus
tcgaGenusKrakenQCFiltWIS_HiSeq <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq),]
tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq <- tcgaGenusKrakenQCFiltWIS_NoHu[rownames(metadataSamplesAllQC_HiSeq),]

#--------------------Subset **genus** raw count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
tcgaGenusKrakenQCFiltWIS_HiSeq_HMS <- tcgaGenusKrakenQCFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_HMS),]
tcgaGenusKrakenQCFiltWIS_HiSeq_BCM <- tcgaGenusKrakenQCFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_BCM),]
tcgaGenusKrakenQCFiltWIS_HiSeq_MDA <- tcgaGenusKrakenQCFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_MDA),]
tcgaGenusKrakenQCFiltWIS_HiSeq_WashU <- tcgaGenusKrakenQCFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WashU),]
tcgaGenusKrakenQCFiltWIS_HiSeq_Broad_WGS <- tcgaGenusKrakenQCFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
tcgaGenusKrakenQCFiltWIS_HiSeq_UNC <- tcgaGenusKrakenQCFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_UNC),]
tcgaGenusKrakenQCFiltWIS_HiSeq_CMS <- tcgaGenusKrakenQCFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_CMS),]
tcgaGenusKrakenQCFiltWIS_HiSeq_Broad_RNA <- tcgaGenusKrakenQCFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------Subset **genus** filtered count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_HMS <- tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_HMS),]
tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_BCM <- tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_BCM),]
tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_MDA <- tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_MDA),]
tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_WashU <- tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WashU),]
tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_Broad_WGS <- tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_UNC <- tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_UNC),]
tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_CMS <- tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_CMS),]
tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_Broad_RNA <- tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset genus raw count data
  tcgaGenusKrakenQCFiltWIS_HiSeq_HMS,
  tcgaGenusKrakenQCFiltWIS_HiSeq_BCM,
  tcgaGenusKrakenQCFiltWIS_HiSeq_MDA,
  tcgaGenusKrakenQCFiltWIS_HiSeq_WashU,
  tcgaGenusKrakenQCFiltWIS_HiSeq_Broad_WGS,
  tcgaGenusKrakenQCFiltWIS_HiSeq_UNC,
  tcgaGenusKrakenQCFiltWIS_HiSeq_CMS,
  tcgaGenusKrakenQCFiltWIS_HiSeq_Broad_RNA,
  
  # Subset genus filtered count data
  tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_HMS,
  tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_BCM,
  tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_MDA,
  tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_WashU,
  tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_Broad_WGS,
  tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_UNC,
  tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_CMS,
  tcgaGenusKrakenQCFiltWIS_NoHu_HiSeq_Broad_RNA,
  
  # Subset metadata
  metadataSamplesAllQC_HiSeq_HMS,
  metadataSamplesAllQC_HiSeq_BCM,
  metadataSamplesAllQC_HiSeq_MDA,
  metadataSamplesAllQC_HiSeq_WashU,
  metadataSamplesAllQC_HiSeq_Broad_WGS,
  metadataSamplesAllQC_HiSeq_UNC,
  metadataSamplesAllQC_HiSeq_CMS,
  metadataSamplesAllQC_HiSeq_Broad_RNA,
  
  file = "Interim_data/kraken_WIS_genus_raw_and_NoHu_data_for_ml_tcga_by_seq_center_rerun_22Nov23.RData")

#----------------------------------------------------------#
# Plot results: Kraken WIS vs WIS_NoHu
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals
source("00-functions.R") # for barplotSummaryPerf() and barplotPerf() functions

## Load data
mlPerfAll10k_Allcancer_Kraken_WIS_NoHu <- read.csv("Supporting_scripts/S22-ML-10k-tcga-Kraken-WIS-raw-vs-NoHu-seqcenter-rerun/rep_perfML_10k_tcga_Kraken_WIS_raw_NoHu_rerun_ALL_22Nov23.csv", 
                                               stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

source("00-functions.R")
#----------------All----------------#
# Note: The results are *not* split by seqcenters, so only use summary plot FXN
barplotSummaryPerfNoHu(inputData = mlPerfAll10k_Allcancer_Kraken_WIS_NoHu,
                       seqCenterAbbrev="All",
                       sampleTypeInput = "Blood Derived Normal",
                       prefix2Remove = "tcgaGenusKrakenQCFiltWIS_",
                       intFlag = FALSE,
                       plotWidthSingle = 8,
                       plotWidthCombined = 8,
                       fileNameString = "Kraken_WIS_huVSnohu_rerun")

barplotSummaryPerfNoHu(inputData = mlPerfAll10k_Allcancer_Kraken_WIS_NoHu,
                       seqCenterAbbrev="All",
                       sampleTypeInput = "Primary Tumor",
                       prefix2Remove = "tcgaGenusKrakenQCFiltWIS_",
                       intFlag = FALSE,
                       outputData = TRUE,
                       plotWidthSingle = 12,
                       plotWidthCombined = 12,
                       fileNameString = "Kraken_WIS_huVSnohu_rerun")

barplotSummaryPerfNoHu(inputData = mlPerfAll10k_Allcancer_Kraken_WIS_NoHu,
                       seqCenterAbbrev="All",
                       sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                       prefix2Remove = "tcgaGenusKrakenQCFiltWIS_",
                       intFlag = FALSE,
                       outputData = TRUE,
                       plotWidthSingle = 5,
                       plotWidthCombined = 5,
                       fileNameString = "Kraken_WIS_huVSnohu_rerun")
