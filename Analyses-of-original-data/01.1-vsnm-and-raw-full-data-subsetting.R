# 01.1-vsnm-and-raw-full-data-subsetting.R
# Author: Greg Poore
# Date: Aug 21, 2023
# Purpose: Intersect Kraken TCGA data with WIS taxa

#-------------------------------#
# Load dependencies
require(doMC)
require(plyr)
require(dplyr)
require(tibble)
require(reshape2)
require(phyloseq)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import TCGA data
#----------------------------------------------------------#

load("Input_data/tcgaVbDataAndMetadataAndSNM.RData", verbose = TRUE)
# # metadataSamplesAll
# # metadataSamplesAllQC
# # metadataSamplesAllQCSurvival
# # vbDataBarnDFReconciled
# # vbDataBarnDFReconciledQC
# # vbDataBarnDFReconciledQCSurvival
# # snmDataSampleType
# 
load("Input_data/snmDataSampleTypeWithExpStrategyFINAL.RData", verbose = TRUE)
# # snmDataSampleTypeWithExpStrategy
# 
save(metadataSamplesAll,
     metadataSamplesAllQC,
     vbDataBarnDFReconciled,
     vbDataBarnDFReconciledQC,
     snmDataSampleTypeWithExpStrategy,
     file = "Input_data/tcgaVbDataAndMetadataAndSNM_consolidated_Nov23.RData")

load("Input_data/tcgaVbDataAndMetadataAndSNM_consolidated_Nov23.RData")

#--------------------------------------------------------------------------------------------------------------------#
# Separate raw data into seq center-experimental strategy groups (to preclude needing batch correction)
#--------------------------------------------------------------------------------------------------------------------#
metadataSamplesAllQC %>% count(data_submitting_center_label, experimental_strategy)
metadataSamplesAllQC %>% count(platform) # HiSeq accounts for 91.27% of samples
metadataSamplesAllQC %>% count(platform, data_submitting_center_label)

# Subset metadata to Illumina HiSeq samples
metadataSamplesAllQC_HiSeq <- metadataSamplesAllQC %>%
  filter(platform == "Illumina HiSeq") %>% droplevels()

# Subset count data to Illumina HiSeq samples
vbDataBarnDFReconciledQC_HiSeq <- vbDataBarnDFReconciledQC[rownames(metadataSamplesAllQC_HiSeq),]
snmDataSampleTypeWithExpStrategy_HiSeq <- snmDataSampleTypeWithExpStrategy[rownames(metadataSamplesAllQC_HiSeq),]

dim(vbDataBarnDFReconciledQC_HiSeq) # 16087  1993
dim(snmDataSampleTypeWithExpStrategy_HiSeq) # 16087  1795

# Optionally intersect features
feats2Keep <- intersect(colnames(vbDataBarnDFReconciledQC_HiSeq),
                        colnames(snmDataSampleTypeWithExpStrategy_HiSeq))

vbDataBarnDFReconciledQC_Int_HiSeq <- vbDataBarnDFReconciledQC_HiSeq[,colnames(vbDataBarnDFReconciledQC_HiSeq) %in% feats2Keep]
dim(vbDataBarnDFReconciledQC_Int_HiSeq) # 16087  1795

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_HMS <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
metadataSamplesAllQC_HiSeq_BCM <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Baylor College of Medicine") %>% droplevels()
metadataSamplesAllQC_HiSeq_MDA <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% droplevels()
metadataSamplesAllQC_HiSeq_WashU <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Washington University School of Medicine") %>% droplevels()
metadataSamplesAllQC_HiSeq_Broad_WGS <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_UNC <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "University of North Carolina") %>% droplevels()
metadataSamplesAllQC_HiSeq_CMS <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% droplevels()
metadataSamplesAllQC_HiSeq_Broad_RNA <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

##--------------------Subset count data by seq center--------------------##
#--------------------vbDataBarnDFReconciledQC_HiSeq--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
vbDataBarnDFReconciledQC_HiSeq_HMS <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_HMS),]
vbDataBarnDFReconciledQC_HiSeq_BCM <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_BCM),]
vbDataBarnDFReconciledQC_HiSeq_MDA <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_MDA),]
vbDataBarnDFReconciledQC_HiSeq_WashU <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WashU),]
vbDataBarnDFReconciledQC_HiSeq_Broad_WGS <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
vbDataBarnDFReconciledQC_HiSeq_UNC <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_UNC),]
vbDataBarnDFReconciledQC_HiSeq_CMS <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_CMS),]
vbDataBarnDFReconciledQC_HiSeq_Broad_RNA <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------vbDataBarnDFReconciledQC_Int_HiSeq--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
vbDataBarnDFReconciledQC_Int_HiSeq_HMS <- vbDataBarnDFReconciledQC_Int_HiSeq[rownames(metadataSamplesAllQC_HiSeq_HMS),]
vbDataBarnDFReconciledQC_Int_HiSeq_BCM <- vbDataBarnDFReconciledQC_Int_HiSeq[rownames(metadataSamplesAllQC_HiSeq_BCM),]
vbDataBarnDFReconciledQC_Int_HiSeq_MDA <- vbDataBarnDFReconciledQC_Int_HiSeq[rownames(metadataSamplesAllQC_HiSeq_MDA),]
vbDataBarnDFReconciledQC_Int_HiSeq_WashU <- vbDataBarnDFReconciledQC_Int_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WashU),]
vbDataBarnDFReconciledQC_Int_HiSeq_Broad_WGS <- vbDataBarnDFReconciledQC_Int_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
vbDataBarnDFReconciledQC_Int_HiSeq_UNC <- vbDataBarnDFReconciledQC_Int_HiSeq[rownames(metadataSamplesAllQC_HiSeq_UNC),]
vbDataBarnDFReconciledQC_Int_HiSeq_CMS <- vbDataBarnDFReconciledQC_Int_HiSeq[rownames(metadataSamplesAllQC_HiSeq_CMS),]
vbDataBarnDFReconciledQC_Int_HiSeq_Broad_RNA <- vbDataBarnDFReconciledQC_Int_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------snmDataSampleTypeWithExpStrategy_HiSeq--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
snmDataSampleTypeWithExpStrategy_HiSeq_HMS <- snmDataSampleTypeWithExpStrategy_HiSeq[rownames(metadataSamplesAllQC_HiSeq_HMS),]
snmDataSampleTypeWithExpStrategy_HiSeq_BCM <- snmDataSampleTypeWithExpStrategy_HiSeq[rownames(metadataSamplesAllQC_HiSeq_BCM),]
snmDataSampleTypeWithExpStrategy_HiSeq_MDA <- snmDataSampleTypeWithExpStrategy_HiSeq[rownames(metadataSamplesAllQC_HiSeq_MDA),]
snmDataSampleTypeWithExpStrategy_HiSeq_WashU <- snmDataSampleTypeWithExpStrategy_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WashU),]
snmDataSampleTypeWithExpStrategy_HiSeq_Broad_WGS <- snmDataSampleTypeWithExpStrategy_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
snmDataSampleTypeWithExpStrategy_HiSeq_UNC <- snmDataSampleTypeWithExpStrategy_HiSeq[rownames(metadataSamplesAllQC_HiSeq_UNC),]
snmDataSampleTypeWithExpStrategy_HiSeq_CMS <- snmDataSampleTypeWithExpStrategy_HiSeq[rownames(metadataSamplesAllQC_HiSeq_CMS),]
snmDataSampleTypeWithExpStrategy_HiSeq_Broad_RNA <- snmDataSampleTypeWithExpStrategy_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------Subset metadata and count data by WGS (for multi-class classification)--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_WGS <- metadataSamplesAllQC_HiSeq %>% filter(experimental_strategy == "WGS") %>% droplevels()

vbDataBarnDFReconciledQC_HiSeq_WGS <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WGS),]
vbDataBarnDFReconciledQC_Int_HiSeq_WGS <- vbDataBarnDFReconciledQC_Int_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WGS),]
snmDataSampleTypeWithExpStrategy_HiSeq_WGS <- snmDataSampleTypeWithExpStrategy_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WGS),]

# Save data for ML
save(metadataSamplesAllQC_HiSeq_WGS,
     vbDataBarnDFReconciledQC_HiSeq_WGS,
     vbDataBarnDFReconciledQC_Int_HiSeq_WGS,
     snmDataSampleTypeWithExpStrategy_HiSeq_WGS,
  file = "Interim_data/data_for_multiclass_ml_tcga_wgs_raw_vsnm_21Aug23.RData")

# Scripts: SXXX

#----------------------------------------------------#
# Save data for ML
#----------------------------------------------------#

save(# Subset raw count data
  vbDataBarnDFReconciledQC_HiSeq_HMS,
  vbDataBarnDFReconciledQC_HiSeq_BCM,
  vbDataBarnDFReconciledQC_HiSeq_MDA,
  vbDataBarnDFReconciledQC_HiSeq_WashU,
  vbDataBarnDFReconciledQC_HiSeq_Broad_WGS,
  vbDataBarnDFReconciledQC_HiSeq_UNC,
  vbDataBarnDFReconciledQC_HiSeq_CMS,
  vbDataBarnDFReconciledQC_HiSeq_Broad_RNA,
  
  # Subset raw count data with same features as VSNM
  vbDataBarnDFReconciledQC_Int_HiSeq_HMS,
  vbDataBarnDFReconciledQC_Int_HiSeq_BCM,
  vbDataBarnDFReconciledQC_Int_HiSeq_MDA,
  vbDataBarnDFReconciledQC_Int_HiSeq_WashU,
  vbDataBarnDFReconciledQC_Int_HiSeq_Broad_WGS,
  vbDataBarnDFReconciledQC_Int_HiSeq_UNC,
  vbDataBarnDFReconciledQC_Int_HiSeq_CMS,
  vbDataBarnDFReconciledQC_Int_HiSeq_Broad_RNA,
  
  # VSNM batch corrected data
  snmDataSampleTypeWithExpStrategy_HiSeq_HMS,
  snmDataSampleTypeWithExpStrategy_HiSeq_BCM,
  snmDataSampleTypeWithExpStrategy_HiSeq_MDA,
  snmDataSampleTypeWithExpStrategy_HiSeq_WashU,
  snmDataSampleTypeWithExpStrategy_HiSeq_Broad_WGS,
  snmDataSampleTypeWithExpStrategy_HiSeq_UNC,
  snmDataSampleTypeWithExpStrategy_HiSeq_CMS,
  snmDataSampleTypeWithExpStrategy_HiSeq_Broad_RNA,
  
  # Subset metadata
  metadataSamplesAllQC_HiSeq_HMS,
  metadataSamplesAllQC_HiSeq_BCM,
  metadataSamplesAllQC_HiSeq_MDA,
  metadataSamplesAllQC_HiSeq_WashU,
  metadataSamplesAllQC_HiSeq_Broad_WGS,
  metadataSamplesAllQC_HiSeq_UNC,
  metadataSamplesAllQC_HiSeq_CMS,
  metadataSamplesAllQC_HiSeq_Broad_RNA,
  file = "Interim_data/raw_and_vsnm_data_for_ml_tcga_by_seq_center_and_experimental_strategy_21Aug23.RData")

# Scripts: SXXX