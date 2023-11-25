#-----------------------------------------------------------------------------
# 04-control-validation-analyses.R
# Copyright (c) 2022--, Greg Poore
# Purposes:
# - Split data into 2 stratified halves, normalize using VSNM, train ML models, and test on each other
# - Scramble metadata and shuffle samples and re-run ML as controls
# - Compare performance of ML controls to actual data for raw data (per seq center) and VSNM (all TCGA)
#-----------------------------------------------------------------------------

#----------------------------------------------------------#
# Load environments
#----------------------------------------------------------#

# Load dependencies
require(devtools)
require(doMC)
require(phyloseq)
require(microbiome)
require(vegan)
require(plyr)
require(dplyr)
require(reshape2)
require(ggpubr)
require(ggsci)
require(ANCOMBC)
require(biomformat)
require(Rhdf5lib)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Load TCGA WIS feature subset data
#----------------------------------------------------------#

load("Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy25July22.RData", verbose = T)

#----------------------------------------------------------#
# Scramble labels and data: raw data - *run once*
#----------------------------------------------------------#

#-----------------------Scramble metadata-----------------------#
# All (VSNM)
metadataSamplesAllQC_scrambled <- metadataSamplesAllQC
set.seed(42)
metadataSamplesAllQC_scrambled$disease_type <- sample(metadataSamplesAllQC_scrambled$disease_type)
set.seed(42)
metadataSamplesAllQC_scrambled$sample_type <- sample(metadataSamplesAllQC_scrambled$sample_type)
# HMS
metadataSamplesAllQC_HiSeq_HMS_scrambled <- metadataSamplesAllQC_HiSeq_HMS
set.seed(42)
metadataSamplesAllQC_HiSeq_HMS_scrambled$disease_type <- sample(metadataSamplesAllQC_HiSeq_HMS_scrambled$disease_type)
set.seed(42)
metadataSamplesAllQC_HiSeq_HMS_scrambled$sample_type <- sample(metadataSamplesAllQC_HiSeq_HMS_scrambled$sample_type)
# BCM
metadataSamplesAllQC_HiSeq_BCM_scrambled <- metadataSamplesAllQC_HiSeq_BCM
set.seed(42)
metadataSamplesAllQC_HiSeq_BCM_scrambled$disease_type <- sample(metadataSamplesAllQC_HiSeq_BCM_scrambled$disease_type)
set.seed(42)
metadataSamplesAllQC_HiSeq_BCM_scrambled$sample_type <- sample(metadataSamplesAllQC_HiSeq_BCM_scrambled$sample_type)
# MDA
metadataSamplesAllQC_HiSeq_MDA_scrambled <- metadataSamplesAllQC_HiSeq_MDA
set.seed(42)
metadataSamplesAllQC_HiSeq_MDA_scrambled$disease_type <- sample(metadataSamplesAllQC_HiSeq_MDA_scrambled$disease_type)
set.seed(42)
metadataSamplesAllQC_HiSeq_MDA_scrambled$sample_type <- sample(metadataSamplesAllQC_HiSeq_MDA_scrambled$sample_type)
# WashU
metadataSamplesAllQC_HiSeq_WashU_scrambled <- metadataSamplesAllQC_HiSeq_WashU
set.seed(42)
metadataSamplesAllQC_HiSeq_WashU_scrambled$disease_type <- sample(metadataSamplesAllQC_HiSeq_WashU_scrambled$disease_type)
set.seed(42)
metadataSamplesAllQC_HiSeq_WashU_scrambled$sample_type <- sample(metadataSamplesAllQC_HiSeq_WashU_scrambled$sample_type)
# Broad_WGS
metadataSamplesAllQC_HiSeq_Broad_WGS_scrambled <- metadataSamplesAllQC_HiSeq_Broad_WGS
set.seed(42)
metadataSamplesAllQC_HiSeq_Broad_WGS_scrambled$disease_type <- sample(metadataSamplesAllQC_HiSeq_Broad_WGS_scrambled$disease_type)
set.seed(42)
metadataSamplesAllQC_HiSeq_Broad_WGS_scrambled$sample_type <- sample(metadataSamplesAllQC_HiSeq_Broad_WGS_scrambled$sample_type)
# UNC
metadataSamplesAllQC_HiSeq_UNC_scrambled <- metadataSamplesAllQC_HiSeq_UNC
set.seed(42)
metadataSamplesAllQC_HiSeq_UNC_scrambled$disease_type <- sample(metadataSamplesAllQC_HiSeq_UNC_scrambled$disease_type)
set.seed(42)
metadataSamplesAllQC_HiSeq_UNC_scrambled$sample_type <- sample(metadataSamplesAllQC_HiSeq_UNC_scrambled$sample_type)
# CMS
metadataSamplesAllQC_HiSeq_CMS_scrambled <- metadataSamplesAllQC_HiSeq_CMS
set.seed(42)
metadataSamplesAllQC_HiSeq_CMS_scrambled$disease_type <- sample(metadataSamplesAllQC_HiSeq_CMS_scrambled$disease_type)
set.seed(42)
metadataSamplesAllQC_HiSeq_CMS_scrambled$sample_type <- sample(metadataSamplesAllQC_HiSeq_CMS_scrambled$sample_type)
# Broad_RNA
metadataSamplesAllQC_HiSeq_Broad_RNA_scrambled <- metadataSamplesAllQC_HiSeq_Broad_RNA
set.seed(42)
metadataSamplesAllQC_HiSeq_Broad_RNA_scrambled$disease_type <- sample(metadataSamplesAllQC_HiSeq_Broad_RNA_scrambled$disease_type)
set.seed(42)
metadataSamplesAllQC_HiSeq_Broad_RNA_scrambled$sample_type <- sample(metadataSamplesAllQC_HiSeq_Broad_RNA_scrambled$sample_type)

#-----------------------Shuffle raw data-----------------------#
# All
set.seed(42)
vsnmDataGenusKrakenQCFiltWIS_shuffled <- vsnmDataGenusKrakenQCFiltWIS[sample(nrow(vsnmDataGenusKrakenQCFiltWIS)),]
rownames(vsnmDataGenusKrakenQCFiltWIS_shuffled) <- rownames(vsnmDataGenusKrakenQCFiltWIS)

# HMS
set.seed(42)
tcgaGenusKrakenQCFiltWIS_HMS_shuffled <- tcgaGenusKrakenQCFiltWIS_HMS[sample(nrow(tcgaGenusKrakenQCFiltWIS_HMS)),]
rownames(tcgaGenusKrakenQCFiltWIS_HMS_shuffled) <- rownames(tcgaGenusKrakenQCFiltWIS_HMS)
# BCM
set.seed(42)
tcgaGenusKrakenQCFiltWIS_BCM_shuffled <- tcgaGenusKrakenQCFiltWIS_BCM[sample(nrow(tcgaGenusKrakenQCFiltWIS_BCM)),]
rownames(tcgaGenusKrakenQCFiltWIS_BCM_shuffled) <- rownames(tcgaGenusKrakenQCFiltWIS_BCM)
# MDA
set.seed(42)
tcgaGenusKrakenQCFiltWIS_MDA_shuffled <- tcgaGenusKrakenQCFiltWIS_MDA[sample(nrow(tcgaGenusKrakenQCFiltWIS_MDA)),]
rownames(tcgaGenusKrakenQCFiltWIS_MDA_shuffled) <- rownames(tcgaGenusKrakenQCFiltWIS_MDA)
# WashU
set.seed(42)
tcgaGenusKrakenQCFiltWIS_WashU_shuffled <- tcgaGenusKrakenQCFiltWIS_WashU[sample(nrow(tcgaGenusKrakenQCFiltWIS_WashU)),]
rownames(tcgaGenusKrakenQCFiltWIS_WashU_shuffled) <- rownames(tcgaGenusKrakenQCFiltWIS_WashU)
# Broad_WGS
set.seed(42)
tcgaGenusKrakenQCFiltWIS_Broad_WGS_shuffled <- tcgaGenusKrakenQCFiltWIS_Broad_WGS[sample(nrow(tcgaGenusKrakenQCFiltWIS_Broad_WGS)),]
rownames(tcgaGenusKrakenQCFiltWIS_Broad_WGS_shuffled) <- rownames(tcgaGenusKrakenQCFiltWIS_Broad_WGS)
# UNC
set.seed(42)
tcgaGenusKrakenQCFiltWIS_UNC_shuffled <- tcgaGenusKrakenQCFiltWIS_UNC[sample(nrow(tcgaGenusKrakenQCFiltWIS_UNC)),]
rownames(tcgaGenusKrakenQCFiltWIS_UNC_shuffled) <- rownames(tcgaGenusKrakenQCFiltWIS_UNC)
# CMS
set.seed(42)
tcgaGenusKrakenQCFiltWIS_CMS_shuffled <- tcgaGenusKrakenQCFiltWIS_CMS[sample(nrow(tcgaGenusKrakenQCFiltWIS_CMS)),]
rownames(tcgaGenusKrakenQCFiltWIS_CMS_shuffled) <- rownames(tcgaGenusKrakenQCFiltWIS_CMS)
# Broad_RNA
set.seed(42)
tcgaGenusKrakenQCFiltWIS_Broad_RNA_shuffled <- tcgaGenusKrakenQCFiltWIS_Broad_RNA[sample(nrow(tcgaGenusKrakenQCFiltWIS_Broad_RNA)),]
rownames(tcgaGenusKrakenQCFiltWIS_Broad_RNA_shuffled) <- rownames(tcgaGenusKrakenQCFiltWIS_Broad_RNA)

save(metadataSamplesAllQC_scrambled,
     metadataSamplesAllQC_HiSeq_HMS_scrambled,
     metadataSamplesAllQC_HiSeq_BCM_scrambled,
     metadataSamplesAllQC_HiSeq_MDA_scrambled,
     metadataSamplesAllQC_HiSeq_WashU_scrambled,
     metadataSamplesAllQC_HiSeq_Broad_WGS_scrambled,
     metadataSamplesAllQC_HiSeq_UNC_scrambled,
     metadataSamplesAllQC_HiSeq_CMS_scrambled,
     metadataSamplesAllQC_HiSeq_Broad_RNA_scrambled,
     # Scrambled meta normal data
     vsnmDataGenusKrakenQCFiltWIS,
     tcgaGenusKrakenQCFiltWIS_HMS,
     tcgaGenusKrakenQCFiltWIS_BCM,
     tcgaGenusKrakenQCFiltWIS_MDA,
     tcgaGenusKrakenQCFiltWIS_WashU,
     tcgaGenusKrakenQCFiltWIS_Broad_WGS,
     tcgaGenusKrakenQCFiltWIS_UNC,
     tcgaGenusKrakenQCFiltWIS_CMS,
     tcgaGenusKrakenQCFiltWIS_Broad_RNA,
     # Shuffled data normal metadata
     metadataSamplesAllQC,
     metadataSamplesAllQC_HiSeq_HMS,
     metadataSamplesAllQC_HiSeq_BCM,
     metadataSamplesAllQC_HiSeq_MDA,
     metadataSamplesAllQC_HiSeq_WashU,
     metadataSamplesAllQC_HiSeq_Broad_WGS,
     metadataSamplesAllQC_HiSeq_UNC,
     metadataSamplesAllQC_HiSeq_CMS,
     metadataSamplesAllQC_HiSeq_Broad_RNA,
     
     vsnmDataGenusKrakenQCFiltWIS_shuffled,
     tcgaGenusKrakenQCFiltWIS_HMS_shuffled,
     tcgaGenusKrakenQCFiltWIS_BCM_shuffled,
     tcgaGenusKrakenQCFiltWIS_MDA_shuffled,
     tcgaGenusKrakenQCFiltWIS_WashU_shuffled,
     tcgaGenusKrakenQCFiltWIS_Broad_WGS_shuffled,
     tcgaGenusKrakenQCFiltWIS_UNC_shuffled,
     tcgaGenusKrakenQCFiltWIS_CMS_shuffled,
     tcgaGenusKrakenQCFiltWIS_Broad_RNA_shuffled,
     file = "Interim_data/shuffled_controls_raw_data_TCGA_25July22.RData")

#----------------------------------------------------------#
# Run ML
#----------------------------------------------------------#

# Script: S03

#----------------------------------------------------------#
# Plot scrambled ML perf
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM <- read.csv("Interim_data/rep_perfML_10k_tcga_wis_feature_subset_scrambled_and_shuffled_controls_ALL_25July22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM <- mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM[,!(colnames(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM) == "X")]
colnames(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$minorityClassName == "SolidTissueNormal",
                                                             yes=mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$majorityClassSize/(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$minorityClassSize+mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$majorityClassSize),
                                                             no=mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$minorityClassSize/(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$minorityClassSize+mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$majorityClassSize))
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
table(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName)
table(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName)
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_HMS_scrambled"] <- "HMS metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_HMS_shuffled"] <- "HMS count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_MDA_scrambled"] <- "MDA metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_MDA_shuffled"] <- "MDA count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_BCM_scrambled"] <- "BCM metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_BCM_shuffled"] <- "BCM count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_WashU_scrambled"] <- "WashU metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_WashU_shuffled"] <- "WashU count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_Broad_WGS_scrambled"] <- "Broad metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_Broad_WGS_shuffled"] <- "Broad count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_UNC_scrambled"] <- "UNC metadata labels scrambled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_UNC_shuffled"] <- "UNC count data shuffled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_CMS_scrambled"] <- "CMS metadata labels scrambled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_HiSeq_CMS_shuffled"] <- "CMS count data shuffled (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_scrambled"] <- "VSNM all metadata labels scrambled (WGS+RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metadataSamplesAllQC_shuffled"] <- "VSNM all count data shuffled (WGS+RNA-Seq)"

mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName <- factor(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName,
                                                                levels = c("VSNM all metadata labels scrambled (WGS+RNA-Seq)",
                                                                           "VSNM all count data shuffled (WGS+RNA-Seq)",
                                                                           "HMS metadata labels scrambled (WGS)",
                                                                           "HMS count data shuffled (WGS)",
                                                                           "MDA metadata labels scrambled (WGS)",
                                                                           "MDA count data shuffled (WGS)",
                                                                           "BCM metadata labels scrambled (WGS)",
                                                                           "BCM count data shuffled (WGS)",
                                                                           "WashU metadata labels scrambled (WGS)",
                                                                           "WashU count data shuffled (WGS)",
                                                                           "Broad metadata labels scrambled (WGS)",
                                                                           "Broad count data shuffled (WGS)",
                                                                           "UNC metadata labels scrambled (RNA-Seq)",
                                                                           "UNC count data shuffled (RNA-Seq)",
                                                                           "CMS metadata labels scrambled (RNA-Seq)",
                                                                           "CMS count data shuffled (RNA-Seq)"))
table(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName)

#----------------------------------------------------------#
# Overlay actual and (scrambled) control performance per seq center
#----------------------------------------------------------#

load("Interim_data/mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS_26July22.RData", verbose = TRUE)

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt <- rbind(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM,
                                                    mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt$datasetName <- factor(mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt$datasetName,
                                                                      levels = rev(levels(mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt$datasetName)))

#----------------------------------------------------------#
# Boxplot summaries of overlay per seq center
#----------------------------------------------------------#

mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short <- mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(grepl("genus|shuffled|scrambled",datasetName)) %>% droplevels()
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub(" \\(WGS\\)| \\(RNA-Seq\\)| \\(WGS\\+RNA-Seq\\)","",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetName)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub(" genus","",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub(" count data","",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub(" metadata labels","",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub("all scrambled","labels scrambled",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub("all shuffled","samples shuffled",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- factor(mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort,
                                                                            levels = c("HMS ∩ WIS", "HMS scrambled", "HMS shuffled",
                                                                                       "BCM ∩ WIS", "BCM scrambled", "BCM shuffled",
                                                                                       "MDA ∩ WIS", "MDA scrambled", "MDA shuffled",
                                                                                       "WashU ∩ WIS", "WashU scrambled", "WashU shuffled",
                                                                                       "Broad ∩ WIS", "Broad scrambled", "Broad shuffled",
                                                                                       "UNC ∩ WIS", "UNC scrambled", "UNC shuffled",
                                                                                       "CMS ∩ WIS", "CMS scrambled", "CMS shuffled",
                                                                                       "VSNM ∩ WIS","VSNM labels scrambled", "VSNM samples shuffled"))
table(mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)

#----------------------------------Plot overlays----------------------------------#

## Write plot function
source("00-functions.R") # for plotControlsRaw() function
#-------------------------Plot primary tumor overlays-------------------------#
plotControlsRaw(seqCenter = "HMS", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "BCM", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "MDA", inputSampleType = "Primary Tumor", qvalSize = 2.5, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "WashU", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "Broad", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "UNC", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "CMS", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
#-------------------------Plot primary tumor vs NAT overlays-------------------------#
plotControlsRaw(seqCenter = "HMS", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                statSpacingPR = 1.8, qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "BCM", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "Broad", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                statSpacingPR = 1, qvalSize = 3, qvalAsterisks = TRUE)
## NOTE: Neither MDA nor WashU had sufficient samples to plot primary tumor vs. NAT performance
plotControlsRaw(seqCenter = "UNC", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                statSpacingPR = 6, statSpacingROC = 0.9, qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "CMS", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                statSpacingPR = 1.8, qvalSize = 3, qvalAsterisks = TRUE)
#-------------------------Plot blood derived normal overlays-------------------------#
plotControlsRaw(seqCenter = "HMS", inputSampleType = "Blood Derived Normal",
                qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "BCM", inputSampleType = "Blood Derived Normal", 
                statSpacingROC=1, qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "MDA", inputSampleType = "Blood Derived Normal", 
                qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "WashU", inputSampleType = "Blood Derived Normal",
                statSpacingPR = 0.8, statSpacingROC = 0.8, qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "Broad", inputSampleType = "Blood Derived Normal",
                qvalSize = 3, qvalAsterisks = TRUE)
