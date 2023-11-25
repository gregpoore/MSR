#-----------------------------------------------------------------------------
# 05.4-KU-T2T-plot-ml-results.R
# Copyright (c) 2023--, Greg Poore
# Purposes:
# - Plot machine learning results
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

source("Supporting_scripts/S00-SummarySE.R")
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

# New color scheme:
# - Raw: Gray (#868686FF, #ADB6B6FF)
# - VSNM: Purple (#7876B1FF, #925E9FFF)
# - ConQuR: Orange (#E18727FF)
#----------------------------------------------------------#
# Calc CIs barplots using:
# Raw data subsets
#----------------------------------------------------------#

##--------------------BIO--------------------##
mlPerf_KUT2T_BIO_RawSeqCenter <- read.csv("Supporting_scripts/S12-KU-T2T-BIO-seqcenter/rep_perfML_10k_tcga_KU_T2T_BIO_seqcenter_ALL_10Oct23.csv", 
                                               stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

source("00-Functions.R")
runSeqCenterFXN(inputDataDf = mlPerf_KUT2T_BIO_RawSeqCenter,
                prefixRawInput = "kuT2TFinalNonzero_BIO_HiSeq_",
                fileString = "KUT2T_Raw_BIO",
                errorbarPlotFlagInput = TRUE,
                factorCQInput = TRUE)

#----------------------------------------------------------#
# Calc CIs barplots using:
# CQ data
#----------------------------------------------------------#

mlPerf_KUT2T_CQ <- read.csv("Supporting_scripts/S27-KU-T2T-CQ/rep_perfML_10k_tcga_KUT2T_CQ_ALL_16Oct23.csv", 
                                stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

##--------------------BIO--------------------##
## WGS
mlPerf_KUT2T_CQ_BIO_WGS <- mlPerf_KUT2T_CQ %>%
  filter(grepl("BIO_",datasetName),
         grepl("WGS",datasetName))

source("00-Functions.R")
runBarplotPerfCQFXN(inputDataDf = mlPerf_KUT2T_CQ_BIO_WGS,
                    prefixRawInput = "kuT2TFinalNonzero_BIO_HiSeq_",
                    errorbarPlotFlagInput = TRUE,
                    fileString = "KUT2T_CQ_WGS_BIO")

## RNA
mlPerf_KUT2T_CQ_BIO_RNA <- mlPerf_KUT2T_CQ %>%
  filter(grepl("BIO_",datasetName),
         grepl("RNA",datasetName))

source("00-Functions.R")
runBarplotPerfCQFXN(inputDataDf = mlPerf_KUT2T_CQ_BIO_RNA,
                    prefixRawInput = "kuT2TFinalNonzero_BIO_HiSeq_",
                    errorbarPlotFlagInput = TRUE,
                    fileString = "KUT2T_CQ_RNA_BIO")

#----------------------------------------------------------#
# Calc CIs barplots using:
# Raw data subsets vs. CQ data subsets
#----------------------------------------------------------#

##--------------------BIO--------------------##
mlPerf_KUT2T_BIO_RawVsCQ_SeqCenter <- read.csv("Supporting_scripts/S18-KU-T2T-RawVsCQ-BIO-seqcenter/rep_perfML_10k_tcga_KU_T2T_BIO_RawVsCQ_seqcenter_ALL_15Oct23.csv", 
                                                       stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  filter(!grepl("_UNC|_CMS",datasetName)) %>% # ***filtering out RNA results d/t low counts of UNC (making CMS the only RNA center)***
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

source("00-Functions.R")
runSeqCenterFXN(inputDataDf = mlPerf_KUT2T_BIO_RawVsCQ_SeqCenter,
                prefixRawInput = "kuT2TFinalNonzero_BIO_HiSeq_",
                fileString = "KUT2T_RawVsCQ_BIO",
                factorCQInput = TRUE)

#----------------------------------------------------------#
# Calc CIs barplots using:
# HG38 vs T2T vs Pangenome data subsets
#----------------------------------------------------------#

##--------------------HG38 vs T2T vs Pangenome--------------------##
mlPerf_KUT2T_BIO_All_SeqCenter <- read.csv("Supporting_scripts/S16-KU-All-BIO-seqcenter/rep_perfML_10k_tcga_KU_All_BIO_seqcenter_ALL_15Oct23.csv", 
                                               stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  filter(!grepl("_UNC",datasetName)) %>% # ***filtering out UNC results d/t low counts (CMS is ok)***
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

source("00-Functions.R")
runSeqCenterFXN(inputDataDf = mlPerf_KUT2T_BIO_All_SeqCenter,
                prefixRawInput = "kuAllFinalNonzero_BIO_",
                fileString = "ku_HG38_T2T_Pan",
                factorHD = TRUE,
                factorCQInput = FALSE)

#----------------------------------------------------------#
# Calc enrichment using:
# HG38 vs T2T vs Pangenome data subsets
#----------------------------------------------------------#

##--------------------HG38 vs T2T vs Pangenome--------------------##
load("Interim_data/KU_All_BIO_data_for_ml_tcga_by_seq_center_15Oct23.RData", verbose = TRUE)

# HG38 vs T2T
source("00-Functions.R") # for enrichmentFxn2() function
runSeqCenterEnrichmentFxn2(pathInput = "Supporting_scripts/S16-KU-All-BIO-seqcenter/features__",
                           totalFeaturesInput = colnames(kuAllFinalNonzero_BIO_HG38_HiSeq_WGS),
                           datasetRawInput = "kuAllFinalNonzero_BIO_HG38_HiSeq",
                           datasetVSNMInput = "kuAllFinalNonzero_BIO_T2T_HiSeq",
                           fileNameStringInput = "ku_HG38_vs_T2T",
                           kendallOnlyFlagInput = FALSE,
                           combineFlag = TRUE)

# HG38 vs Pan
source("00-Functions.R") # for enrichmentFxn2() function
runSeqCenterEnrichmentFxn2(pathInput = "Supporting_scripts/S16-KU-All-BIO-seqcenter/features__",
                           totalFeaturesInput = colnames(kuAllFinalNonzero_BIO_HG38_HiSeq_WGS),
                           datasetRawInput = "kuAllFinalNonzero_BIO_HG38_HiSeq",
                           datasetVSNMInput = "kuAllFinalNonzero_BIO_Pan_HiSeq",
                           fileNameStringInput = "ku_HG38_vs_Pan",
                           kendallOnlyFlagInput = FALSE,
                           combineFlag = TRUE)

# T2T vs Pan
source("00-Functions.R") # for enrichmentFxn2() function
runSeqCenterEnrichmentFxn2(pathInput = "Supporting_scripts/S16-KU-All-BIO-seqcenter/features__",
                           totalFeaturesInput = colnames(kuAllFinalNonzero_BIO_T2T_HiSeq_WGS),
                           datasetRawInput = "kuAllFinalNonzero_BIO_T2T_HiSeq",
                           datasetVSNMInput = "kuAllFinalNonzero_BIO_Pan_HiSeq",
                           fileNameStringInput = "ku_T2T_vs_Pan",
                           kendallOnlyFlagInput = FALSE,
                           combineFlag = TRUE)

#----------------------------------------------------------#
# Calc enrichment using:
# Raw data subsets vs. CQ data subsets
# NOTE: If certain ML subsets are missing (eg UNC tumor vs normal) due to sample dropout,
# the function will flag a warning and continue
#----------------------------------------------------------#

##--------------------BIO--------------------##
load("Interim_data/KU_T2T_BIO_data_for_ml_tcga_by_seq_center_14Oct23.RData")

source("00-Functions.R") # for enrichmentFxn2() function
runSeqCenterEnrichmentFxn2(pathInput = "Supporting_scripts/S18-KU-T2T-RawVsCQ-BIO-seqcenter/features__",
                           totalFeaturesInput = colnames(kuT2TFinalNonzero_BIO_HiSeq_WGS),
                           datasetRawInput = "kuT2TFinalNonzero_BIO_HiSeq",
                           datasetVSNMInput = "kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ",
                           fileNameStringInput = "kuT2T_BIO",
                           kendallOnlyFlagInput = FALSE,
                           combineFlag = TRUE,
                           combinedWGSonly = TRUE)

#----------------------------------------------------------#
# Plot multiclass ML
#----------------------------------------------------------#

## Add CQ plot data
load("Interim_data/multiclass_KU_T2T_CQ_WGS_PT_WIS_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_WIS_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_RNA_PT_WIS_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_WGS_PT_BIO_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_BIO_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_RNA_PT_BIO_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_WGS_PT_Filt_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_Filt_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_RNA_PT_Filt_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_WGS_PT_BIOFilt_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_BIOFilt_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_RNA_PT_BIOFilt_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_WGS_PT_Full_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_Full_15Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_KU_T2T_CQ_RNA_PT_Full_15Oct23.RData", verbose = TRUE)

# CMS only
load("Interim_data/multiclass_KU_T2T_CMSonly_RNA_PT_BIO_21Oct23.RData", verbose = TRUE)

#-------------------BIO-------------------#
plot_ku_cqWGS_BIO <- cmPredHeatmap(ku_cqWGS_PT_BIO$resPredAll, ku_cqWGS_BDN_BIO$resPredAll, 
                                         plotString = "KUT2T_BIO_WGS", midPT = 110, midBDN = 80)
plot_ku_cqRNA_BIO <- cmPredHeatmapTCGA(ku_cqRNA_PT_BIO$resPredAll, ptOnlyFlag = TRUE, 
                                         plotString = "KUT2T_BIO_RNA", midPT = 700, dimPT = 10)

plot_ku_CMS_PT_BIO <- cmPredHeatmapTCGA(ku_CMS_PT_BIO$resPredAll, ptOnlyFlag = TRUE, 
                                       plotString = "KUT2T_BIO_CMS", midPT = 700, dimPT = 4)
