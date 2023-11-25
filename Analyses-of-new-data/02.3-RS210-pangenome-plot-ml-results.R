#-----------------------------------------------------------------------------
# 02.3-RS210-pangenome-plot-ml-results.R
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

##--------------------Filt5050--------------------##
mlPerf_RS210Pan_Filt5050_Raw_SeqCenter <- read.csv("Supporting_scripts/S05-RS210-Filt5050-seqcenter/rep_perfML_10k_tcga_RS_Filt5050_seqcenter_ALL_7Oct23.csv", 
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
runSeqCenterFXN(inputDataDf = mlPerf_RS210Pan_Filt5050_Raw_SeqCenter,
                prefixRawInput = "rs210PanFinal5050_Nonzero_HiSeq_",
                fileString = "rs210Pan_Raw_Filt5050",
                errorbarPlotFlagInput = TRUE,
                factorCQInput = FALSE)

#----------------------------------------------------------#
# Calc CIs barplots using:
# CQ data
#----------------------------------------------------------#

mlPerf_RS210Pan_CQ <- read.csv("Supporting_scripts/S26-RS210-Pan-CQ/rep_perfML_10k_tcga_RS210Pan_CQ_ALL_16Oct23.csv", 
                            stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

##--------------------Filt5050--------------------##
## WGS
mlPerf_RS210Pan_CQ_Filt5050_WGS <- mlPerf_RS210Pan_CQ %>%
  filter(grepl("5050_",datasetName),
         grepl("WGS",datasetName))

source("00-Functions.R")
runBarplotPerfCQFXN(inputDataDf = mlPerf_RS210Pan_CQ_Filt5050_WGS,
                    prefixRawInput = "rs210PanFinal5050_Nonzero_HiSeq_",
                    errorbarPlotFlagInput = TRUE,
                    fileString = "RS210Pan_CQ_WGS_Filt5050")

## RNA
mlPerf_RS210Pan_CQ_Filt5050_RNA <- mlPerf_RS210Pan_CQ %>%
  filter(grepl("5050_",datasetName),
         grepl("RNA",datasetName))

source("00-Functions.R")
runBarplotPerfCQFXN(inputDataDf = mlPerf_RS210Pan_CQ_Filt5050_RNA,
                    prefixRawInput = "rs210PanFinal5050_Nonzero_HiSeq_",
                    errorbarPlotFlagInput = TRUE,
                    fileString = "RS210Pan_CQ_RNA_Filt5050")

#----------------------------------------------------------#
# Calc CIs barplots using:
# Raw data subsets vs. CQ data subsets
#----------------------------------------------------------#

##--------------------Filt5050--------------------##
mlPerf_RS210Pan_Filt5050_RawVsCQ_SeqCenter <- read.csv("Supporting_scripts/S25-RS210-Pan-RawVsCQ-Filt5050-seqcenter/rep_perfML_10k_tcga_RS210Pan_RawVsCQ_Filt5050_seqcenter_ALL_15Oct23.csv", 
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
runSeqCenterFXN(inputDataDf = mlPerf_RS210Pan_Filt5050_RawVsCQ_SeqCenter,
                prefixRawInput = "rs210PanFinal5050_Nonzero_HiSeq_",
                fileString = "rs210Pan_RawVsCQ_Filt5050",
                factorCQInput = TRUE)

#----------------------------------------------------------#
# Calc enrichment using:
# Raw data subsets vs. CQ data subsets
#----------------------------------------------------------#

##--------------------Filt5050--------------------##
load("Interim_data/RS210_Filt5050_data_for_ml_tcga_by_seq_center_13Oct23.RData")

source("00-Functions.R") # for enrichmentFxn2() function
runSeqCenterEnrichmentFxn2(pathInput = "Supporting_scripts/S25-RS210-Pan-RawVsCQ-Filt5050-seqcenter/features__",
                           totalFeaturesInput = colnames(rs210PanFinal5050_Nonzero_HiSeq_WGS),
                           datasetRawInput = "rs210PanFinal5050_Nonzero_HiSeq",
                           datasetVSNMInput = "rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ",
                           fileNameStringInput = "rs210Pan_Filt5050",
                           kendallOnlyFlagInput = FALSE,
                           combineFlag = TRUE, 
                           combinedWGSonly = TRUE)

#----------------------------------------------------------#
# Plot multiclass ML
#----------------------------------------------------------#

## Add CQ plot data
load("Interim_data/multiclass_CQ_WGS_PT_WIS_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_WGS_BDN_WIS_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_RNA_PT_WIS_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_WGS_PT_Filt9090_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_WGS_BDN_Filt9090_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_RNA_PT_Filt9090_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_WGS_PT_Filt7575_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_WGS_BDN_Filt7575_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_RNA_PT_Filt7575_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_WGS_PT_Filt5050_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_WGS_BDN_Filt5050_13Oct23.RData", verbose = TRUE)
load("Interim_data/multiclass_CQ_RNA_PT_Filt5050_13Oct23.RData", verbose = TRUE)

# CMS only
load("Interim_data/multiclass_CMSonly_RNA_PT_Filt5050_21Oct23.RData", verbose = TRUE)

#-------------------Filt5050-------------------#
plot_cqWGS_Filt5050 <- cmPredHeatmapTCGA(cqWGS_PT_Filt5050$resPredAll, cqWGS_BDN_Filt5050$resPredAll, 
                                    plotString = "RS210Pan_Filt5050_WGS", midPT = 90, midBDN = 70)
plot_cqRNA_Filt5050 <- cmPredHeatmapTCGA(cqRNA_PT_Filt5050$resPredAll, ptOnlyFlag = TRUE, 
                                    plotString = "RS210Pan_Filt5050_RNA", midPT = 700, dimPT = 10)

plot_cqCMS_PT_Filt5050 <- cmPredHeatmapTCGA(cqCMS_PT_Filt5050$resPredAll, ptOnlyFlag = TRUE, 
                                         plotString = "RS210Pan_Filt5050_CMS", midPT = 700, dimPT = 4)

