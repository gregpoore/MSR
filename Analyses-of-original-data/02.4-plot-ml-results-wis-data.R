#-----------------------------------------------------------------------------
# 02.4-plot-ml-results-wis-data.R
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

abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

# New color scheme:
# - Raw: Gray (#868686FF, #ADB6B6FF)
# - VSNM: Purple (#7876B1FF, #925E9FFF)
# - ConQuR: Orange (#E18727FF)
#----------------------------------------------------------#
# Calc CIs barplots using:
# WIS-overlapping data
# Raw data subsets vs. VSNM data subsets vs. ConQuR data subsets
# Plus scrambled and shuffled versions thereof
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals
source("00-functions.R") # for barplotSummaryPerf() and barplotPerf() functions

## Load data
mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter <- read.csv("Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/rep_perfML_10k_tcga_WIS_raw_vsnm_cq_subsets_ALL_1Sept23.csv", 
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
runSeqCenterFXN(inputDataDf = mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter,
                intFlagInput = FALSE,
                fileString = "rawVSvsnmVSconqur",
                factorCQInput = TRUE)

## Scrambled metadata
mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter_Scr <- read.csv("Supporting_scripts/S15-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter-scrambled-controls/rep_perfML_10k_tcga_WIS_raw_vsnm_cq_subsets_scrambled_ALL_1Sept23.csv", 
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
runSeqCenterFXN(inputDataDf = mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter_Scr,
                intFlagInput = FALSE,
                fileString = "rawVSvsnmVSconqur_SCR",
                factorCQInput = TRUE)

## Shuffled count data
mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter_Shu <- read.csv("Supporting_scripts/S16-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter-shuffled-controls/rep_perfML_10k_tcga_WIS_raw_vsnm_cq_subsets_shuffled_ALL_1Sept23.csv", 
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
runSeqCenterFXN(inputDataDf = mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter_Shu,
                intFlagInput = FALSE,
                fileString = "rawVSvsnmVSconqur_SHU",
                factorCQInput = TRUE)

#----------------------------------------------------------#
# Compare VSNM data subsets to scrambled and shuffled raw data subsets
#----------------------------------------------------------#
require(rstatix)

mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter_Full_Scr_Shu <- rbind(mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter,
                                                                       mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter_Scr,
                                                                       mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter_Shu)

plotFullScrShuCQ2 <- function(inputData = mlPerfAll10k_Allcancer_WIS_Raw_VSNM_CQ_SeqCenter_Full_Scr_Shu,
                              sampleTypeInput = "Primary Tumor",
                              plot1Width = 10,
                              plot2Width = 3,
                              plot1Height = 5,
                              plot2Height = 5,
                              plot2Start = 0.5,
                              plot2Stop = 1.5,
                              var2Plot = "AUROC",
                              plotPrefix = "auroc_Raw_VSNM_CQ_PT_full_scr_shu"){
  source("Other_scripts/ggsci_adaptive_colors.R")
  # Running the above script imports functions to interpolate an arbitrary number
  # of colors pulling from ggsci's palettes. In this case, we have 9 colors,
  # which can be observed by running: show_col(pal_nejm_adaptive("default")(9))
  # We want to keep the first three colors (Raw, VSNM, ConQuR; in that order) consistent,
  # so we will load the 9 colors and then substitute the first 3 values.
  
  plotColors <- pal_nejm_adaptive("default")(9)
  plotColors[1:3] <- c("#ADB6B6FF","#925E9FFF","#E18727FF")
  # Swap the purple color since it is already being used for VSNM
  plotColors[which(plotColors=="#7483AFFF")] <- "#FBBC92FF" # "#296695FF"
  # print(plotColors)
  
  inputData_grouped <- inputData %>%
    group_by(datasetName, minorityClassName, diseaseType, sampleType, abbrev, metadataName) %>%
    summarize(AUROC=mean(AUROC), 
              AUPR=mean(AUPR),
              nullAUPR=mean(nullAUPR), 
              nullAUROC=mean(nullAUROC)) %>%
    ungroup()
  
  plotData <- inputData %>%
    filter(sampleType == sampleTypeInput) %>%
    mutate(datasetName = case_when(
      grepl("_CQ_",datasetName) ~ "ConQuR",
      grepl("_VSNM_",datasetName) ~ "VSNM",
      (!grepl("_VSNM_|_CQ_",datasetName)) ~ "Raw"
    )) %>%
    mutate(datasetName = gsub("vb.+","Raw",datasetName)) %>%
    mutate(datasetName = gsub("snm.+","VSNM",datasetName)) %>%
    select(AUROC, AUPR, abbrev,diseaseType,
           sampleType,datasetName,metadataName) %>%
    mutate(datasetName = factor(case_when(
      grepl("scrambled",metadataName) & (datasetName=="Raw") ~ "Raw_Scr",
      grepl("shuffled",metadataName) & (datasetName=="Raw") ~ "Raw_Shu",
      grepl("scrambled",metadataName) & (datasetName=="VSNM") ~ "VSNM_Scr",
      grepl("shuffled",metadataName) & (datasetName=="VSNM") ~ "VSNM_Shu",
      grepl("scrambled",metadataName) & (datasetName=="ConQuR") ~ "ConQuR_Scr",
      grepl("shuffled",metadataName) & (datasetName=="ConQuR") ~ "ConQuR_Shu",
      datasetName=="Raw" ~ "Raw",
      datasetName=="VSNM" ~ "VSNM",
      datasetName=="ConQuR" ~ "ConQuR"),
      levels = c("Raw","VSNM","ConQuR",
                 "Raw_Scr","VSNM_Scr","ConQuR_Scr",
                 "Raw_Shu","VSNM_Shu","ConQuR_Shu")
    )) %>%
    mutate(metadataName = "All") %>%
    reshape2::melt(id.vars = c("abbrev","diseaseType","sampleType","datasetName",
                               "metadataName"))
  
  plotData_grouped <- inputData_grouped %>%
    filter(sampleType == sampleTypeInput) %>%
    mutate(datasetName = case_when(
      grepl("_CQ_",datasetName) ~ "ConQuR",
      grepl("_VSNM_",datasetName) ~ "VSNM",
      (!grepl("_VSNM_|_CQ_",datasetName)) ~ "Raw"
    )) %>%
    mutate(datasetName = gsub("vb.+","Raw",datasetName)) %>%
    mutate(datasetName = gsub("snm.+","VSNM",datasetName)) %>%
    select(AUROC, AUPR, abbrev,diseaseType,
           sampleType,datasetName,metadataName) %>%
    mutate(datasetName = factor(case_when(
      grepl("scrambled",metadataName) & (datasetName=="Raw") ~ "Raw_Scr",
      grepl("shuffled",metadataName) & (datasetName=="Raw") ~ "Raw_Shu",
      grepl("scrambled",metadataName) & (datasetName=="VSNM") ~ "VSNM_Scr",
      grepl("shuffled",metadataName) & (datasetName=="VSNM") ~ "VSNM_Shu",
      grepl("scrambled",metadataName) & (datasetName=="ConQuR") ~ "ConQuR_Scr",
      grepl("shuffled",metadataName) & (datasetName=="ConQuR") ~ "ConQuR_Shu",
      datasetName=="Raw" ~ "Raw",
      datasetName=="VSNM" ~ "VSNM",
      datasetName=="ConQuR" ~ "ConQuR"),
      levels = c("Raw","VSNM","ConQuR",
                 "Raw_Scr","VSNM_Scr","ConQuR_Scr",
                 "Raw_Shu","VSNM_Shu","ConQuR_Shu")
    )) %>%
    mutate(metadataName = "All") %>%
    reshape2::melt(id.vars = c("abbrev","diseaseType","sampleType","datasetName",
                               "metadataName"))
  
  # Barplot
  plotData %>%
    summarySE(measurevar = "value", 
              groupvars = c("datasetName", "metadataName","variable","abbrev"),
              conf.interval = 0.99) %>%
    filter(variable == var2Plot) %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=0,
                  position=position_dodge(.9)) +
    scale_fill_manual(values = plotColors, name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = var2Plot,
         title = "") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right")
  ggsave(filename = paste0("Figures/",plotPrefix,"_wide.jpeg"),
         units = "in", width = plot1Width, height = plot1Height)
  
  # Boxplot
  plotData_groupedWilcox <- plotData_grouped %>%
    filter(variable == var2Plot) %>%
    wilcox_test(value ~ datasetName, comparisons = list( c("Raw","VSNM"),
                                                         c("Raw","ConQuR"),
                                                         c("VSNM","ConQuR"),
                                                         c("Raw","Raw_Scr"),
                                                         c("Raw","Raw_Shu"),
                                                         c("VSNM","VSNM_Scr"),
                                                         c("VSNM","VSNM_Shu"),
                                                         c("ConQuR","ConQuR_Scr"),
                                                         c("ConQuR","ConQuR_Shu"))) %>%
    adjust_pvalue(method = "BH") %>%
    p_round()
  
  plotData_grouped %>%
    filter(variable == var2Plot) %>%
    ggboxplot(x = "datasetName",
              y = "value",
              fill = "datasetName",
              xlab = "",
              ylab = var2Plot,
              legend = "none") +
    rotate_x_text(30) +
    scale_fill_manual(values = plotColors) +
    theme(axis.text.x = element_text(size=8)) +
    scale_y_continuous(breaks = seq(plot2Start, plot2Stop, by = 0.1)) +
    stat_pvalue_manual(plotData_groupedWilcox, label = "p.adj.signif",
                       y.position = c(seq(1.05,1.5,0.45/(nrow(plotData_groupedWilcox)-1) )),
                       size = 2, tip.length = 0.01)
  ggsave(filename = paste0("Figures/",plotPrefix,"_narrow_GROUPED.jpeg"),
         units = "in", width = plot2Width, height = plot2Height)
}

#--------------------Grouped--------------------#

# AUROCs
plotFullScrShuCQ2(sampleTypeInput = "Primary Tumor",
                  plot1Width = 18,
                  plot1Height = 4.5,
                  plot2Width = 3,
                  plot2Height = 4,
                  plot2Start = 0.4,
                  plot2Stop = 1.5,
                  plotPrefix = "auroc_Raw_VSNM_CQ_PT_full_scr_shu_v2")

plotFullScrShuCQ2(sampleTypeInput = "Blood Derived Normal",
                  plot1Width = 14,
                  plot1Height = 3.5,
                  plot2Width = 3,
                  plot2Height = 4,
                  plot2Start = 0.4,
                  plot2Stop = 1.5,
                  plotPrefix = "auroc_Raw_VSNM_CQ_BDN_full_scr_shu_v2")

plotFullScrShuCQ2(sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                  plot1Width = 14,
                  plot1Height = 3.5,
                  plot2Width = 3,
                  plot2Height = 4,
                  plot2Start = 0.4,
                  plot2Stop = 1.5,
                  plotPrefix = "auroc_Raw_VSNM_CQ_PTvsSTN_full_scr_shu_v2")

# AUPRs
plotFullScrShuCQ2(sampleTypeInput = "Primary Tumor",
                  plot1Width = 18,
                  plot1Height = 4.5,
                  plot2Width = 3,
                  plot2Height = 4,
                  plot2Start = 0,
                  plot2Stop = 1.5,
                  var2Plot = "AUPR",
                  plotPrefix = "aupr_Raw_VSNM_CQ_PT_full_scr_shu_v2")

plotFullScrShuCQ2(sampleTypeInput = "Blood Derived Normal",
                  plot1Width = 14,
                  plot1Height = 3.5,
                  plot2Width = 3,
                  plot2Height = 4,
                  plot2Start = 0,
                  plot2Stop = 1.5,
                  var2Plot = "AUPR",
                  plotPrefix = "aupr_Raw_VSNM_CQ_BDN_full_scr_shu_v2")

plotFullScrShuCQ2(sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                  plot1Width = 14,
                  plot1Height = 3.5,
                  plot2Width = 3,
                  plot2Height = 4,
                  plot2Start = 0,
                  plot2Stop = 1.5,
                  var2Plot = "AUPR",
                  plotPrefix = "aupr_Raw_VSNM_CQ_PTvsSTN_full_scr_shu_v2")

#----------------------------------------------------------#
# Calc CIs barplots using:
# WIS-overlapping *batch corrected* data
# VSNM data subsets vs. ConQuR data subsets
# Plus scrambled and shuffled versions thereof
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals
source("00-functions.R") # for barplotSummaryPerf() and barplotPerf() functions

## Load data
mlPerfAll10k_Allcancer_WIS_VSNM_CQ <- read.csv("Supporting_scripts/S13-ML-10k-tcga-vsnm-cbs-cq-and-controls/rep_perfML_10k_tcga_WGS_RNA_cq_cbs_vsnm_ALL_21Aug23.csv", 
                                               stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  filter(grepl("_CQ|_VSNM",datasetName)) %>%
  filter(!grepl("_scrambled|_shuffled",metadataName)) %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

source("00-functions.R")
#----------------All----------------#
# Note: The results are *not* split by seqcenters, so only use summary plot FXN
barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_WIS_VSNM_CQ,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Blood Derived Normal",
                   intFlag = FALSE,
                   plotWidthSingle = 8,
                   plotWidthCombined = 8,
                   fileNameString = "WIS_vsnmVSconqur")

barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_WIS_VSNM_CQ,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor",
                   intFlag = FALSE,
                   plotWidthSingle = 12,
                   plotWidthCombined = 12,
                   fileNameString = "WIS_vsnmVSconqur")

barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_WIS_VSNM_CQ,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                   intFlag = FALSE,
                   plotWidthSingle = 6,
                   plotWidthCombined = 6,
                   fileNameString = "WIS_vsnmVSconqur")

## Scrambled metadata
mlPerfAll10k_Allcancer_WIS_VSNM_CQ_Scr <- read.csv("Supporting_scripts/S13-ML-10k-tcga-vsnm-cbs-cq-and-controls/rep_perfML_10k_tcga_WGS_RNA_cq_cbs_vsnm_ALL_21Aug23.csv", 
                                                   stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  filter(grepl("_CQ|_VSNM",datasetName)) %>%
  filter(grepl("_scrambled",metadataName)) %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

source("00-functions.R")
#----------------All----------------#
# Note: The results are *not* split by seqcenters, so only use summary plot FXN
barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_WIS_VSNM_CQ_Scr,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Blood Derived Normal",
                   intFlag = FALSE,
                   plotWidthSingle = 8,
                   plotWidthCombined = 8,
                   fileNameString = "WIS_vsnmVSconqur_SCR")

barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_WIS_VSNM_CQ_Scr,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor",
                   intFlag = FALSE,
                   plotWidthSingle = 12,
                   plotWidthCombined = 12,
                   fileNameString = "WIS_vsnmVSconqur_SCR")

barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_WIS_VSNM_CQ_Scr,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                   intFlag = FALSE,
                   plotWidthSingle = 6,
                   plotWidthCombined = 6,
                   fileNameString = "WIS_vsnmVSconqur_SCR")

## Shuffled count data
mlPerfAll10k_Allcancer_WIS_VSNM_CQ_Shu <- read.csv("Supporting_scripts/S13-ML-10k-tcga-vsnm-cbs-cq-and-controls/rep_perfML_10k_tcga_WGS_RNA_cq_cbs_vsnm_ALL_21Aug23.csv", 
                                                   stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  filter(grepl("_CQ|_VSNM",datasetName)) %>%
  filter(grepl("_shuffled",metadataName)) %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

source("00-functions.R")
#----------------All----------------#
# Note: The results are *not* split by seqcenters, so only use summary plot FXN
barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_WIS_VSNM_CQ_Shu,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Blood Derived Normal",
                   intFlag = FALSE,
                   plotWidthSingle = 8,
                   plotWidthCombined = 8,
                   fileNameString = "WIS_vsnmVSconqur_SHU")

barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_WIS_VSNM_CQ_Shu,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor",
                   intFlag = FALSE,
                   plotWidthSingle = 12,
                   plotWidthCombined = 12,
                   fileNameString = "WIS_vsnmVSconqur_SHU")

barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_WIS_VSNM_CQ_Shu,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                   intFlag = FALSE,
                   plotWidthSingle = 6,
                   plotWidthCombined = 6,
                   fileNameString = "WIS_vsnmVSconqur_SHU")

#----------------------------------------------------------#
# Calc feature overlap / enrichment / rank similarity
# WIS-overlapping data
#----------------------------------------------------------#

## Load total features
load("Interim_data/data_vsnm_tcga_wis_features_subset_25July22.RData", verbose = TRUE)

source("00-functions.R") # for enrichmentFxn2() function
## With feature intersection
dataSetNormWGS <- "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ" # can swap with VSNM
dataSetNormRNA <- "tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ" # can swap with VSNM
fileNameStringInput <- "WIS_rawVSconqur" # can swap with VSNM
kendallOnlyFlagInput <- FALSE
#----------------WGS----------------#
# HMS
kf_WIS_HMS_BDN <- enrichmentFxn2(seqCenter="HMS",
                                 sampleType = "Blood Derived Normal",
                                 plotWidth = 5,
                                 myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                 totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                 dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                 dataSetVSNM = dataSetNormWGS,
                                 fileNameString = fileNameStringInput,
                                 kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_HMS_PT <- enrichmentFxn2(seqCenter="HMS",
                                sampleType = "Primary Tumor",
                                plotWidth = 5,
                                myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                dataSetVSNM = dataSetNormWGS,
                                fileNameString = fileNameStringInput,
                                kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_HMS_STN <- enrichmentFxn2(seqCenter="HMS",
                                 sampleType = "Primary Tumor vs Solid Tissue Normal",
                                 plotWidth = 2,
                                 myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                 totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                 dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                 dataSetVSNM = dataSetNormWGS,
                                 fileNameString = fileNameStringInput,
                                 kendallOnlyFlag = kendallOnlyFlagInput)

# BCM
source("00-functions.R")
kf_WIS_BCM_BDN <- enrichmentFxn2(seqCenter="BCM",
                                 showCM = TRUE,
                                 sampleType = "Blood Derived Normal",
                                 plotWidth = 4,
                                 myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                 totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                 dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                 dataSetVSNM = dataSetNormWGS,
                                 fileNameString = fileNameStringInput,
                                 kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_BCM_PT <- enrichmentFxn2(seqCenter="BCM",
                                sampleType = "Primary Tumor",
                                plotWidth = 4,
                                myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                dataSetVSNM = dataSetNormWGS,
                                fileNameString = fileNameStringInput,
                                kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_BCM_STN <- enrichmentFxn2(seqCenter="BCM",
                                 sampleType = "Primary Tumor vs Solid Tissue Normal",
                                 plotWidth = 2,
                                 myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                 totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                 dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                 dataSetVSNM = dataSetNormWGS,
                                 fileNameString = fileNameStringInput,
                                 kendallOnlyFlag = kendallOnlyFlagInput)

# MDA
kf_WIS_MDA_BDN <- enrichmentFxn2(seqCenter="MDA",
                                 sampleType = "Blood Derived Normal",
                                 plotWidth = 4,
                                 myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                 totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                 dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                 dataSetVSNM = dataSetNormWGS,
                                 fileNameString = fileNameStringInput,
                                 kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_MDA_PT <- enrichmentFxn2(seqCenter="MDA",
                                sampleType = "Primary Tumor",
                                plotWidth = 4,
                                myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                dataSetVSNM = dataSetNormWGS,
                                fileNameString = fileNameStringInput,
                                kendallOnlyFlag = kendallOnlyFlagInput)

# WashU
kf_WIS_WashU_BDN <- enrichmentFxn2(seqCenter="WashU",
                                   sampleType = "Blood Derived Normal",
                                   plotWidth = 3,
                                   myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                   totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                   dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                   dataSetVSNM = dataSetNormWGS,
                                   fileNameString = fileNameStringInput,
                                   kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_WashU_PT <- enrichmentFxn2(seqCenter="WashU",
                                  sampleType = "Primary Tumor",
                                  plotWidth = 3,
                                  myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                  totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                  dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                  dataSetVSNM = dataSetNormWGS,
                                  fileNameString = fileNameStringInput,
                                  kendallOnlyFlag = kendallOnlyFlagInput)

# Broad_WGS
kf_WIS_Broad_WGS_BDN <- enrichmentFxn2(seqCenter="Broad_WGS",
                                       sampleType = "Blood Derived Normal",
                                       plotWidth = 4,
                                       myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                       totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                       dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                       dataSetVSNM = dataSetNormWGS,
                                       fileNameString = fileNameStringInput,
                                       kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_Broad_WGS_PT <- enrichmentFxn2(seqCenter="Broad_WGS",
                                      sampleType = "Primary Tumor",
                                      plotWidth = 5,
                                      myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                      totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                      dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                      dataSetVSNM = dataSetNormWGS,
                                      fileNameString = fileNameStringInput,
                                      kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_Broad_WGS_STN <- enrichmentFxn2(seqCenter="Broad_WGS",
                                       sampleType = "Primary Tumor vs Solid Tissue Normal",
                                       plotWidth = 2,
                                       myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                       totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                       dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_WGS",
                                       dataSetVSNM = dataSetNormWGS,
                                       fileNameString = fileNameStringInput,
                                       kendallOnlyFlag = kendallOnlyFlagInput)

#----------------RNA----------------#

# UNC
kf_WIS_UNC_PT <- enrichmentFxn2(seqCenter="UNC",
                                sampleType = "Primary Tumor",
                                plotWidth = 10,
                                myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_RNA",
                                dataSetVSNM = dataSetNormRNA,
                                fileNameString = fileNameStringInput,
                                kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_UNC_STN <- enrichmentFxn2(seqCenter="UNC",
                                 sampleType = "Primary Tumor vs Solid Tissue Normal",
                                 plotWidth = 6,
                                 myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                 totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                 dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_RNA",
                                 dataSetVSNM = dataSetNormRNA,
                                 fileNameString = fileNameStringInput,
                                 kendallOnlyFlag = kendallOnlyFlagInput)

# CMS
kf_WIS_CMS_PT <- enrichmentFxn2(seqCenter="CMS",
                                sampleType = "Primary Tumor",
                                plotWidth = 3,
                                myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_RNA",
                                dataSetVSNM = dataSetNormRNA,
                                fileNameString = fileNameStringInput,
                                kendallOnlyFlag = kendallOnlyFlagInput)

kf_WIS_CMS_STN <- enrichmentFxn2(seqCenter="CMS",
                                 sampleType = "Primary Tumor vs Solid Tissue Normal",
                                 plotWidth = 2,
                                 myPath = "Supporting_scripts/S14-ML-10k-tcga-WIS-raw-vsnm-cq-seqcenter/features__",
                                 totalFeatures = colnames(vsnmDataGenusKrakenQCFiltWIS),
                                 dataSetRaw = "tcgaGenusKrakenAllFiltWIS_HiSeq_RNA",
                                 dataSetVSNM = dataSetNormRNA,
                                 fileNameString = fileNameStringInput,
                                 kendallOnlyFlag = kendallOnlyFlagInput)

#----------------Combine outputs----------------#
kf_WIS_Comb <- rbind(# HMS
  kf_WIS_HMS_BDN$fisherKendallCombinedDf,
  kf_WIS_HMS_PT$fisherKendallCombinedDf,
  kf_WIS_HMS_STN$fisherKendallCombinedDf,
  # BCM
  kf_WIS_BCM_BDN$fisherKendallCombinedDf,
  kf_WIS_BCM_PT$fisherKendallCombinedDf,
  kf_WIS_BCM_STN$fisherKendallCombinedDf,
  # MDA
  kf_WIS_MDA_BDN$fisherKendallCombinedDf,
  kf_WIS_MDA_PT$fisherKendallCombinedDf,
  # WashU
  kf_WIS_WashU_BDN$fisherKendallCombinedDf,
  kf_WIS_WashU_PT$fisherKendallCombinedDf,
  # Broad_WGS
  kf_WIS_Broad_WGS_BDN$fisherKendallCombinedDf,
  kf_WIS_Broad_WGS_PT$fisherKendallCombinedDf,
  kf_WIS_Broad_WGS_STN$fisherKendallCombinedDf,
  # UNC
  kf_WIS_UNC_PT$fisherKendallCombinedDf,
  kf_WIS_UNC_STN$fisherKendallCombinedDf,
  # CMS
  kf_WIS_CMS_PT$fisherKendallCombinedDf,
  kf_WIS_CMS_STN$fisherKendallCombinedDf)

kf_WIS_CombP <- kf_WIS_Comb %>%
  group_by(abbrev, ST) %>%
  summarise(p.fisher.comb = survcomp::combine.test(p.fisher),
            p.kendall.comb = survcomp::combine.test(p.kendall),
            tau.comb = median(tau),
            tau.se = ifelse(is.na(sd(tau)/n()),0,sd(tau)/n()),
            OR.comb = mean(OR),
            OR.se = ifelse(is.na(sd(OR)/n()),0,sd(OR)/n()),
            count = n()) %>%
  rstatix::adjust_pvalue("p.fisher.comb", "p.fisher.comb.adj") %>%
  rstatix::adjust_pvalue("p.kendall.comb", "p.kendall.comb.adj") %>%
  rstatix::add_significance(p.col = "p.fisher.comb.adj") %>%
  rstatix::add_significance(p.col = "p.kendall.comb.adj")

fileNameString <- "WIS_rawVSconqur" # "WIS_rawVSvsnm" # "WIS_rawVSconqur"
seqCenter <- "All"
plotWidthInput <- c(10,8,5)
sampleTypeInput <- c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal")

# Kendall and Fisher
for(ii in 1:length(sampleTypeInput)){
  sampleType <- sampleTypeInput[ii]
  plotWidth <- plotWidthInput[ii]
  kf_WIS_CombP_Filt <- kf_WIS_CombP %>%
    filter(ST == sampleType)
  
  p.adj.signif.fisher.kendall <- c(kf_WIS_CombP_Filt$p.fisher.comb.adj.signif,
                                   kf_WIS_CombP_Filt$p.kendall.comb.adj.signif)
  p.adj.signif.kendall.fisher <- c(kf_WIS_CombP_Filt$p.kendall.comb.adj.signif,
                                   kf_WIS_CombP_Filt$p.fisher.comb.adj.signif)
  
  kf_WIS_CombP_Filt %>%
    mutate(log.p.fisher.comb.adj = -log10(p.fisher.comb.adj),
           log.p.kendall.comb.adj = -log10(p.kendall.comb.adj)) %>%
    select(abbrev, log.p.fisher.comb.adj, log.p.kendall.comb.adj,
           p.fisher.comb.adj.signif, p.kendall.comb.adj.signif, count) %>%
    reshape2::melt(id.vars = c("abbrev","p.fisher.comb.adj.signif","p.kendall.comb.adj.signif","count")) %>%
    mutate(variable = factor(case_when(
      variable == "log.p.fisher.comb.adj" ~ "Fisher",
      variable == "log.p.kendall.comb.adj" ~ "Kendall",
    ), levels = c("Kendall","Fisher"))) %>%
    mutate(p.adj.signif.combined = `if`(head(as.character(variable),1)=="Fisher",
                                        p.adj.signif.fisher.kendall,
                                        p.adj.signif.kendall.fisher)) %>%
    ggplot(aes(x = reorder(abbrev,-value,median), y = value, fill=variable)) +
    geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
    geom_text(aes(label = p.adj.signif.combined, y = value),
              position = position_dodge(0.9), 
              hjust = -0.2,
              vjust = 0.8,
              size = 3,
              angle = 90) +
    geom_text(aes(label = count, y = value),
              position = position_dodge(0.9), 
              hjust = 0.5,
              vjust = 1.5,
              size = 2,
              angle = 0,
              color = "white") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
    ylim(c(0,1.1*max( c(-log10(kf_WIS_CombP_Filt$p.fisher.comb.adj),
                        -log10(kf_WIS_CombP_Filt$p.kendall.comb.adj)) ))) +
    scale_fill_nejm(name = "Test type") +
    theme_pubr() +
    rotate_x_text(30) +
    theme(axis.text.x = element_text(size=10)) +
    labs(x = "TCGA Cancer Type", 
         y = "-Log(p-adjust)",
         title = paste(seqCenter,sampleType,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> pvalPlot
  
  print(pvalPlot)
  fileNamePval <- paste0("Figures/pvalue_combined_Fisher_Kendall_barplot_tcga_",fileNameString,"_",seqCenter,"_",
                         gsub('([[:punct:]])|\\s+','',sampleType),
                         ".jpeg")
  ggsave(filename = fileNamePval,
         plot = pvalPlot,
         dpi = "retina", units = "in", height = 3, width = plotWidth)
  
  kf_WIS_CombP_Filt %>%
    ggplot(aes(x = reorder(abbrev, -tau.comb,median), y = tau.comb)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    # geom_errorbar(aes(ymin=tau.comb-tau.se, ymax=tau.comb+tau.se), width=.2) +
    theme_pubr() +
    geom_text(aes(label = count, y = tau.comb),
              position = position_dodge(0.9), 
              hjust = 0.5,
              vjust = 1.5,
              size = 2,
              angle = 0,
              color = "white") +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "Kendall tau",
         title = paste(seqCenter,sampleType,sep = " | ")) +
    geom_text(aes(label = p.kendall.comb.adj.signif, y = tau.comb), vjust = -0.4) +
    ylim(c(0,1.1*max( (kf_WIS_CombP_Filt$tau.comb) ))) +
    # ylim(c(0,1.1*max( (kf_WIS_CombP_Filt$tau.comb+kf_WIS_CombP_Filt$tau.se) ))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotKendallTau
  
  print(barPlotKendallTau)
  fileNameTau <- paste0("Figures/kendall_combined_tau_barplot_tcga_",fileNameString,"_",seqCenter,"_",
                        gsub('([[:punct:]])|\\s+','',sampleType),
                        ".jpeg")
  ggsave(filename = fileNameTau,
         plot = barPlotKendallTau,
         dpi = "retina", units = "in", height = 3.5, width = plotWidth)
  
  kf_WIS_CombP_Filt %>%
    ggplot(aes(x = reorder(abbrev, -OR.comb,median), y = OR.comb)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    # geom_errorbar(aes(ymin=OR.comb-OR.se, ymax=OR.comb+OR.se), width=.2) +
    theme_pubr() +
    geom_text(aes(label = count, y = OR.comb),
              position = position_dodge(0.9), 
              hjust = 0.5,
              vjust = 1.5,
              size = 2,
              angle = 0,
              color = "white") +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "Odds ratio feature enrichment",
         title = paste(seqCenter,sampleType,sep = " | ")) +
    geom_text(aes(label = p.fisher.comb.adj.signif, y = OR.comb), vjust = -0.4) +
    ylim(c(0,1.1*max( (kf_WIS_CombP_Filt$OR.comb) ))) +
    # ylim(c(0,1.1*max( (kf_WIS_CombP_Filt$OR.comb+kf_WIS_CombP_Filt$OR.se) ))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotOR
  
  print(barPlotOR)
  fileNameOR <- paste0("Figures/odds_ratio_combined_barplot_tcga_",fileNameString,"_",seqCenter,"_",
                       gsub('([[:punct:]])|\\s+','',sampleType),
                       ".jpeg")
  ggsave(filename = fileNameOR,
         plot = barPlotOR,
         dpi = "retina", units = "in", height = 3.5, width = plotWidth)
}


#-----------------------------------------------------------------------------------------------------#
# Format machine learning performances for all TCGA cancers on raw data split by sequencing center
#-----------------------------------------------------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS <- read.csv("Interim_data/rep_perfML_10k_tcga_wis_feature_subset_seqcenter_vsnm_ALL_25July22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS <- mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS[,!(colnames(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS) == "X")]
colnames(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$minorityClassName == "SolidTissueNormal",
                                                                 yes=mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$majorityClassSize/(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$minorityClassSize+mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$majorityClassSize),
                                                                 no=mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$minorityClassSize/(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$minorityClassSize+mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$majorityClassSize))
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$nullAUROC <- 0.5

# Rename entries in the "datasetName" column

# All (VSNM) - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName[mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName == "vsnmDataGenusKrakenQCFiltWIS"] <- "VSNM genus ∩ WIS (WGS)"
# HMS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName[mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName == "tcgaGenusKrakenQCFiltWIS_HMS"] <- "HMS genus ∩ WIS (WGS)"
# BCM - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName[mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName == "tcgaGenusKrakenQCFiltWIS_BCM"] <- "BCM genus ∩ WIS (WGS)"
# MDA - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName[mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName == "tcgaGenusKrakenQCFiltWIS_MDA"] <- "MDA genus ∩ WIS (WGS)"
# WashU - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName[mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName == "tcgaGenusKrakenQCFiltWIS_WashU"] <- "WashU genus ∩ WIS (WGS)"
# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName[mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName == "tcgaGenusKrakenQCFiltWIS_Broad_WGS"] <- "Broad genus ∩ WIS (WGS)"
# UNC - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName[mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName == "tcgaGenusKrakenQCFiltWIS_UNC"] <- "UNC genus ∩ WIS (RNA-Seq)"
# CMS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName[mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName == "tcgaGenusKrakenQCFiltWIS_CMS"] <- "CMS genus ∩ WIS (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName <- factor(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName,
                                                                    levels = c("VSNM genus ∩ WIS (WGS)",
                                                                               "HMS genus ∩ WIS (WGS)",
                                                                               "BCM genus ∩ WIS (WGS)",
                                                                               "MDA genus ∩ WIS (WGS)",
                                                                               "WashU genus ∩ WIS (WGS)",
                                                                               "Broad genus ∩ WIS (WGS)",
                                                                               "UNC genus ∩ WIS (RNA-Seq)",
                                                                               "CMS genus ∩ WIS (RNA-Seq)"))
table(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS$datasetName)

save(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS,
     file = "Interim_data/mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS_26July22.RData")

#-------------------------------------------#
# Plot ML perf
#-------------------------------------------#

#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data ∩ WIS Features\nPrimary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/mlPerfAll10k_rep1_HMS_PT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# BCM - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/mlPerfAll10k_rep1_BCM_PT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# MDA - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_MDA_PT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# WashU - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_WashU_PT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_Broad_WGS_PT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# UNC - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_UNC_PT.svg", dpi = "retina",
       width = 10, height = 3.5, units = "in")

# CMS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_CMS_PT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/mlPerfAll10k_rep1_HMS_PT_vs_NAT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# BCM - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_BCM_PT_vs_NAT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

## NOTE: Neither MDA nor WashU has enough tumor vs. normal samples to plot

# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# UNC - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_UNC_PT_vs_NAT.svg", dpi = "retina",
       width = 10, height = 3.5, units = "in")

# CMS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_CMS_PT_vs_NAT.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# HMS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_HMS_BDN.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# BCM - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_BCM_BDN.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# MDA - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_MDA_BDN.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# WashU - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_WashU_BDN.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_WIS %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) +
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/mlPerfAll10k_rep1_Broad_WGS_BDN.svg", dpi = "retina",
       width = 4, height = 3.5, units = "in")

#-------------------------------------------#
# Plot ML multiclass pred/perf
#-------------------------------------------#
require(tibble)

## Add CQ plot data
load("Interim_data/multiclass_PT_results_wgs_CQ_28Sept23.RData", verbose = TRUE)
load("Interim_data/multiclass_BDN_results_wgs_CQ_28Sept23.RData", verbose = TRUE)
load("Interim_data/multiclass_PT_results_RNA_CQ_30Sept23.RData", verbose = TRUE)

source("00-functions.R") # for cmPredHeatmapTCGA() function
cqCM_WGS <- cmPredHeatmapTCGA(cq_PT$resPredAll, cq_BDN$resPredAll, plotString = "WGS_CQ", midPT = 90, midBDN = 60)
cqCM_RNA <- cmPredHeatmapTCGA(cqRNA_PT$resPredAll, ptOnlyFlag = TRUE, 
                              plotString = "RNA_CQ", midPT = 700, dimPT = 10)

## Add VSNM (paired to CQ) plot data
load("Interim_data/multiclass_PT_results_wgs_VSNM_28Sept23.RData", verbose = TRUE)
load("Interim_data/multiclass_BDN_results_wgs_VSNM_28Sept23.RData", verbose = TRUE)
load("Interim_data/multiclass_PT_results_RNA_VSNM_30Sept23.RData", verbose = TRUE)

vsnmCM_WGS <- cmPredHeatmapTCGA(vsnmMulti_PT$resPredAll, vsnmMulti_BDN$resPredAll, 
                                plotString = "WGS_VSNM", midPT = 90, midBDN = 60)
vsnmCM_RNA <- cmPredHeatmapTCGA(vsnmMultiRNA_PT$resPredAll, ptOnlyFlag = TRUE, 
                                plotString = "RNA_VSNM", midPT = 700, dimPT = 10)



