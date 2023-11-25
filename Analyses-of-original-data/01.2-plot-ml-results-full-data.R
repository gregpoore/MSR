#-----------------------------------------------------------------------------
# 01.2-plot-ml-results-full-data.R
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
# Raw data subsets vs. VSNM data subsets
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals
source("00-functions.R") # for barplotSummaryPerf() and barplotPerf() functions

## With feature intersection
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int <- read.csv("Supporting_scripts/S10-ML-10k-tcga-raw-vsnm-subsets/rep_perfML_10k_tcga_raw_vsnm_subsets_ALL_21Aug23.csv", 
                                                          stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  filter(grepl("_Int_|^snm",datasetName)) %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

## With feature intersection
#----------------All----------------#
source("00-functions.R")
barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Blood Derived Normal",
                   intFlag = TRUE,
                   plotWidthSingle = 8,
                   plotWidthCombined = 8)

barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor",
                   intFlag = TRUE,
                   plotWidthSingle = 12,
                   plotWidthCombined = 12)

barplotSummaryPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
                   seqCenterAbbrev="All",
                   sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                   intFlag = TRUE,
                   plotWidthSingle = 5,
                   plotWidthCombined = 5)

#----------------WGS----------------#
source("00-functions.R")
# HMS
barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="HMS",
            sampleTypeInput = "Blood Derived Normal",
            intFlag = TRUE,
            plotWidthSingle = 5,
            singlePlotHeight = 3,
            plotWidthCombined = 8)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="HMS",
            sampleTypeInput = "Primary Tumor",
            intFlag = TRUE,
            plotWidthSingle = 4.5,
            singlePlotHeight = 3,
            plotWidthCombined = 8)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="HMS",
            sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
            intFlag = TRUE,
            plotWidthSingle = 2.5,
            singlePlotHeight = 3,
            plotWidthCombined = 5)

# BCM
source("00-functions.R")
barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="BCM",
            sampleTypeInput = "Blood Derived Normal",
            intFlag = TRUE,
            plotWidthSingle = 4,
            singlePlotHeight = 3,
            plotWidthCombined = 6)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="BCM",
            sampleTypeInput = "Primary Tumor",
            intFlag = TRUE,
            plotWidthSingle = 4,
            singlePlotHeight = 3,
            plotWidthCombined = 6)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="BCM",
            sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
            intFlag = TRUE,
            plotWidthSingle = 3,
            singlePlotHeight = 3,
            plotWidthCombined = 5)

# MDA
barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="MDA",
            sampleTypeInput = "Blood Derived Normal",
            intFlag = TRUE,
            plotWidthSingle = 4,
            singlePlotHeight = 3,
            plotWidthCombined = 6)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="MDA",
            sampleTypeInput = "Primary Tumor",
            intFlag = TRUE,
            plotWidthSingle = 4,
            singlePlotHeight = 3,
            plotWidthCombined = 6)

# WashU
barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="WashU",
            sampleTypeInput = "Blood Derived Normal",
            intFlag = TRUE,
            plotWidthSingle = 3,
            singlePlotHeight = 3,
            plotWidthCombined = 5)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="WashU",
            sampleTypeInput = "Primary Tumor",
            intFlag = TRUE,
            plotWidthSingle = 3,
            singlePlotHeight = 3,
            plotWidthCombined = 5)

# Broad_WGS
barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="Broad_WGS",
            sampleTypeInput = "Blood Derived Normal",
            intFlag = TRUE,
            plotWidthSingle = 4,
            singlePlotHeight = 3,
            plotWidthCombined = 6)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="Broad_WGS",
            sampleTypeInput = "Primary Tumor",
            intFlag = TRUE,
            plotWidthSingle = 6,
            singlePlotHeight = 3,
            plotWidthCombined = 8)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="Broad_WGS",
            sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
            intFlag = TRUE,
            plotWidthSingle = 2.5,
            singlePlotHeight = 3,
            plotWidthCombined = 5)

#----------------RNA----------------#

# UNC
barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="UNC",
            sampleTypeInput = "Primary Tumor",
            intFlag = TRUE,
            plotWidthSingle = 10,
            singlePlotHeight = 3,
            plotWidthCombined = 18)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="UNC",
            sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
            intFlag = TRUE,
            plotWidthSingle = 6,
            singlePlotHeight = 3,
            plotWidthCombined = 10)

# CMS
barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="CMS",
            sampleTypeInput = "Primary Tumor",
            intFlag = TRUE,
            plotWidthSingle = 3,
            singlePlotHeight = 3,
            plotWidthCombined = 5)

barplotPerf(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
            seqCenterAbbrev="CMS",
            sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
            intFlag = TRUE,
            plotWidthSingle = 2.5,
            singlePlotHeight = 3,
            plotWidthCombined = 3)

#----------------------------------------------------------#
# Calc CIs barplots using:
# Scrambled and shuffled raw data subsets vs. VSNM data subsets
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals
source("00-functions.R") # for barplotSummaryPerf() and barplotPerf() functions

## Scrambled metadata labels
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Scr <- read.csv("Supporting_scripts/S11-ML-10k-tcga-raw-vsnm-subsets-scrambled-controls/rep_perfML_10k_tcga_raw_vsnm_subsets_scrambled_ALL_31Aug23.csv", 
                                                              stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  filter(grepl("_Int_|^snm",datasetName)) %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

runSeqCenterFXN(inputDataDf = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Scr,
                intFlagInput = TRUE,
                fileString = "vsnmVsRaw_SCR")

## Shuffled count data
mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Shu <- read.csv("Supporting_scripts/S12-ML-10k-tcga-raw-vsnm-subsets-shuffled-controls/rep_perfML_10k_tcga_raw_vsnm_subsets_shuffled_ALL_31Aug23.csv", 
                                                              stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  filter(grepl("_Int_|^snm",datasetName)) %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

runSeqCenterFXN(inputDataDf = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Shu,
                intFlagInput = TRUE,
                fileString = "vsnmVsRaw_SHU")

#----------------------------------------------------------#
# Compare VSNM data subsets to scrambled and shuffled raw data subsets
#----------------------------------------------------------#
require(rstatix)

mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Full_Scr_Shu <- rbind(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int,
                                                                    mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Scr,
                                                                    mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Shu)

# save(mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Full_Scr_Shu,
#      file = "Interim_data/mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Full_Scr_Shu_4Nov23.RData")

plotFullScrShu2 <- function(inputData = mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int_Full_Scr_Shu,
                            sampleTypeInput = "Primary Tumor",
                            plot1Width = 10,
                            plot2Width = 3,
                            plot1Height = 5,
                            plot2Height = 5,
                            plot2Start = 0.4,
                            plot2Stop = 1.5,
                            var2Plot = "AUROC",
                            plotPrefix = "auroc_Raw_VSNM_PT_full_scr_shu"){
  
  plotColors <- pal_nejm("default")(6)
  plotColors[1:2] <- c("#ADB6B6FF","#925E9FFF")
  # Swap the purple color since it is already being used for VSNM
  plotColors[which(plotColors=="#7876B1FF")] <- "#FBBC92FF" # "#296695FF"
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
    mutate(datasetName = gsub("vb.+","Raw",datasetName)) %>%
    mutate(datasetName = gsub("snm.+","VSNM",datasetName)) %>%
    select(AUROC, AUPR, abbrev,diseaseType,
           sampleType,datasetName,metadataName) %>%
    mutate(datasetName = factor(case_when(
      grepl("scrambled",metadataName) & (datasetName=="Raw") ~ "Raw_Scr",
      grepl("shuffled",metadataName) & (datasetName=="Raw") ~ "Raw_Shu",
      grepl("scrambled",metadataName) & (datasetName=="VSNM") ~ "VSNM_Scr",
      grepl("shuffled",metadataName) & (datasetName=="VSNM") ~ "VSNM_Shu",
      datasetName=="Raw" ~ "Raw",
      datasetName=="VSNM" ~ "VSNM"),
      levels = c("Raw","VSNM","Raw_Scr","VSNM_Scr","Raw_Shu","VSNM_Shu")
    )) %>%
    mutate(metadataName = "All") %>%
    reshape2::melt(id.vars = c("abbrev","diseaseType","sampleType","datasetName",
                               "metadataName"))
  
  plotData_grouped <- inputData_grouped %>%
    filter(sampleType == sampleTypeInput) %>%
    mutate(datasetName = gsub("vb.+","Raw",datasetName)) %>%
    mutate(datasetName = gsub("snm.+","VSNM",datasetName)) %>%
    select(AUROC, AUPR, abbrev,diseaseType,
           sampleType,datasetName,metadataName) %>%
    mutate(datasetName = factor(case_when(
      grepl("scrambled",metadataName) & (datasetName=="Raw") ~ "Raw_Scr",
      grepl("shuffled",metadataName) & (datasetName=="Raw") ~ "Raw_Shu",
      grepl("scrambled",metadataName) & (datasetName=="VSNM") ~ "VSNM_Scr",
      grepl("shuffled",metadataName) & (datasetName=="VSNM") ~ "VSNM_Shu",
      datasetName=="Raw" ~ "Raw",
      datasetName=="VSNM" ~ "VSNM"),
      levels = c("Raw","VSNM","Raw_Scr","VSNM_Scr","Raw_Shu","VSNM_Shu")
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
                                                         c("Raw","Raw_Scr"),
                                                         c("Raw","Raw_Shu"),
                                                         c("VSNM","VSNM_Scr"),
                                                         c("VSNM","VSNM_Shu"))) %>%
    adjust_pvalue(method = "BH") %>%
    p_round()
  
  plotData_grouped %>%
    filter(variable == var2Plot) %>%
    ggboxplot(x = "datasetName",
              y = "value",
              fill = "datasetName",
              xlab = "",
              ylab = var2Plot,
              legend = "none",
              palette = plotColors) +
    rotate_x_text(30) +
    theme(axis.text.x = element_text(size=10)) +
    scale_y_continuous(breaks = seq(plot2Start, plot2Stop, by = 0.1)) +
    stat_pvalue_manual(plotData_groupedWilcox, label = "p.adj.signif",
                       y.position = c(seq(1.1,1.5,0.4/(nrow(plotData_groupedWilcox)-1) )),
                       size = 3)
  ggsave(filename = paste0("Figures/",plotPrefix,"_narrow_GROUPED.jpeg"),
         units = "in", width = plot2Width, height = plot2Height)
}

#--------------------Grouped--------------------#

# AUROCs
plotFullScrShu2(sampleTypeInput = "Primary Tumor",
                plot1Width = 18,
                plot1Height = 4.5,
                plot2Width = 3,
                plot2Height = 4,
                plot2Start = 0,
                plot2Stop = 1.5,
                plotPrefix = "auroc_Raw_VSNM_PT_full_scr_shu_v2")

plotFullScrShu2(sampleTypeInput = "Blood Derived Normal",
                plot1Width = 14,
                plot1Height = 3.5,
                plot2Width = 3,
                plot2Height = 4,
                plot2Start = 0,
                plot2Stop = 1.5,
                plotPrefix = "auroc_Raw_VSNM_BDN_full_scr_shu_v2")

plotFullScrShu2(sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                plot1Width = 14,
                plot1Height = 3.5,
                plot2Width = 3,
                plot2Height = 4,
                plot2Start = 0,
                plot2Stop = 1.5,
                plotPrefix = "auroc_Raw_VSNM_PTvsSTN_full_scr_shu_v2")

# AUPRs
plotFullScrShu2(sampleTypeInput = "Primary Tumor",
                plot1Width = 18,
                plot1Height = 4.5,
                plot2Width = 3,
                plot2Height = 4,
                plot2Start = 0,
                plot2Stop = 1.5,
                var2Plot = "AUPR",
                plotPrefix = "aupr_Raw_VSNM_PT_full_scr_shu_v2")

plotFullScrShu2(sampleTypeInput = "Blood Derived Normal",
                plot1Width = 14,
                plot1Height = 3.5,
                plot2Width = 3,
                plot2Height = 4,
                plot2Start = 0,
                plot2Stop = 1.5,
                var2Plot = "AUPR",
                plotPrefix = "aupr_Raw_VSNM_BDN_full_scr_shu_v2")

plotFullScrShu2(sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                plot1Width = 14,
                plot1Height = 3.5,
                plot2Width = 3,
                plot2Height = 4,
                plot2Start = 0,
                plot2Stop = 1.5,
                var2Plot = "AUPR",
                plotPrefix = "aupr_Raw_VSNM_PTvsSTN_full_scr_shu_v2")

#----------------------------------------------------------#
# Calc feature overlap / enrichment / rank similarity
# Full data
#----------------------------------------------------------#

## Load total features
load("Input_data/snmDataSampleTypeWithExpStrategyFINAL.RData", verbose = TRUE)

source("00-functions.R") # for enrichmentFxn2() function
## With feature intersection
#----------------WGS----------------#
# HMS
kf_HMS_BDN <- enrichmentFxn2(seqCenter="HMS",
                             sampleType = "Blood Derived Normal",
                             totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                             plotWidth = 5)

kf_HMS_PT <-  enrichmentFxn2(seqCenter="HMS",
                             sampleType = "Primary Tumor",
                             totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                             plotWidth = 5)

kf_HMS_STN <- enrichmentFxn2(seqCenter="HMS",
                             sampleType = "Primary Tumor vs Solid Tissue Normal",
                             totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                             plotWidth = 2.5)

# BCM
kf_BCM_BDN <- enrichmentFxn2(seqCenter="BCM",
                             sampleType = "Blood Derived Normal",
                             totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                             plotWidth = 4)

kf_BCM_PT <- enrichmentFxn2(seqCenter="BCM",
                            sampleType = "Primary Tumor",
                            totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                            plotWidth = 4)

kf_BCM_STN <- enrichmentFxn2(seqCenter="BCM",
                             sampleType = "Primary Tumor vs Solid Tissue Normal",
                             totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                             plotWidth = 3)

# MDA
kf_MDA_BDN <- enrichmentFxn2(seqCenter="MDA",
                             sampleType = "Blood Derived Normal",
                             totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                             plotWidth = 4)

kf_MDA_PT <- enrichmentFxn2(seqCenter="MDA",
                            sampleType = "Primary Tumor",
                            totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                            plotWidth = 4)

# WashU
kf_WashU_BDN <- enrichmentFxn2(seqCenter="WashU",
                               sampleType = "Blood Derived Normal",
                               totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                               plotWidth = 3)

kf_WashU_PT <- enrichmentFxn2(seqCenter="WashU",
                              sampleType = "Primary Tumor",
                              plotWidth = 3)

# Broad_WGS
kf_Broad_WGS_BDN <- enrichmentFxn2(seqCenter="Broad_WGS",
                                   sampleType = "Blood Derived Normal",
                                   totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                                   plotWidth = 4)

kf_Broad_WGS_PT <- enrichmentFxn2(seqCenter="Broad_WGS",
                                  sampleType = "Primary Tumor",
                                  totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                                  plotWidth = 5)

kf_Broad_WGS_STN <- enrichmentFxn2(seqCenter="Broad_WGS",
                                   sampleType = "Primary Tumor vs Solid Tissue Normal",
                                   totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                                   plotWidth = 2.5)

#----------------RNA----------------#

# UNC
kf_UNC_PT <- enrichmentFxn2(seqCenter="UNC",
                            sampleType = "Primary Tumor",
                            totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                            plotWidth = 10)

kf_UNC_STN <- enrichmentFxn2(seqCenter="UNC",
                             sampleType = "Primary Tumor vs Solid Tissue Normal",
                             totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                             plotWidth = 6)

# CMS
kf_CMS_PT <- enrichmentFxn2(seqCenter="CMS",
                            sampleType = "Primary Tumor",
                            totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                            plotWidth = 3)

kf_CMS_STN <- enrichmentFxn2(seqCenter="CMS",
                             sampleType = "Primary Tumor vs Solid Tissue Normal",
                             totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                             plotWidth = 2.5)

#----------------Combine outputs----------------#
kf_Comb <- rbind(# HMS
  kf_HMS_BDN$fisherKendallCombinedDf,
  kf_HMS_PT$fisherKendallCombinedDf,
  kf_HMS_STN$fisherKendallCombinedDf,
  # BCM
  kf_BCM_BDN$fisherKendallCombinedDf,
  kf_BCM_PT$fisherKendallCombinedDf,
  kf_BCM_STN$fisherKendallCombinedDf,
  # MDA
  kf_MDA_BDN$fisherKendallCombinedDf,
  kf_MDA_PT$fisherKendallCombinedDf,
  # WashU
  kf_WashU_BDN$fisherKendallCombinedDf,
  kf_WashU_PT$fisherKendallCombinedDf,
  # Broad_WGS
  kf_Broad_WGS_BDN$fisherKendallCombinedDf,
  kf_Broad_WGS_PT$fisherKendallCombinedDf,
  kf_Broad_WGS_STN$fisherKendallCombinedDf,
  # UNC
  kf_UNC_PT$fisherKendallCombinedDf,
  kf_UNC_STN$fisherKendallCombinedDf,
  # CMS
  kf_CMS_PT$fisherKendallCombinedDf,
  kf_CMS_STN$fisherKendallCombinedDf)

kf_CombP <- kf_Comb %>%
  group_by(abbrev, ST) %>%
  summarise(p.fisher.comb = survcomp::combine.test(p.fisher),
            p.kendall.comb = survcomp::combine.test(p.kendall),
            tau.comb = median(tau),
            tau.se = ifelse(is.na(sd(tau)/n()),0,sd(tau)/n()),
            OR.comb = median(OR),
            OR.se = ifelse(is.na(sd(OR)/n()),0,sd(OR)/n()),
            count = n()) %>%
  rstatix::adjust_pvalue("p.fisher.comb", "p.fisher.comb.adj") %>%
  rstatix::adjust_pvalue("p.kendall.comb", "p.kendall.comb.adj") %>%
  rstatix::add_significance(p.col = "p.fisher.comb.adj") %>%
  rstatix::add_significance(p.col = "p.kendall.comb.adj")

fileNameString <- "rawVSvsnm"
seqCenter <- "All"
plotWidthInput <- c(10,8,5)
sampleTypeInput <- c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal")
for(ii in 1:length(sampleTypeInput)){
  sampleType <- sampleTypeInput[ii]
  plotWidth <- plotWidthInput[ii]
  kf_CombP_Filt <- kf_CombP %>%
    filter(ST == sampleType)
  
  p.adj.signif.fisher.kendall <- c(kf_CombP_Filt$p.fisher.comb.adj.signif,
                                   kf_CombP_Filt$p.kendall.comb.adj.signif)
  p.adj.signif.kendall.fisher <- c(kf_CombP_Filt$p.kendall.comb.adj.signif,
                                   kf_CombP_Filt$p.fisher.comb.adj.signif)
  
  kf_CombP_Filt %>%
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
    ylim(c(0,1.1*max( c(-log10(kf_CombP_Filt$p.fisher.comb.adj),
                        -log10(kf_CombP_Filt$p.kendall.comb.adj)) ))) +
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
  
  kf_CombP_Filt %>%
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
    ylim(c(0,1.1*max( (kf_CombP_Filt$tau.comb) ))) +
    # ylim(c(0,1.1*max( (kf_CombP_Filt$tau.comb+kf_CombP_Filt$tau.se) ))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotKendallTau
  
  print(barPlotKendallTau)
  fileNameTau <- paste0("Figures/kendall_combined_tau_barplot_tcga_",fileNameString,"_",seqCenter,"_",
                        gsub('([[:punct:]])|\\s+','',sampleType),
                        ".jpeg")
  ggsave(filename = fileNameTau,
         plot = barPlotKendallTau,
         dpi = "retina", units = "in", height = 3.5, width = plotWidth)
  
  kf_CombP_Filt %>%
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
    ylim(c(0,1.1*max( (kf_CombP_Filt$OR.comb) ))) +
    # ylim(c(0,1.1*max( (kf_CombP_Filt$OR.comb+kf_CombP_Filt$OR.se) ))) +
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


