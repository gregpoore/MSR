# S14-pvca.R
# Author: Greg Poore
# Date: Oct 10, 2023
# Purpose: Evaluate batch correction efficacy

#-------------------------------#
require(doMC) # for parallel computing
require(lme4)
require(tibble)

numCores <- detectCores()
registerDoMC(cores=numCores)

source("../00-Functions.R") # for PVCA() function

#----------------------------------------------------------#
# RS210-Pan data
#----------------------------------------------------------#

## Load raw datasets
load("../Interim_data/RS210_WIS_data_for_ml_tcga_by_seq_center_7Oct23.RData")
load("../Interim_data/RS210_Filt9090_data_for_ml_tcga_by_seq_center_7Oct23.RData", verbose = TRUE)
load("../Interim_data/RS210_Filt7575_data_for_ml_tcga_by_seq_center_7Oct23.RData")
load("../Interim_data/RS210_Filt5050_data_for_ml_tcga_by_seq_center_7Oct23.RData")

## Load ConQuR-corrected datasets
load("../Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_WGS_CQ_BCM.RData")
load("../Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_RNA_CQ_CMS.RData")
load("../Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("../Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ_CMS.RData")
load("../Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("../Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ_CMS.RData")
load("../Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("../Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ_CMS.RData")

##----------------------------Write function----------------------------##
require(compositions)

pvcaRun <- function(metaData_WGS = metaRSFinalWISNonzero_HiSeq_WGS,
                    countDataRaw_WGS = rs210PanFinalWISNonzero_HiSeq_WGS,
                    countDataCQ_WGS = rs210PanFinalWISNonzero_HiSeq_WGS_CQ,
                    metaData_RNA = NA,
                    countDataRaw_RNA = NA,
                    countDataCQ_RNA = NA,
                    pvcaThreshold = 0.7,
                    wgsOnlyFlag = FALSE){
  
  if(wgsOnlyFlag){
    metaData_WGS_Filt <- metaData_WGS[,c("sample_type",
                                         "disease_type",
                                         "data_submitting_center_label")]
    
    print("Working on raw WGS data...")
    pvca_Raw_WGS <- PVCA(counts = t(countDataRaw_WGS),
                         meta = metaData_WGS_Filt,
                         threshold = pvcaThreshold,
                         inter = FALSE)
    print("Working on CQ WGS data...")
    pvca_CQ_WGS <- PVCA(counts = t(countDataCQ_WGS),
                        meta = metaData_WGS_Filt,
                        threshold = pvcaThreshold,
                        inter = FALSE)
    
    pvcaData <- as.data.frame(rbind(pvca_Raw_WGS,
                                    pvca_CQ_WGS)) %>%
      tibble::rownames_to_column("group")
    print(pvcaData)
    
  } else{
    metaData_WGS_Filt <- metaData_WGS[,c("sample_type",
                                         "disease_type",
                                         "data_submitting_center_label")]
    
    print("Working on raw WGS data...")
    pvca_Raw_WGS <- PVCA(counts = t(countDataRaw_WGS),
                         meta = metaData_WGS_Filt,
                         threshold = pvcaThreshold,
                         inter = FALSE)
    print("Working on CQ WGS data...")
    pvca_CQ_WGS <- PVCA(counts = t(countDataCQ_WGS),
                        meta = metaData_WGS_Filt,
                        threshold = pvcaThreshold,
                        inter = FALSE)
    
    metaData_RNA_Filt <- metaData_RNA[,c("sample_type",
                                         "disease_type",
                                         "data_submitting_center_label")]
    print("Working on raw RNA data...")
    pvca_Raw_RNA <- PVCA(counts = t(countDataRaw_RNA),
                         meta = metaData_RNA_Filt,
                         threshold = pvcaThreshold,
                         inter = FALSE)
    print("Working on CQ RNA data...")
    pvca_CQ_RNA <- PVCA(counts = t(countDataCQ_RNA),
                        meta = metaData_RNA_Filt,
                        threshold = pvcaThreshold,
                        inter = FALSE)
    
    
    pvcaData <- as.data.frame(rbind(pvca_Raw_WGS,
                                    pvca_CQ_WGS,
                                    pvca_Raw_RNA,
                                    pvca_CQ_RNA)) %>%
      tibble::rownames_to_column("group")
    print(pvcaData)
  }
  
  return(pvcaData)
  
}

##----------------------------Apply to RS210-Pan data----------------------------##

pvcaRSPan_WIS <- pvcaRun(metaData_WGS = metaRSFinalWISNonzero_HiSeq_WGS,
                        countDataRaw_WGS = rs210PanFinalWISNonzero_HiSeq_WGS,
                        countDataCQ_WGS = rs210PanFinalWISNonzero_HiSeq_WGS_CQ,
                        metaData_RNA = metaRSFinalWISNonzero_HiSeq_RNA,
                        countDataRaw_RNA = rs210PanFinalWISNonzero_HiSeq_RNA,
                        countDataCQ_RNA = rs210PanFinalWISNonzero_HiSeq_RNA_CQ,
                        pvcaThreshold = 0.8,
                        wgsOnlyFlag = FALSE)
# group           sample_type disease_type data_submitting_center_label     resid
# 1 pvca_Raw_WGS 0.017509069    0.1556433                   0.08254249 0.7443051
# 2  pvca_CQ_WGS 0.001674509    0.1774843                   0.02370120 0.7971400
# 3 pvca_Raw_RNA 0.073382893    0.1305501                   0.14994204 0.6461250
# 4  pvca_CQ_RNA 0.101958160    0.1124170                   0.06678444 0.7188404

pvcaRSPan_Filt9090 <- pvcaRun(metaData_WGS = metaRSFinal9090_Nonzero_HiSeq_WGS,
                        countDataRaw_WGS = rs210PanFinal9090_Nonzero_HiSeq_WGS,
                        countDataCQ_WGS = rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ,
                        metaData_RNA = metaRSFinal9090_Nonzero_HiSeq_RNA,
                        countDataRaw_RNA = rs210PanFinal9090_Nonzero_HiSeq_RNA,
                        countDataCQ_RNA = rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ,
                        pvcaThreshold = 0.8,
                        wgsOnlyFlag = FALSE)
# group           sample_type disease_type data_submitting_center_label     resid
# 1 pvca_Raw_WGS 1.936330e-02   0.28122389                   0.16205966 0.5373532
# 2  pvca_CQ_WGS 2.140615e-03   0.18505813                   0.02298778 0.7898135
# 3 pvca_Raw_RNA 4.124828e-10   0.08598358                   0.01443317 0.8995833
# 4  pvca_CQ_RNA 0.000000e+00   0.06235723                   0.00000000 0.9376428

pvcaRSPan_Filt7575 <- pvcaRun(metaData_WGS = metaRSFinal7575_Nonzero_HiSeq_WGS,
                        countDataRaw_WGS = rs210PanFinal7575_Nonzero_HiSeq_WGS,
                        countDataCQ_WGS = rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ,
                        metaData_RNA = metaRSFinal7575_Nonzero_HiSeq_RNA,
                        countDataRaw_RNA = rs210PanFinal7575_Nonzero_HiSeq_RNA,
                        countDataCQ_RNA = rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ,
                        pvcaThreshold = 0.8,
                        wgsOnlyFlag = FALSE)
# group          sample_type disease_type data_submitting_center_label     resid
# 1 pvca_Raw_WGS 0.013323157   0.26609636                   0.10543822 0.6151423
# 2  pvca_CQ_WGS 0.001767009   0.18703015                   0.02160637 0.7895965
# 3 pvca_Raw_RNA 0.000000000   0.08394740                   0.01942273 0.8966299
# 4  pvca_CQ_RNA 0.000000000   0.06185953                   0.00000000 0.9381405

pvcaRSPan_Filt5050 <- pvcaRun(metaData_WGS = metaRSFinal5050_Nonzero_HiSeq_WGS,
                        countDataRaw_WGS = rs210PanFinal5050_Nonzero_HiSeq_WGS,
                        countDataCQ_WGS = rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ,
                        metaData_RNA = metaRSFinal5050_Nonzero_HiSeq_RNA,
                        countDataRaw_RNA = rs210PanFinal5050_Nonzero_HiSeq_RNA,
                        countDataCQ_RNA = rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ,
                        pvcaThreshold = 0.8,
                        wgsOnlyFlag = FALSE)
# group          sample_type disease_type data_submitting_center_label     resid
# 1 pvca_Raw_WGS 0.013038005   0.26258618                   0.10702286 0.6173530
# 2  pvca_CQ_WGS 0.001726637   0.18643362                   0.02156677 0.7902730
# 3 pvca_Raw_RNA 0.000000000   0.08251285                   0.02378653 0.8937006
# 4  pvca_CQ_RNA 0.000000000   0.06138158                   0.00000000 0.9386184

save(pvcaRSPan_WIS,
     pvcaRSPan_Filt9090,
     pvcaRSPan_Filt7575,
     pvcaRSPan_Filt5050,
     file = "../Interim_data/pvcaRSPan_21Oct23.RData")

##----------------------------Apply to KU-T2T data----------------------------##

## Load raw datasets
load("../Interim_data/KU_T2T_WIS_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../Interim_data/KU_T2T_BIO_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../Interim_data/KU_T2T_Filt_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../Interim_data/KU_T2T_BIOFilt_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../Interim_data/KU_T2T_Full_data_for_ml_tcga_by_seq_center_14Oct23.RData")

## Load ConQuR-corrected datasets
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_BCM.RData")
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ_CMS.RData")
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_BCM.RData")
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ_CMS.RData")
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_BCM.RData")
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ_CMS.RData")
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_BCM.RData")
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ_CMS.RData")
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_BCM.RData")
load("../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_RNA_CQ_CMS.RData")

## WIS
pvcaKUT2T_WIS <- pvcaRun(metaData_WGS = metaKUT2TFinalNonzero_WIS_HiSeq_WGS,
                    countDataRaw_WGS = kuT2TFinalNonzero_WIS_HiSeq_WGS,
                    countDataCQ_WGS = kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ,
                    metaData_RNA = metaKUT2TFinalNonzero_WIS_HiSeq_RNA,
                    countDataRaw_RNA = kuT2TFinalNonzero_WIS_HiSeq_RNA,
                    countDataCQ_RNA = kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ,
                    pvcaThreshold = 0.8,
                    wgsOnlyFlag = FALSE)
# group          sample_type disease_type data_submitting_center_label     resid
# 1 pvca_Raw_WGS 0.024385120   0.19348562                   0.07196838 0.7101609
# 2  pvca_CQ_WGS 0.001625107   0.21314347                   0.05635787 0.7288736
# 3 pvca_Raw_RNA 0.046275089   0.08926960                   0.17397321 0.6904821
# 4  pvca_CQ_RNA 0.681746216   0.02813698                   0.02168271 0.2684341

## BIO
pvcaKUT2T_BIO <- pvcaRun(metaData_WGS = metaKUT2TFinalNonzero_BIO_HiSeq_WGS,
                  countDataRaw_WGS = kuT2TFinalNonzero_BIO_HiSeq_WGS,
                  countDataCQ_WGS = kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ,
                  metaData_RNA = metaKUT2TFinalNonzero_BIO_HiSeq_RNA,
                  countDataRaw_RNA = kuT2TFinalNonzero_BIO_HiSeq_RNA,
                  countDataCQ_RNA = kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ,
                  pvcaThreshold = 0.8,
                  wgsOnlyFlag = FALSE)
# group          sample_type disease_type data_submitting_center_label     resid
# 1 pvca_Raw_WGS 0.013482792   0.18893357                 8.831029e-02 0.7092733
# 2  pvca_CQ_WGS 0.001472363   0.21516080                 5.434995e-02 0.7290169
# 3 pvca_Raw_RNA 0.046340187   0.08291897                 2.243347e-01 0.6464061
# 4  pvca_CQ_RNA 0.685570047   0.03923459                 2.334351e-11 0.2751954


## Filt
pvcaKUT2T_Filt <- pvcaRun(metaData_WGS = metaKUT2TFinalNonzero_Filt_HiSeq_WGS,
                    countDataRaw_WGS = kuT2TFinalNonzero_Filt_HiSeq_WGS,
                    countDataCQ_WGS = kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ,
                    metaData_RNA = metaKUT2TFinalNonzero_Filt_HiSeq_RNA,
                    countDataRaw_RNA = kuT2TFinalNonzero_Filt_HiSeq_RNA,
                    countDataCQ_RNA = kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ,
                    pvcaThreshold = 0.8,
                    wgsOnlyFlag = FALSE)
# group          sample_type disease_type data_submitting_center_label     resid
# 1 pvca_Raw_WGS  0.01560763    0.1584299                  0.064726861 0.7612356
# 2  pvca_CQ_WGS  0.48196426    0.0886015                  0.036757622 0.3926766
# 3 pvca_Raw_RNA  0.10268246    0.3805463                  0.001221717 0.5155496
# 4  pvca_CQ_RNA  0.08801384    0.2277962                  0.000000000 0.6841900

## BIOFilt
pvcaKUT2T_BIOFilt <- pvcaRun(metaData_WGS = metaKUT2TFinalNonzero_BIOFilt_HiSeq_WGS,
                          countDataRaw_WGS = kuT2TFinalNonzero_BIOFilt_HiSeq_WGS,
                          countDataCQ_WGS = kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ,
                          metaData_RNA = metaKUT2TFinalNonzero_BIOFilt_HiSeq_RNA,
                          countDataRaw_RNA = kuT2TFinalNonzero_BIOFilt_HiSeq_RNA,
                          countDataCQ_RNA = kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ,
                          pvcaThreshold = 0.8,
                          wgsOnlyFlag = FALSE)
# group          sample_type disease_type data_submitting_center_label     resid
# 1 pvca_Raw_WGS  0.01560763    0.1584299                  0.064726861 0.7612356
# 2  pvca_CQ_WGS  0.48196426    0.0886015                  0.036757622 0.3926766
# 3 pvca_Raw_RNA  0.10268246    0.3805463                  0.001221717 0.5155496
# 4  pvca_CQ_RNA  0.08801384    0.2277962                  0.000000000 0.6841900

## Full
pvcaKUT2T_Full <- pvcaRun(metaData_WGS = metaKUT2TFinalNonzero_Full_HiSeq_WGS,
                    countDataRaw_WGS = kuT2TFinalNonzero_Full_HiSeq_WGS,
                    countDataCQ_WGS = kuT2TFinalNonzero_Full_HiSeq_WGS_CQ,
                    metaData_RNA = metaKUT2TFinalNonzero_Full_HiSeq_RNA,
                    countDataRaw_RNA = kuT2TFinalNonzero_Full_HiSeq_RNA,
                    countDataCQ_RNA = kuT2TFinalNonzero_Full_HiSeq_RNA_CQ,
                    pvcaThreshold = 0.8,
                    wgsOnlyFlag = FALSE)
# group          sample_type disease_type data_submitting_center_label     resid
# 1 pvca_Raw_WGS 0.014003591   0.19294021                   0.09016826 0.7028879
# 2  pvca_CQ_WGS 0.004688471   0.21740229                   0.05122950 0.7266797
# 3 pvca_Raw_RNA 0.031769176   0.07645448                   0.26841747 0.6233589
# 4  pvca_CQ_RNA 0.691967455   0.02707709                   0.02983489 0.2511206

save(pvcaKUT2T_WIS,
     pvcaKUT2T_BIO,
     pvcaKUT2T_Filt,
     pvcaKUT2T_BIOFilt,
     pvcaKUT2T_Full,
     file = "../Interim_data/pvcaKUT2T_21Oct23.RData")

#----------------------------------------------------------#
# PVCA plots
#----------------------------------------------------------#

load("../Interim_data/pvcaRSPan_21Oct23.RData")
load("../Interim_data/pvcaKUT2T_21Oct23.RData")

##----------------------------Apply to RS210-Pan data----------------------------##

require(ggpubr)
require(ggsci)
# WGS

pvcaPlot <- function(pvcaData=pvcaRSPan_WIS, 
                     dataString="rsPan_WIS"){
  pvcaData$group <- c("Raw WGS", "ConQuR WGS", "Raw RNA","ConQuR RNA")
  pvcaData.melted <- reshape2::melt(pvcaData, id.vars = "group")
  pvcaData.melted$group <- factor(pvcaData.melted$group, 
                                  levels =  c("Raw WGS", "ConQuR WGS",
                                              "Raw RNA", "ConQuR RNA"))
  
  pvcaP <- ggplot(pvcaData.melted, aes(x = variable, y = value, fill = group)) + 
    geom_bar(aes(fill = group), position = "dodge", stat = "identity") + 
    geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), 
              vjust=-0.25, size = 3) + 
    labs(x = "Technical & Biological Effects",
         y = "Weighted average proportion variance",
         title = "PVCA of batch effect correction") +
    theme_pubr() + ylim(c(0,1)) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.5, 0.8)) + 
    guides(fill = guide_legend(nrow = 2)) +
    scale_x_discrete(labels=c("sample_type" = "ST", 
                              "disease_type" = "CT", 
                              "data_submitting_center_label" = "Seqcenter",
                              "resid" = "Residual")) +
    scale_fill_nejm(name = "Data types", labels = c("Raw WGS", "ConQuR WGS",
                                                    "Raw RNA", "ConQuR RNA"))
  print(pvcaP)
  ggsave(plot = pvcaP,
         filename = paste0("../Figures/pvca_",dataString,".jpeg"), 
         dpi = "retina",
         width = 8, height = 4, units = "in")
}

pvcaPlotNoRNA <- function(pvcaData=pvcaRSPan_WIS, 
                          dataString="rsPan_WIS"){
  pvcaData$group <- c("Raw WGS", "ConQuR WGS")
  pvcaData.melted <- reshape2::melt(pvcaData, id.vars = "group")
  pvcaData.melted$group <- factor(pvcaData.melted$group, 
                                  levels =  c("Raw WGS", "ConQuR WGS"))
  
  pvcaP <- ggplot(pvcaData.melted, aes(x = variable, y = value, fill = group)) + 
    geom_bar(aes(fill = group), position = "dodge", stat = "identity") + 
    geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), 
              vjust=-0.25, size = 3) + 
    labs(x = "Technical & Biological Effects",
         y = "Weighted average proportion variance",
         title = "PVCA of batch effect correction") +
    theme_pubr() + ylim(c(0,1)) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.5, 0.8)) + 
    guides(fill = guide_legend(nrow = 2)) +
    scale_x_discrete(labels=c("sample_type" = "ST", 
                              "disease_type" = "CT", 
                              "data_submitting_center_label" = "Seqcenter",
                              "resid" = "Residual")) +
    theme(text = element_text(size = 14)) +
    rotate_x_text(30) +
    scale_fill_nejm(name = "Data types", labels = c("Raw WGS", "ConQuR WGS"))
  print(pvcaP)
  ggsave(plot = pvcaP,
         filename = paste0("../Figures/pvcaNoRNA_",dataString,".jpeg"), 
         dpi = "retina",
         width = 4, height = 6, units = "in")
}

# RS210Pan
pvcaPlot(pvcaData=pvcaRSPan_WIS, 
         dataString="rsPan_WIS")
pvcaPlot(pvcaData=pvcaRSPan_Filt9090, 
         dataString="rsPan_Filt9090")
pvcaPlot(pvcaData=pvcaRSPan_Filt7575, 
         dataString="rsPan_Filt7575")
pvcaPlot(pvcaData=pvcaRSPan_Filt5050, 
         dataString="rsPan_Filt5050")

pvcaRSPan_Filt5050_noRNA <- pvcaRSPan_Filt5050 %>%
  filter(!grepl("RNA",group)) %>% droplevels()
pvcaPlotNoRNA(pvcaData=pvcaRSPan_Filt5050_noRNA, 
              dataString="RSPan_Filt5050")

# KUT2T
pvcaPlot(pvcaData=pvcaKUT2T_WIS, 
         dataString="KUT2T_WIS")
pvcaPlot(pvcaData=pvcaKUT2T_BIO, 
         dataString="KUT2T_BIO")
pvcaPlot(pvcaData=pvcaKUT2T_Filt, 
         dataString="KUT2T_Filt")
pvcaPlot(pvcaData=pvcaKUT2T_BIOFilt, 
         dataString="KUT2T_BIOFilt")
pvcaPlot(pvcaData=pvcaKUT2T_Full, 
         dataString="KUT2T_Full")

# KUT2T
pvcaKUT2T_BIO_noRNA <- pvcaKUT2T_BIO %>%
  filter(!grepl("RNA",group)) %>% droplevels()
pvcaPlotNoRNA(pvcaData=pvcaKUT2T_BIO_noRNA, 
              dataString="KUT2T_BIO")


# pvcaPlotNoRNA(pvcaData=pvcaKUT2T_WIS, 
#          dataString="KUT2T_WIS")
# pvcaPlotNoRNA(pvcaData=pvcaKUT2T_BIO_noRNA, 
#          dataString="KUT2T_BIO")
# pvcaPlotNoRNA(pvcaData=pvcaKUT2T_Filt, 
#          dataString="KUT2T_Filt")
# pvcaPlotNoRNA(pvcaData=pvcaKUT2T_BIOFilt, 
#          dataString="KUT2T_BIOFilt")
# pvcaPlotNoRNA(pvcaData=pvcaKUT2T_Full, 
#          dataString="KUT2T_Full")
