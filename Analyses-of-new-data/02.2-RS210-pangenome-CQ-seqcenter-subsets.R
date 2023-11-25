# 02.2-RS210-pangenome-CQ-seqcenter-subsets.R
# Author: Greg Poore
# Date: Oct 15, 2023
# Purposes:
# - Load RS210 data and format

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
# Load seqcenter count data
load("Interim_data/RS210_WIS_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/RS210_Filt9090_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/RS210_Filt7575_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/RS210_Filt5050_data_for_ml_tcga_by_seq_center_13Oct23.RData")

# Load ConQuR-corrected data
load("Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_WGS_CQ_BCM.RData", verbose=T)
load("Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ_CMS.RData")

#----------------------------------------------------------#
# Subset WIS data
#----------------------------------------------------------#

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rs210PanFinalWISNonzero_HiSeq_WGS_CQ_HMS <- rs210PanFinalWISNonzero_HiSeq_WGS_CQ[rownames(metaRSFinalWISNonzero_HiSeq_HMS),]
rs210PanFinalWISNonzero_HiSeq_WGS_CQ_BCM <- rs210PanFinalWISNonzero_HiSeq_WGS_CQ[rownames(metaRSFinalWISNonzero_HiSeq_BCM),]
rs210PanFinalWISNonzero_HiSeq_WGS_CQ_MDA <- rs210PanFinalWISNonzero_HiSeq_WGS_CQ[rownames(metaRSFinalWISNonzero_HiSeq_MDA),]
rs210PanFinalWISNonzero_HiSeq_WGS_CQ_WashU <- rs210PanFinalWISNonzero_HiSeq_WGS_CQ[rownames(metaRSFinalWISNonzero_HiSeq_WashU),]
rs210PanFinalWISNonzero_HiSeq_WGS_CQ_Broad_WGS <- rs210PanFinalWISNonzero_HiSeq_WGS_CQ[rownames(metaRSFinalWISNonzero_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both RNA and RNA-Seq, so separate objects are made for both)
rs210PanFinalWISNonzero_HiSeq_RNA_CQ_UNC <- rs210PanFinalWISNonzero_HiSeq_RNA_CQ[rownames(metaRSFinalWISNonzero_HiSeq_UNC),]
rs210PanFinalWISNonzero_HiSeq_RNA_CQ_CMS <- rs210PanFinalWISNonzero_HiSeq_RNA_CQ[rownames(metaRSFinalWISNonzero_HiSeq_CMS),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  rs210PanFinalWISNonzero_HiSeq_HMS,
  rs210PanFinalWISNonzero_HiSeq_BCM,
  rs210PanFinalWISNonzero_HiSeq_MDA,
  rs210PanFinalWISNonzero_HiSeq_WashU,
  rs210PanFinalWISNonzero_HiSeq_Broad_WGS,
  rs210PanFinalWISNonzero_HiSeq_UNC,
  rs210PanFinalWISNonzero_HiSeq_CMS,
  
  # ConQuR-corrected data
  rs210PanFinalWISNonzero_HiSeq_WGS_CQ_HMS,
  rs210PanFinalWISNonzero_HiSeq_WGS_CQ_BCM,
  rs210PanFinalWISNonzero_HiSeq_WGS_CQ_MDA,
  rs210PanFinalWISNonzero_HiSeq_WGS_CQ_WashU,
  rs210PanFinalWISNonzero_HiSeq_WGS_CQ_Broad_WGS,
  rs210PanFinalWISNonzero_HiSeq_RNA_CQ_UNC,
  rs210PanFinalWISNonzero_HiSeq_RNA_CQ_CMS,
  
  # WGS and RNA subset data
  rs210PanFinalWISNonzero_HiSeq_WGS,
  rs210PanFinalWISNonzero_HiSeq_RNA,
  rs210PanFinalWISNonzero_HiSeq_WGS_CQ,
  rs210PanFinalWISNonzero_HiSeq_RNA_CQ,
  
  # Subset metadata
  metaRSFinalWISNonzero_HiSeq_HMS,
  metaRSFinalWISNonzero_HiSeq_BCM,
  metaRSFinalWISNonzero_HiSeq_MDA,
  metaRSFinalWISNonzero_HiSeq_WashU,
  metaRSFinalWISNonzero_HiSeq_Broad_WGS,
  metaRSFinalWISNonzero_HiSeq_UNC,
  metaRSFinalWISNonzero_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaRSFinalWISNonzero_HiSeq_WGS,
  metaRSFinalWISNonzero_HiSeq_RNA,
  
  # WIS metadata
  metaRSFinalWISNonzero_HiSeq,
  file = "Interim_data/RS210_Pan_WIS_RawVsCQ_data_for_ml_tcga_by_seq_center_15Oct23.RData")

#----------------------------------------------------------#
# Subset Filt5050 data
#----------------------------------------------------------#

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_HMS <- rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal5050_Nonzero_HiSeq_HMS),]
rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_BCM <- rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal5050_Nonzero_HiSeq_BCM),]
rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_MDA <- rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal5050_Nonzero_HiSeq_MDA),]
rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_WashU <- rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal5050_Nonzero_HiSeq_WashU),]
rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_Broad_WGS <- rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal5050_Nonzero_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both RNA and RNA-Seq, so separate objects are made for both)
rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ_UNC <- rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ[rownames(metaRSFinal5050_Nonzero_HiSeq_UNC),]
rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ_CMS <- rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ[rownames(metaRSFinal5050_Nonzero_HiSeq_CMS),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  rs210PanFinal5050_Nonzero_HiSeq_HMS,
  rs210PanFinal5050_Nonzero_HiSeq_BCM,
  rs210PanFinal5050_Nonzero_HiSeq_MDA,
  rs210PanFinal5050_Nonzero_HiSeq_WashU,
  rs210PanFinal5050_Nonzero_HiSeq_Broad_WGS,
  rs210PanFinal5050_Nonzero_HiSeq_UNC,
  rs210PanFinal5050_Nonzero_HiSeq_CMS,
  
  # ConQuR-corrected data
  rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_HMS,
  rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_BCM,
  rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_MDA,
  rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_WashU,
  rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_Broad_WGS,
  rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ_UNC,
  rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ_CMS,
  
  # WGS and RNA subset data
  rs210PanFinal5050_Nonzero_HiSeq_WGS,
  rs210PanFinal5050_Nonzero_HiSeq_RNA,
  rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ,
  rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ,
  
  # Subset metadata
  metaRSFinal5050_Nonzero_HiSeq_HMS,
  metaRSFinal5050_Nonzero_HiSeq_BCM,
  metaRSFinal5050_Nonzero_HiSeq_MDA,
  metaRSFinal5050_Nonzero_HiSeq_WashU,
  metaRSFinal5050_Nonzero_HiSeq_Broad_WGS,
  metaRSFinal5050_Nonzero_HiSeq_UNC,
  metaRSFinal5050_Nonzero_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaRSFinal5050_Nonzero_HiSeq_WGS,
  metaRSFinal5050_Nonzero_HiSeq_RNA,
  
  # Full metadata
  metaRSFinal5050_Nonzero_HiSeq,
  file = "Interim_data/RS210_Pan_Filt5050_RawVsCQ_data_for_ml_tcga_by_seq_center_15Oct23.RData")

#----------------------------------------------------------#
# Subset Filt7575 data
#----------------------------------------------------------#

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_HMS <- rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal7575_Nonzero_HiSeq_HMS),]
rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_BCM <- rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal7575_Nonzero_HiSeq_BCM),]
rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_MDA <- rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal7575_Nonzero_HiSeq_MDA),]
rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_WashU <- rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal7575_Nonzero_HiSeq_WashU),]
rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_Broad_WGS <- rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal7575_Nonzero_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both RNA and RNA-Seq, so separate objects are made for both)
rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ_UNC <- rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ[rownames(metaRSFinal7575_Nonzero_HiSeq_UNC),]
rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ_CMS <- rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ[rownames(metaRSFinal7575_Nonzero_HiSeq_CMS),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  rs210PanFinal7575_Nonzero_HiSeq_HMS,
  rs210PanFinal7575_Nonzero_HiSeq_BCM,
  rs210PanFinal7575_Nonzero_HiSeq_MDA,
  rs210PanFinal7575_Nonzero_HiSeq_WashU,
  rs210PanFinal7575_Nonzero_HiSeq_Broad_WGS,
  rs210PanFinal7575_Nonzero_HiSeq_UNC,
  rs210PanFinal7575_Nonzero_HiSeq_CMS,
  
  # ConQuR-corrected data
  rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_HMS,
  rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_BCM,
  rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_MDA,
  rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_WashU,
  rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_Broad_WGS,
  rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ_UNC,
  rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ_CMS,
  
  # WGS and RNA subset data
  rs210PanFinal7575_Nonzero_HiSeq_WGS,
  rs210PanFinal7575_Nonzero_HiSeq_RNA,
  rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ,
  rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ,
  
  # Subset metadata
  metaRSFinal7575_Nonzero_HiSeq_HMS,
  metaRSFinal7575_Nonzero_HiSeq_BCM,
  metaRSFinal7575_Nonzero_HiSeq_MDA,
  metaRSFinal7575_Nonzero_HiSeq_WashU,
  metaRSFinal7575_Nonzero_HiSeq_Broad_WGS,
  metaRSFinal7575_Nonzero_HiSeq_UNC,
  metaRSFinal7575_Nonzero_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaRSFinal7575_Nonzero_HiSeq_WGS,
  metaRSFinal7575_Nonzero_HiSeq_RNA,
  
  # Full metadata
  metaRSFinal7575_Nonzero_HiSeq,
  file = "Interim_data/RS210_Pan_Filt7575_RawVsCQ_data_for_ml_tcga_by_seq_center_15Oct23.RData")

#----------------------------------------------------------#
# Subset Filt9090 data
#----------------------------------------------------------#

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_HMS <- rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal9090_Nonzero_HiSeq_HMS),]
rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_BCM <- rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal9090_Nonzero_HiSeq_BCM),]
rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_MDA <- rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal9090_Nonzero_HiSeq_MDA),]
rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_WashU <- rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal9090_Nonzero_HiSeq_WashU),]
rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_Broad_WGS <- rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ[rownames(metaRSFinal9090_Nonzero_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both RNA and RNA-Seq, so separate objects are made for both)
rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ_UNC <- rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ[rownames(metaRSFinal9090_Nonzero_HiSeq_UNC),]
rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ_CMS <- rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ[rownames(metaRSFinal9090_Nonzero_HiSeq_CMS),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  rs210PanFinal9090_Nonzero_HiSeq_HMS,
  rs210PanFinal9090_Nonzero_HiSeq_BCM,
  rs210PanFinal9090_Nonzero_HiSeq_MDA,
  rs210PanFinal9090_Nonzero_HiSeq_WashU,
  rs210PanFinal9090_Nonzero_HiSeq_Broad_WGS,
  rs210PanFinal9090_Nonzero_HiSeq_UNC,
  rs210PanFinal9090_Nonzero_HiSeq_CMS,
  
  # ConQuR-corrected data
  rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_HMS,
  rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_BCM,
  rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_MDA,
  rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_WashU,
  rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_Broad_WGS,
  rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ_UNC,
  rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ_CMS,
  
  # WGS and RNA subset data
  rs210PanFinal9090_Nonzero_HiSeq_WGS,
  rs210PanFinal9090_Nonzero_HiSeq_RNA,
  rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ,
  rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ,
  
  # Subset metadata
  metaRSFinal9090_Nonzero_HiSeq_HMS,
  metaRSFinal9090_Nonzero_HiSeq_BCM,
  metaRSFinal9090_Nonzero_HiSeq_MDA,
  metaRSFinal9090_Nonzero_HiSeq_WashU,
  metaRSFinal9090_Nonzero_HiSeq_Broad_WGS,
  metaRSFinal9090_Nonzero_HiSeq_UNC,
  metaRSFinal9090_Nonzero_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaRSFinal9090_Nonzero_HiSeq_WGS,
  metaRSFinal9090_Nonzero_HiSeq_RNA,
  
  # Full metadata
  metaRSFinal9090_Nonzero_HiSeq,
  file = "Interim_data/RS210_Pan_Filt9090_RawVsCQ_data_for_ml_tcga_by_seq_center_15Oct23.RData")

