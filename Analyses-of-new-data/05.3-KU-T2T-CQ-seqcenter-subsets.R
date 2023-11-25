# 05.3-KU-T2T-CQ-seqcenter-subsets.R
# Author: Greg Poore
# Date: Oct 10, 2023
# Purposes:
# - Load KrakenUniq data and format

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
load("Interim_data/KU_T2T_Full_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("Interim_data/KU_T2T_WIS_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("Interim_data/KU_T2T_BIO_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("Interim_data/KU_T2T_Filt_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("Interim_data/KU_T2T_BIOFilt_data_for_ml_tcga_by_seq_center_14Oct23.RData")

# # Load ConQuR-corrected data
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_RNA_CQ_CMS.RData")

#----------------------------------------------------------#
# Subset WIS data
#----------------------------------------------------------#

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_HMS <- kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_HMS),]
kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_BCM <- kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_BCM),]
kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_MDA <- kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_MDA),]
kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_WashU <- kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_WashU),]
kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_Broad_WGS <- kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both RNA and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ_UNC <- kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_UNC),]
kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ_CMS <- kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_CMS),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_WIS_HiSeq_HMS,
  kuT2TFinalNonzero_WIS_HiSeq_BCM,
  kuT2TFinalNonzero_WIS_HiSeq_MDA,
  kuT2TFinalNonzero_WIS_HiSeq_WashU,
  kuT2TFinalNonzero_WIS_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_WIS_HiSeq_UNC,
  kuT2TFinalNonzero_WIS_HiSeq_CMS,
  
  # ConQuR-corrected data
  kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_HMS,
  kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_BCM,
  kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_MDA,
  kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_WashU,
  kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_Broad_WGS,
  kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ_UNC,
  kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_WIS_HiSeq_WGS,
  kuT2TFinalNonzero_WIS_HiSeq_RNA,
  kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ,
  kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ,
  
  # Subset metadata
  metaKUT2TFinalNonzero_WIS_HiSeq_HMS,
  metaKUT2TFinalNonzero_WIS_HiSeq_BCM,
  metaKUT2TFinalNonzero_WIS_HiSeq_MDA,
  metaKUT2TFinalNonzero_WIS_HiSeq_WashU,
  metaKUT2TFinalNonzero_WIS_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_WIS_HiSeq_UNC,
  metaKUT2TFinalNonzero_WIS_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_WIS_HiSeq_WGS,
  metaKUT2TFinalNonzero_WIS_HiSeq_RNA,
  
  # WIS metadata
  metaKUT2TFinalNonzero_WIS_HiSeq,
  file = "Interim_data/KU_T2T_WIS_RawVsCQ_data_for_ml_tcga_by_seq_center_15Oct23.RData")

#----------------------------------------------------------#
# Subset BIO data
#----------------------------------------------------------#

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_HMS <- kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_HMS),]
kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_BCM <- kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_BCM),]
kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_MDA <- kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_MDA),]
kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_WashU <- kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_WashU),]
kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_Broad_WGS <- kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both RNA and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ_UNC <- kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_UNC),]
kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ_CMS <- kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_CMS),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_BIO_HiSeq_HMS,
  kuT2TFinalNonzero_BIO_HiSeq_BCM,
  kuT2TFinalNonzero_BIO_HiSeq_MDA,
  kuT2TFinalNonzero_BIO_HiSeq_WashU,
  kuT2TFinalNonzero_BIO_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_BIO_HiSeq_UNC,
  kuT2TFinalNonzero_BIO_HiSeq_CMS,
  
  # ConQuR-corrected data
  kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_HMS,
  kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_BCM,
  kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_MDA,
  kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_WashU,
  kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_Broad_WGS,
  kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ_UNC,
  kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_BIO_HiSeq_WGS,
  kuT2TFinalNonzero_BIO_HiSeq_RNA,
  kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ,
  kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ,
  
  # Subset metadata
  metaKUT2TFinalNonzero_BIO_HiSeq_HMS,
  metaKUT2TFinalNonzero_BIO_HiSeq_BCM,
  metaKUT2TFinalNonzero_BIO_HiSeq_MDA,
  metaKUT2TFinalNonzero_BIO_HiSeq_WashU,
  metaKUT2TFinalNonzero_BIO_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_BIO_HiSeq_UNC,
  metaKUT2TFinalNonzero_BIO_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_BIO_HiSeq_WGS,
  metaKUT2TFinalNonzero_BIO_HiSeq_RNA,
  
  # BIO metadata
  metaKUT2TFinalNonzero_BIO_HiSeq,
  file = "Interim_data/KU_T2T_BIO_RawVsCQ_data_for_ml_tcga_by_seq_center_15Oct23.RData")

#----------------------------------------------------------#
# Subset Filt data
#----------------------------------------------------------#

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_HMS <- kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_HMS),]
kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_BCM <- kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_BCM),]
kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_MDA <- kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_MDA),]
kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_WashU <- kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_WashU),]
kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_Broad_WGS <- kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both RNA and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ_UNC <- kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_UNC),]
kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ_CMS <- kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_CMS),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_Filt_HiSeq_HMS,
  kuT2TFinalNonzero_Filt_HiSeq_BCM,
  kuT2TFinalNonzero_Filt_HiSeq_MDA,
  kuT2TFinalNonzero_Filt_HiSeq_WashU,
  kuT2TFinalNonzero_Filt_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_Filt_HiSeq_UNC,
  kuT2TFinalNonzero_Filt_HiSeq_CMS,
  
  # ConQuR-corrected data
  kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_HMS,
  kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_BCM,
  kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_MDA,
  kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_WashU,
  kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_Broad_WGS,
  kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ_UNC,
  kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_Filt_HiSeq_WGS,
  kuT2TFinalNonzero_Filt_HiSeq_RNA,
  kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ,
  kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ,
  
  # Subset metadata
  metaKUT2TFinalNonzero_Filt_HiSeq_HMS,
  metaKUT2TFinalNonzero_Filt_HiSeq_BCM,
  metaKUT2TFinalNonzero_Filt_HiSeq_MDA,
  metaKUT2TFinalNonzero_Filt_HiSeq_WashU,
  metaKUT2TFinalNonzero_Filt_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_Filt_HiSeq_UNC,
  metaKUT2TFinalNonzero_Filt_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_Filt_HiSeq_WGS,
  metaKUT2TFinalNonzero_Filt_HiSeq_RNA,
  
  # Filt metadata
  metaKUT2TFinalNonzero_Filt_HiSeq,
  file = "Interim_data/KU_T2T_Filt_RawVsCQ_data_for_ml_tcga_by_seq_center_15Oct23.RData")

#----------------------------------------------------------#
# Subset BIOFilt data
#----------------------------------------------------------#

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_HMS <- kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_HMS),]
kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_BCM <- kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_BCM),]
kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_MDA <- kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_MDA),]
kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_WashU <- kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_WashU),]
kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_Broad_WGS <- kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both RNA and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ_UNC <- kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_UNC),]
kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ_CMS <- kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_CMS),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_BIOFilt_HiSeq_HMS,
  kuT2TFinalNonzero_BIOFilt_HiSeq_BCM,
  kuT2TFinalNonzero_BIOFilt_HiSeq_MDA,
  kuT2TFinalNonzero_BIOFilt_HiSeq_WashU,
  kuT2TFinalNonzero_BIOFilt_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_BIOFilt_HiSeq_UNC,
  kuT2TFinalNonzero_BIOFilt_HiSeq_CMS,
  
  # ConQuR-corrected data
  kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_HMS,
  kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_BCM,
  kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_MDA,
  kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_WashU,
  kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_Broad_WGS,
  kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ_UNC,
  kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_BIOFilt_HiSeq_WGS,
  kuT2TFinalNonzero_BIOFilt_HiSeq_RNA,
  kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ,
  kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ,
  
  # Subset metadata
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_HMS,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_BCM,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_MDA,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_WashU,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_UNC,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_WGS,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_RNA,
  
  # BIOFilt metadata
  metaKUT2TFinalNonzero_BIOFilt_HiSeq,
  file = "Interim_data/KU_T2T_BIOFilt_RawVsCQ_data_for_ml_tcga_by_seq_center_15Oct23.RData")

#----------------------------------------------------------#
# Subset Full data
#----------------------------------------------------------#

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_HMS <- kuT2TFinalNonzero_Full_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Full_HiSeq_HMS),]
kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_BCM <- kuT2TFinalNonzero_Full_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Full_HiSeq_BCM),]
kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_MDA <- kuT2TFinalNonzero_Full_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Full_HiSeq_MDA),]
kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_WashU <- kuT2TFinalNonzero_Full_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Full_HiSeq_WashU),]
kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_Broad_WGS <- kuT2TFinalNonzero_Full_HiSeq_WGS_CQ[rownames(metaKUT2TFinalNonzero_Full_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both RNA and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_Full_HiSeq_RNA_CQ_UNC <- kuT2TFinalNonzero_Full_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_Full_HiSeq_UNC),]
kuT2TFinalNonzero_Full_HiSeq_RNA_CQ_CMS <- kuT2TFinalNonzero_Full_HiSeq_RNA_CQ[rownames(metaKUT2TFinalNonzero_Full_HiSeq_CMS),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_Full_HiSeq_HMS,
  kuT2TFinalNonzero_Full_HiSeq_BCM,
  kuT2TFinalNonzero_Full_HiSeq_MDA,
  kuT2TFinalNonzero_Full_HiSeq_WashU,
  kuT2TFinalNonzero_Full_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_Full_HiSeq_UNC,
  kuT2TFinalNonzero_Full_HiSeq_CMS,
  
  # ConQuR-corrected data
  kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_HMS,
  kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_BCM,
  kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_MDA,
  kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_WashU,
  kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_Broad_WGS,
  kuT2TFinalNonzero_Full_HiSeq_RNA_CQ_UNC,
  kuT2TFinalNonzero_Full_HiSeq_RNA_CQ_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_Full_HiSeq_WGS,
  kuT2TFinalNonzero_Full_HiSeq_RNA,
  kuT2TFinalNonzero_Full_HiSeq_WGS_CQ,
  kuT2TFinalNonzero_Full_HiSeq_RNA_CQ,
  
  # Subset metadata
  metaKUT2TFinalNonzero_Full_HiSeq_HMS,
  metaKUT2TFinalNonzero_Full_HiSeq_BCM,
  metaKUT2TFinalNonzero_Full_HiSeq_MDA,
  metaKUT2TFinalNonzero_Full_HiSeq_WashU,
  metaKUT2TFinalNonzero_Full_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_Full_HiSeq_UNC,
  metaKUT2TFinalNonzero_Full_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_Full_HiSeq_WGS,
  metaKUT2TFinalNonzero_Full_HiSeq_RNA,
  
  # Full metadata
  metaKUT2TFinalNonzero_Full_HiSeq,
  file = "Interim_data/KU_T2T_Full_RawVsCQ_data_for_ml_tcga_by_seq_center_15Oct23.RData")
