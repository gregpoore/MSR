# 02.3-wis-raw-vsnm-cq-data-subsetting.R
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
load("Interim_data/tcga_WIS_hiseq_VSNM_CBS_CQ_31Aug23.RData", verbose = TRUE)
# metadataSamplesAllQC_HiSeq_WGS
# metadataSamplesAllQC_HiSeq_RNA
# tcgaGenusKrakenAllFiltWIS_HiSeq_WGS
# tcgaGenusKrakenAllFiltWIS_HiSeq_RNA
# tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM
# tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM
# tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CBS
# tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CBS
# tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ
# tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ

##--------------------Subset count data by seq center--------------------##
#--------------------Raw--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_HMS <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS[rownames(metadataSamplesAllQC_HiSeq_HMS),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_BCM <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS[rownames(metadataSamplesAllQC_HiSeq_BCM),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_MDA <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS[rownames(metadataSamplesAllQC_HiSeq_MDA),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_WashU <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS[rownames(metadataSamplesAllQC_HiSeq_WashU),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_Broad_WGS <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_UNC <- tcgaGenusKrakenAllFiltWIS_HiSeq_RNA[rownames(metadataSamplesAllQC_HiSeq_UNC),]
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CMS <- tcgaGenusKrakenAllFiltWIS_HiSeq_RNA[rownames(metadataSamplesAllQC_HiSeq_CMS),]
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_Broad_RNA <- tcgaGenusKrakenAllFiltWIS_HiSeq_RNA[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------VSNM--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_HMS <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM[rownames(metadataSamplesAllQC_HiSeq_HMS),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_BCM <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM[rownames(metadataSamplesAllQC_HiSeq_BCM),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_MDA <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM[rownames(metadataSamplesAllQC_HiSeq_MDA),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_WashU <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM[rownames(metadataSamplesAllQC_HiSeq_WashU),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_Broad_WGS <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM_UNC <- tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM[rownames(metadataSamplesAllQC_HiSeq_UNC),]
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM_CMS <- tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM[rownames(metadataSamplesAllQC_HiSeq_CMS),]
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM_Broad_RNA <- tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------CQ--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_HMS <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ[rownames(metadataSamplesAllQC_HiSeq_HMS),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_BCM <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ[rownames(metadataSamplesAllQC_HiSeq_BCM),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_MDA <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ[rownames(metadataSamplesAllQC_HiSeq_MDA),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_WashU <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ[rownames(metadataSamplesAllQC_HiSeq_WashU),]
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_Broad_WGS <- tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_UNC <- tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ[rownames(metadataSamplesAllQC_HiSeq_UNC),]
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_CMS <- tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ[rownames(metadataSamplesAllQC_HiSeq_CMS),]
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_Broad_RNA <- tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

## Save data
save(# Subset raw count data
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_HMS,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_BCM,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_MDA,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_WashU,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_Broad_WGS,
  tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_UNC,
  tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CMS,
  tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_Broad_RNA,
  
  # VSNM
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_HMS,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_BCM,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_MDA,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_WashU,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM_Broad_WGS,
  tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM_UNC,
  tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM_CMS,
  tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM_Broad_RNA,
  
  # CQ
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_HMS,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_BCM,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_MDA,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_WashU,
  tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_Broad_WGS,
  tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_UNC,
  tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_CMS,
  tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_Broad_RNA,
  
  # Subset metadata
  metadataSamplesAllQC_HiSeq_HMS,
  metadataSamplesAllQC_HiSeq_BCM,
  metadataSamplesAllQC_HiSeq_MDA,
  metadataSamplesAllQC_HiSeq_WashU,
  metadataSamplesAllQC_HiSeq_Broad_WGS,
  metadataSamplesAllQC_HiSeq_UNC,
  metadataSamplesAllQC_HiSeq_CMS,
  metadataSamplesAllQC_HiSeq_Broad_RNA,
  file = "Interim_data/raw_vsnm_cq_data_for_ml_tcga_by_seq_center_and_experimental_strategy_1Sept23.RData")



