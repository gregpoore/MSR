#-----------------------------------------------------------------------------
# conqur-rs210-wis-wgs.R
# Copyright (c) 2023--, Greg Poore
# Purpose: ConQuR batch correction by seq center
#-----------------------------------------------------------------------------

#-------------------------------#
# Load dependencies
library(ConQuR)
library(doParallel)

## Import data
load("../../Interim_data/RS210_WIS_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("../../Interim_data/RS210_Filt9090_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("../../Interim_data/RS210_Filt7575_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("../../Interim_data/RS210_Filt5050_data_for_ml_tcga_by_seq_center_13Oct23.RData")

#-----------------------------------Filt5050-----------------------------------#

## WGS
batchid_TCGA_WGS <- factor(metaRSFinal5050_Nonzero_HiSeq_WGS[,"data_submitting_center_label"])
covar_TCGA_WGS <- metaRSFinal5050_Nonzero_HiSeq_WGS[,"sample_type",drop = FALSE]
covar_TCGA_WGS$sample_type <- factor(covar_TCGA_WGS$sample_type)
rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ <- as.data.frame(ConQuR(tax_tab = rs210PanFinal5050_Nonzero_HiSeq_WGS,
                                                               batchid = batchid_TCGA_WGS,
                                                               covariates = covar_TCGA_WGS,
                                                               num_core = 32,
                                                               batch_ref = "Baylor College of Medicine"))
save(rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ,
     file = "../../Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_BCM.RData")

## RNA
batchid_TCGA_RNA <- factor(metaRSFinal5050_Nonzero_HiSeq_RNA[,"data_submitting_center_label"])
covar_TCGA_RNA <- metaRSFinal5050_Nonzero_HiSeq_RNA[,"sample_type",drop = FALSE]
covar_TCGA_RNA$sample_type <- factor(covar_TCGA_RNA$sample_type)
rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ <- as.data.frame(ConQuR(tax_tab = rs210PanFinal5050_Nonzero_HiSeq_RNA,
                                                               batchid = batchid_TCGA_RNA,
                                                               covariates = covar_TCGA_RNA,
                                                               num_core = 32,
                                                               batch_ref = "Canada's Michael Smith Genome Sciences Centre"))
save(rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ,
     file = "../../Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ_CMS.RData")

#-----------------------------------WIS-----------------------------------#
## WGS
batchid_TCGA_WGS <- factor(metaRSFinalWISNonzero_HiSeq_WGS[,"data_submitting_center_label"])
covar_TCGA_WGS <- metaRSFinalWISNonzero_HiSeq_WGS[,"sample_type",drop = FALSE]
covar_TCGA_WGS$sample_type <- factor(covar_TCGA_WGS$sample_type)
rs210PanFinalWISNonzero_HiSeq_WGS_CQ <- as.data.frame(ConQuR(tax_tab = rs210PanFinalWISNonzero_HiSeq_WGS,
                                                 batchid = batchid_TCGA_WGS,
                                                 covariates = covar_TCGA_WGS,
                                                 num_core = 32,
                                                 batch_ref = "Baylor College of Medicine"))
save(rs210PanFinalWISNonzero_HiSeq_WGS_CQ,
     file = "../../Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_WGS_CQ_BCM.RData")

## RNA
batchid_TCGA_RNA <- factor(metaRSFinalWISNonzero_HiSeq_RNA[,"data_submitting_center_label"])
covar_TCGA_RNA <- metaRSFinalWISNonzero_HiSeq_RNA[,"sample_type",drop = FALSE]
covar_TCGA_RNA$sample_type <- factor(covar_TCGA_RNA$sample_type)
rs210PanFinalWISNonzero_HiSeq_RNA_CQ <- as.data.frame(ConQuR(tax_tab = rs210PanFinalWISNonzero_HiSeq_RNA,
                                                             batchid = batchid_TCGA_RNA,
                                                             covariates = covar_TCGA_RNA,
                                                             num_core = 32,
                                                             batch_ref = "Canada's Michael Smith Genome Sciences Centre"))
save(rs210PanFinalWISNonzero_HiSeq_RNA_CQ,
     file = "../../Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_RNA_CQ_CMS.RData")

#-----------------------------------Filt7575-----------------------------------#

## WGS
batchid_TCGA_WGS <- factor(metaRSFinal7575_Nonzero_HiSeq_WGS[,"data_submitting_center_label"])
covar_TCGA_WGS <- metaRSFinal7575_Nonzero_HiSeq_WGS[,"sample_type",drop = FALSE]
covar_TCGA_WGS$sample_type <- factor(covar_TCGA_WGS$sample_type)
rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ <- as.data.frame(ConQuR(tax_tab = rs210PanFinal7575_Nonzero_HiSeq_WGS,
                                                               batchid = batchid_TCGA_WGS,
                                                               covariates = covar_TCGA_WGS,
                                                               num_core = 32,
                                                               batch_ref = "Baylor College of Medicine"))
save(rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ,
     file = "../../Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_BCM.RData")

## RNA
batchid_TCGA_RNA <- factor(metaRSFinal7575_Nonzero_HiSeq_RNA[,"data_submitting_center_label"])
covar_TCGA_RNA <- metaRSFinal7575_Nonzero_HiSeq_RNA[,"sample_type",drop = FALSE]
covar_TCGA_RNA$sample_type <- factor(covar_TCGA_RNA$sample_type)
rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ <- as.data.frame(ConQuR(tax_tab = rs210PanFinal7575_Nonzero_HiSeq_RNA,
                                                               batchid = batchid_TCGA_RNA,
                                                               covariates = covar_TCGA_RNA,
                                                               num_core = 32,
                                                               batch_ref = "Canada's Michael Smith Genome Sciences Centre"))
save(rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ,
     file = "../../Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ_CMS.RData")

#-----------------------------------Filt9090-----------------------------------#

## WGS
batchid_TCGA_WGS <- factor(metaRSFinal9090_Nonzero_HiSeq_WGS[,"data_submitting_center_label"])
covar_TCGA_WGS <- metaRSFinal9090_Nonzero_HiSeq_WGS[,"sample_type",drop = FALSE]
covar_TCGA_WGS$sample_type <- factor(covar_TCGA_WGS$sample_type)
rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ <- as.data.frame(ConQuR(tax_tab = rs210PanFinal9090_Nonzero_HiSeq_WGS,
                                                               batchid = batchid_TCGA_WGS,
                                                               covariates = covar_TCGA_WGS,
                                                               num_core = 32,
                                                               batch_ref = "Baylor College of Medicine"))
save(rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ,
     file = "../../Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_BCM.RData")

## RNA
batchid_TCGA_RNA <- factor(metaRSFinal9090_Nonzero_HiSeq_RNA[,"data_submitting_center_label"])
covar_TCGA_RNA <- metaRSFinal9090_Nonzero_HiSeq_RNA[,"sample_type",drop = FALSE]
covar_TCGA_RNA$sample_type <- factor(covar_TCGA_RNA$sample_type)
rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ <- as.data.frame(ConQuR(tax_tab = rs210PanFinal9090_Nonzero_HiSeq_RNA,
                                                               batchid = batchid_TCGA_RNA,
                                                               covariates = covar_TCGA_RNA,
                                                               num_core = 32,
                                                               batch_ref = "Canada's Michael Smith Genome Sciences Centre"))
save(rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ,
     file = "../../Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ_CMS.RData")



