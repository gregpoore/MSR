#-----------------------------------------------------------------------------
# conqur-ku-t2t-wis-wgs-rna.R
# Copyright (c) 2023--, Greg Poore
# Purpose: ConQuR batch correction by seq center
#-----------------------------------------------------------------------------

#-------------------------------#
# Load dependencies
library(ConQuR)
library(doParallel)

## Import data
load("../../Interim_data/KU_T2T_Full_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../../Interim_data/KU_T2T_WIS_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../../Interim_data/KU_T2T_BIO_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../../Interim_data/KU_T2T_Filt_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../../Interim_data/KU_T2T_BIOFilt_data_for_ml_tcga_by_seq_center_14Oct23.RData")

#-----------------------------------BIO-----------------------------------#
## WGS
batchid_TCGA_WGS <- factor(metaKUT2TFinalNonzero_BIO_HiSeq_WGS[,"data_submitting_center_label"])
covar_TCGA_WGS <- metaKUT2TFinalNonzero_BIO_HiSeq_WGS[,"sample_type",drop = FALSE]
covar_TCGA_WGS$sample_type <- factor(covar_TCGA_WGS$sample_type)
kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_BIO_HiSeq_WGS,
                                                           batchid = batchid_TCGA_WGS,
                                                           covariates = covar_TCGA_WGS,
                                                           num_core = 32,
                                                           batch_ref = "Baylor College of Medicine"))
save(kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_BCM.RData")

## RNA
batchid_TCGA_RNA <- factor(metaKUT2TFinalNonzero_BIO_HiSeq_RNA[,"data_submitting_center_label"])
covar_TCGA_RNA <- metaKUT2TFinalNonzero_BIO_HiSeq_RNA[,"sample_type",drop = FALSE]
covar_TCGA_RNA$sample_type <- factor(covar_TCGA_RNA$sample_type)
kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_BIO_HiSeq_RNA,
                                                           batchid = batchid_TCGA_RNA,
                                                           covariates = covar_TCGA_RNA,
                                                           num_core = 32,
                                                           batch_ref = "Canada's Michael Smith Genome Sciences Centre"))
save(kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ_CMS.RData")

#-----------------------------------WIS-----------------------------------#
## WGS
batchid_TCGA_WGS <- factor(metaKUT2TFinalNonzero_WIS_HiSeq_WGS[,"data_submitting_center_label"])
covar_TCGA_WGS <- metaKUT2TFinalNonzero_WIS_HiSeq_WGS[,"sample_type",drop = FALSE]
covar_TCGA_WGS$sample_type <- factor(covar_TCGA_WGS$sample_type)
kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_WIS_HiSeq_WGS,
                                                           batchid = batchid_TCGA_WGS,
                                                           covariates = covar_TCGA_WGS,
                                                           num_core = 32,
                                                           batch_ref = "Baylor College of Medicine"))
save(kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_BCM.RData")

## RNA
batchid_TCGA_RNA <- factor(metaKUT2TFinalNonzero_WIS_HiSeq_RNA[,"data_submitting_center_label"])
covar_TCGA_RNA <- metaKUT2TFinalNonzero_WIS_HiSeq_RNA[,"sample_type",drop = FALSE]
covar_TCGA_RNA$sample_type <- factor(covar_TCGA_RNA$sample_type)
kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_WIS_HiSeq_RNA,
                                                           batchid = batchid_TCGA_RNA,
                                                           covariates = covar_TCGA_RNA,
                                                           num_core = 32,
                                                           batch_ref = "Canada's Michael Smith Genome Sciences Centre"))
save(kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ_CMS.RData")

#-----------------------------------Filt-----------------------------------#
## WGS
batchid_TCGA_WGS <- factor(metaKUT2TFinalNonzero_Filt_HiSeq_WGS[,"data_submitting_center_label"])
covar_TCGA_WGS <- metaKUT2TFinalNonzero_Filt_HiSeq_WGS[,"sample_type",drop = FALSE]
covar_TCGA_WGS$sample_type <- factor(covar_TCGA_WGS$sample_type)
kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_Filt_HiSeq_WGS,
                                                           batchid = batchid_TCGA_WGS,
                                                           covariates = covar_TCGA_WGS,
                                                           num_core = 32,
                                                           batch_ref = "Baylor College of Medicine"))
save(kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_BCM.RData")

## RNA
batchid_TCGA_RNA <- factor(metaKUT2TFinalNonzero_Filt_HiSeq_RNA[,"data_submitting_center_label"])
covar_TCGA_RNA <- metaKUT2TFinalNonzero_Filt_HiSeq_RNA[,"sample_type",drop = FALSE]
covar_TCGA_RNA$sample_type <- factor(covar_TCGA_RNA$sample_type)
kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_Filt_HiSeq_RNA,
                                                           batchid = batchid_TCGA_RNA,
                                                           covariates = covar_TCGA_RNA,
                                                           num_core = 32,
                                                           batch_ref = "Canada's Michael Smith Genome Sciences Centre"))
save(kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ_CMS.RData")

#-----------------------------------BIOFilt-----------------------------------#
## WGS
batchid_TCGA_WGS <- factor(metaKUT2TFinalNonzero_BIOFilt_HiSeq_WGS[,"data_submitting_center_label"])
covar_TCGA_WGS <- metaKUT2TFinalNonzero_BIOFilt_HiSeq_WGS[,"sample_type",drop = FALSE]
covar_TCGA_WGS$sample_type <- factor(covar_TCGA_WGS$sample_type)
kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_BIOFilt_HiSeq_WGS,
                                                            batchid = batchid_TCGA_WGS,
                                                            covariates = covar_TCGA_WGS,
                                                            num_core = 32,
                                                            batch_ref = "Baylor College of Medicine"))
save(kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_BCM.RData")

## RNA
batchid_TCGA_RNA <- factor(metaKUT2TFinalNonzero_BIOFilt_HiSeq_RNA[,"data_submitting_center_label"])
covar_TCGA_RNA <- metaKUT2TFinalNonzero_BIOFilt_HiSeq_RNA[,"sample_type",drop = FALSE]
covar_TCGA_RNA$sample_type <- factor(covar_TCGA_RNA$sample_type)
kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_BIOFilt_HiSeq_RNA,
                                                            batchid = batchid_TCGA_RNA,
                                                            covariates = covar_TCGA_RNA,
                                                            num_core = 32,
                                                            batch_ref = "Canada's Michael Smith Genome Sciences Centre"))
save(kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ_CMS.RData")

#-----------------------------------Full-----------------------------------#
## WGS
batchid_TCGA_WGS <- factor(metaKUT2TFinalNonzero_Full_HiSeq_WGS[,"data_submitting_center_label"])
covar_TCGA_WGS <- metaKUT2TFinalNonzero_Full_HiSeq_WGS[,"sample_type",drop = FALSE]
covar_TCGA_WGS$sample_type <- factor(covar_TCGA_WGS$sample_type)
kuT2TFinalNonzero_Full_HiSeq_WGS_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_Full_HiSeq_WGS,
                                                           batchid = batchid_TCGA_WGS,
                                                           covariates = covar_TCGA_WGS,
                                                           num_core = 32,
                                                           batch_ref = "Baylor College of Medicine"))
save(kuT2TFinalNonzero_Full_HiSeq_WGS_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_BCM.RData")

## RNA
batchid_TCGA_RNA <- factor(metaKUT2TFinalNonzero_Full_HiSeq_RNA[,"data_submitting_center_label"])
covar_TCGA_RNA <- metaKUT2TFinalNonzero_Full_HiSeq_RNA[,"sample_type",drop = FALSE]
covar_TCGA_RNA$sample_type <- factor(covar_TCGA_RNA$sample_type)
kuT2TFinalNonzero_Full_HiSeq_RNA_CQ <- as.data.frame(ConQuR(tax_tab = kuT2TFinalNonzero_Full_HiSeq_RNA,
                                                           batchid = batchid_TCGA_RNA,
                                                           covariates = covar_TCGA_RNA,
                                                           num_core = 32,
                                                           batch_ref = "Canada's Michael Smith Genome Sciences Centre"))
save(kuT2TFinalNonzero_Full_HiSeq_RNA_CQ,
     file = "../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_RNA_CQ_CMS.RData")
