# 02.2-wis-batch-correction-comparisons.R
# Author: Greg Poore
# Date: Aug 18, 2023
# Purpose: Apply multiple batch correction options

#-------------------------------#
# Load dependencies
# require(doMC)
require(plyr)
require(dplyr)
require(tibble)
require(reshape2)
require(phyloseq)

# numCores <- detectCores()
# registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import data
#----------------------------------------------------------#

# Full data
load("Input_data/tcgaVbDataAndMetadataAndSNM_consolidated_Nov23.RData", verbose = TRUE)
# metadataSamplesAll
# metadataSamplesAllQC
# metadataSamplesAllQCSurvival
# vbDataBarnDFReconciled
# vbDataBarnDFReconciledQC
# vbDataBarnDFReconciledQCSurvival
# snmDataSampleType

# WIS data
load("Interim_data/tcga-hiseq-wis-overlapping-data-and-metadata-subset-25July22.RData", verbose = T)
# tcgaGenusKrakenAllFiltWIS_HiSeq
# tcgaGenusKrakenQCFiltWIS_HiSeq
# metadataSamplesAll_HiSeq
# metadataSamplesAllQC_HiSeq
#----------------------------------------------------------#
# Subset metadata and count data
#----------------------------------------------------------#

## Subset metadata
metadataSamplesAllQC_HiSeq_WGS <- metadataSamplesAllQC_HiSeq %>%
  filter(sample_type %in% c("Primary Tumor","Blood Derived Normal","Solid Tissue Normal")) %>%
  filter(experimental_strategy == "WGS") %>% droplevels()

metadataSamplesAllQC_HiSeq_RNA <- metadataSamplesAllQC_HiSeq %>%
  filter(sample_type %in% c("Primary Tumor","Blood Derived Normal","Solid Tissue Normal")) %>%
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

## Examine sequencing center sample distributions
metadataSamplesAllQC_HiSeq_WGS %>%
  count(data_submitting_center_label)

# data_submitting_center_label    n
# 1                         Baylor College of Medicine  641
# 2                 Broad Institute of MIT and Harvard  810
# 3                             Harvard Medical School 2007
# 4 MD Anderson - Institute for Applied Cancer Science  526
# 5      MD Anderson - RPPA Core Facility (Proteomics)   11
# 6           Washington University School of Medicine  414

metadataSamplesAllQC_HiSeq_RNA %>%
  count(data_submitting_center_label)

# data_submitting_center_label    n
# 1            Broad Institute of MIT and Harvard  337
# 2 Canada's Michael Smith Genome Sciences Centre 1732
# 3                  University of North Carolina 9609

## Subset full count data
vbDataBarnDFReconciledQC_HiSeq_WGS <- vbDataBarnDFReconciledQC[rownames(metadataSamplesAllQC_HiSeq_WGS),]
vbDataBarnDFReconciledQC_HiSeq_RNA <- vbDataBarnDFReconciledQC[rownames(metadataSamplesAllQC_HiSeq_RNA),]

dim(vbDataBarnDFReconciledQC_HiSeq_WGS) # 4231 1993
dim(vbDataBarnDFReconciledQC_HiSeq_RNA) # 11202  1993

# Sanity checks
sum(rowSums(vbDataBarnDFReconciledQC_HiSeq_WGS)==0) # 0
sum(rowSums(vbDataBarnDFReconciledQC_HiSeq_RNA)==0) # 0
sum(colSums(vbDataBarnDFReconciledQC_HiSeq_WGS)==0) # 44
sum(colSums(vbDataBarnDFReconciledQC_HiSeq_RNA)==0) # 58

# Remove zero-sum features
vbDataBarnDFReconciledQC_Nonzero_HiSeq_WGS <- vbDataBarnDFReconciledQC_HiSeq_WGS[,(colSums(vbDataBarnDFReconciledQC_HiSeq_WGS)!=0)]
vbDataBarnDFReconciledQC_Nonzero_HiSeq_RNA <- vbDataBarnDFReconciledQC_HiSeq_RNA[,(colSums(vbDataBarnDFReconciledQC_HiSeq_RNA)!=0)]

dim(vbDataBarnDFReconciledQC_Nonzero_HiSeq_WGS) # 4231 1949
dim(vbDataBarnDFReconciledQC_Nonzero_HiSeq_RNA) # 11202  1935

## Subset WIS count data
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS <- tcgaGenusKrakenAllFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WGS),]
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA <- tcgaGenusKrakenAllFiltWIS_HiSeq[rownames(metadataSamplesAllQC_HiSeq_RNA),]

dim(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS) # 4231  184
dim(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA) # 11202   184

# Sanity checks
sum(rowSums(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS)==0) # 0
sum(rowSums(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA)==0) # 0
sum(colSums(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS)==0) # 0
sum(colSums(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA)==0) # 0

#----------------------------------------------------------#
#----------------------------------------------------------#
# WIS data batch corrections
#----------------------------------------------------------#
#----------------------------------------------------------#

#----------------------------------------------------------#
# VSNM WIS-overlapping data
#----------------------------------------------------------#
require(limma)
require(edgeR)
require(snm)

##------------------WGS------------------##
bioFactorsVec <- c("sample_type")
techFactorsVec <- c("data_submitting_center_label")

# Set up design matrix
formulaString <- paste("~0",paste(bioFactorsVec, collapse = " + "),
                       paste(techFactorsVec, collapse = " + "), sep = " + ")
covDesignNorm <- model.matrix(stats::as.formula(formulaString),
                              data = metadataSamplesAllQC_HiSeq_WGS)
colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))

# Set up counts matrix
countsWGS <- t(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS) # DGEList object from a table of counts (rows=features, columns=samples)

# Normalize using edgeR and then plug into voom
print("Running Voom...")
dge <- DGEList(counts = countsWGS)
vdge_data <- voom(dge, design = covDesignNorm, plot = FALSE, save.plot = FALSE,
                  normalize.method="quantile")

# Define biological factors
bioVarFormula <- paste("~",paste(bioFactorsVec, collapse = " + "), sep = "")
bio.var <- model.matrix(stats::as.formula(bioVarFormula),
                        data=metadataSamplesAllQC_HiSeq_WGS)
# Define technical factors
techVarFormula <- paste("~",paste(techFactorsVec, collapse = " + "), sep = "")
adj.var <- model.matrix(stats::as.formula(techVarFormula),
                        data=metadataSamplesAllQC_HiSeq_WGS)
# Reformat column names
colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))

# Calculate SNM
print("Running SNM...")
snmDataObjOnlyWGS <- snm(raw.dat = vdge_data$E, 
                      bio.var = bio.var, 
                      adj.var = adj.var, 
                      rm.adj=TRUE,
                      verbose = TRUE,
                      diagnose = TRUE)
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM <- t(snmDataObjOnlyWGS$norm.dat)

##------------------RNA------------------##
bioFactorsVec <- c("sample_type")
techFactorsVec <- c("data_submitting_center_label")

# Set up design matrix
formulaString <- paste("~0",paste(bioFactorsVec, collapse = " + "),
                       paste(techFactorsVec, collapse = " + "), sep = " + ")
covDesignNorm <- model.matrix(stats::as.formula(formulaString),
                              data = metadataSamplesAllQC_HiSeq_RNA)
colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))

# Set up counts matrix
countsRNA <- t(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA) # DGEList object from a table of counts (rows=features, columns=samples)

# Normalize using edgeR and then plug into voom
print("Running Voom...")
dge <- DGEList(counts = countsRNA)
vdge_data <- voom(dge, design = covDesignNorm, plot = FALSE, save.plot = FALSE,
                  normalize.method="quantile")

# Define biological factors
bioVarFormula <- paste("~",paste(bioFactorsVec, collapse = " + "), sep = "")
bio.var <- model.matrix(stats::as.formula(bioVarFormula),
                        data=metadataSamplesAllQC_HiSeq_RNA)
# Define technical factors
techVarFormula <- paste("~",paste(techFactorsVec, collapse = " + "), sep = "")
adj.var <- model.matrix(stats::as.formula(techVarFormula),
                        data=metadataSamplesAllQC_HiSeq_RNA)
# Reformat column names
colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))

# Calculate SNM
print("Running SNM...")
snmDataObjOnlyRNA <- snm(raw.dat = vdge_data$E, 
                         bio.var = bio.var, 
                         adj.var = adj.var, 
                         rm.adj=TRUE,
                         verbose = TRUE,
                         diagnose = TRUE)
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM <- t(snmDataObjOnlyRNA$norm.dat)

#----------------------------------------------------------#
# Combat-seq WIS-overlapping data -- not used in final analysis
#----------------------------------------------------------#
require(sva)

tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CBS <- data.frame(t(ComBat_seq(t(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS), 
                                                        batch=metadataSamplesAllQC_HiSeq_WGS[,"data_submitting_center_label"], 
                                                        group=metadataSamplesAllQC_HiSeq_WGS[,"sample_type"],
                                                        full_mod = TRUE)))

tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CBS <- data.frame(t(ComBat_seq(t(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA), 
                                                                   batch=metadataSamplesAllQC_HiSeq_RNA[,"data_submitting_center_label"], 
                                                                   group=metadataSamplesAllQC_HiSeq_RNA[,"sample_type"],
                                                                   full_mod = TRUE)))

#----------------------------------------------------------#
# ConQur WIS-overlapping data
#----------------------------------------------------------#
library(ConQuR)
library(doParallel) 

##----------------------------WGS----------------------------##
batchid_TCGA_WGS <- metadataSamplesAllQC_HiSeq_WGS[,"data_submitting_center_label"]
covar_TCGA_WGS <- metadataSamplesAllQC_HiSeq_WGS[,"sample_type",drop = FALSE]
# The following line takes ~25 min to run
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ <- ConQuR(tax_tab = tcgaGenusKrakenAllFiltWIS_HiSeq_WGS,
                                                batchid=batchid_TCGA_WGS,
                                                covariates=covar_TCGA_WGS,
                                                batch_ref = "Harvard Medical School")
save(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ,
     file = "Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_19Aug23.RData")
# The following line takes >30 min to run
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP <- ConQuR(tax_tab = tcgaGenusKrakenAllFiltWIS_HiSeq_WGS,
                                                 batchid=batchid_TCGA_WGS,
                                                 covariates=covar_TCGA_WGS,
                                                 batch_ref = "Harvard Medical School",
                                                 logistic_lasso=T, 
                                                 quantile_type="lasso", 
                                                 interplt=T)
save(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP,
     file = "Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP_19Aug23.RData")

# The following line takes ~25 min to run
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_libsize <- ConQuR_libsize(tax_tab = tcgaGenusKrakenAllFiltWIS_HiSeq_WGS,
                                                 batchid=batchid_TCGA_WGS,
                                                 covariates=covar_TCGA_WGS,
                                                 batch_ref = "Harvard Medical School")
save(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_libsize,
     file = "Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_libsize_19Aug23.RData")
# The following line takes >30 min to run
tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP_libsize <- ConQuR_libsize(tax_tab = tcgaGenusKrakenAllFiltWIS_HiSeq_WGS,
                                                  batchid=batchid_TCGA_WGS,
                                                  covariates=covar_TCGA_WGS,
                                                  batch_ref = "Harvard Medical School",
                                                  logistic_lasso=T, 
                                                  quantile_type="lasso", 
                                                  interplt=T)
save(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP_libsize,
     file = "Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP_libsize_19Aug23.RData")

##----------------------------RNA----------------------------##

batchid_TCGA_RNA <- metadataSamplesAllQC_HiSeq_RNA[,"data_submitting_center_label"]
covar_TCGA_RNA <- metadataSamplesAllQC_HiSeq_RNA[,"sample_type",drop = FALSE]
# The following line takes an estimated 2 hr to run
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ <- ConQuR(tax_tab = tcgaGenusKrakenAllFiltWIS_HiSeq_RNA,
                                                 batchid=batchid_TCGA_RNA,
                                                 covariates=covar_TCGA_RNA,
                                                 batch_ref = "University of North Carolina")
save(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ,
     file = "Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_19Aug23.RData")

# The following line takes an estimated 2 hr to run
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP <- ConQuR(tax_tab = tcgaGenusKrakenAllFiltWIS_HiSeq_RNA,
                                                  batchid=batchid_TCGA_RNA,
                                                  covariates=covar_TCGA_RNA,
                                                  batch_ref = "University of North Carolina",
                                                  logistic_lasso=T, 
                                                  quantile_type="lasso", 
                                                  interplt=T)
save(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP,
     file = "Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP_19Aug23.RData")

# The following line takes an estimated 2 hr to run
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_libsize <- ConQuR_libsize(tax_tab = tcgaGenusKrakenAllFiltWIS_HiSeq_RNA,
                                                                 batchid=batchid_TCGA_RNA,
                                                                 covariates=covar_TCGA_RNA,
                                                                 batch_ref = "University of North Carolina")
save(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_libsize,
     file = "Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_libsize_19Aug23.RData")
# The following line takes an estimated 2 hr to run
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP_libsize <- ConQuR_libsize(tax_tab = tcgaGenusKrakenAllFiltWIS_HiSeq_RNA,
                                                                  batchid=batchid_TCGA_RNA,
                                                                  covariates=covar_TCGA_RNA,
                                                                  batch_ref = "University of North Carolina",
                                                                  logistic_lasso=T, 
                                                                  quantile_type="lasso", 
                                                                  interplt=T)
save(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP_libsize,
     file = "Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP_libsize_19Aug23.RData")

#-------------------------------------------#
##-------> Save using CQ only data (based on PVCA results below)

load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_19Aug23.RData", verbose = TRUE)
load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_19Aug23.RData", verbose = TRUE)

tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ <- as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ)
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ <- as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ)

tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM <- as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM)
tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM <- as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM)

save(metadataSamplesAllQC_HiSeq_WGS,
     metadataSamplesAllQC_HiSeq_RNA,
     tcgaGenusKrakenAllFiltWIS_HiSeq_WGS,
     tcgaGenusKrakenAllFiltWIS_HiSeq_RNA,
     tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM,
     tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM,
     tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CBS,
     tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CBS,
     tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ,
     tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ,
     file = "Interim_data/tcga_WIS_hiseq_VSNM_CBS_CQ_31Aug23.RData")

#----------------------------------------------------------#
# PVCA WIS-overlapping data
# NOTE: These analyses explore the batch correcting ability of the 
# above approaches -- the final plot is below this section
#----------------------------------------------------------#
require(doMC) # for parallel computing
require(lme4)

numCores <- detectCores()
registerDoMC(cores=numCores)

## Load raw, VSNM, and CBS datasets
load("Interim_data/tcga_WIS_hiseq_VSNM_CBS_18Aug23.RData")
# Load CQ data
load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_19Aug23.RData")
load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP_19Aug23.RData")
load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_libsize_19Aug23.RData")
load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP_libsize_19Aug23.RData")

##----------------------------WGS----------------------------##

metaFiltered_WGS <- metadataSamplesAllQC_HiSeq_WGS[,c("sample_type",
                                                      "disease_type",
                                                      "data_submitting_center_label")]
pvcaThreshold <- 0.7

source("00-functions.R") # for PVCA() function
pvca_WIS_WGS_Raw <- PVCA(counts = t(log(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS+0.5)),
                        meta = metaFiltered_WGS,
                        threshold = pvcaThreshold,
                        inter = FALSE)

pvca_WIS_WGS_VSNM <- PVCA(counts = t(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM),
                 meta = metaFiltered_WGS,
                 threshold = pvcaThreshold,
                 inter = FALSE)

pvca_WIS_WGS_CBS <- PVCA(counts = t(log(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CBS+0.5)),
                 meta = metaFiltered_WGS,
                 threshold = pvcaThreshold,
                 inter = FALSE)

pvcaWIS_WGS_CQ <- PVCA(counts = t(log(as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ)+0.5)),
                        meta = metaFiltered_WGS,
                        threshold = pvcaThreshold,
                        inter = FALSE)

pvca_WIS_WGS_CQP <- PVCA(counts = t(log(as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP)+0.5)),
                       meta = metaFiltered_WGS,
                       threshold = pvcaThreshold,
                       inter = FALSE)

pvca_WIS_WGS_CQ_libsize <- PVCA(counts = t(log(as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ_libsize)+0.5)),
                        meta = metaFiltered_WGS,
                        threshold = pvcaThreshold,
                        inter = FALSE)

pvca_WIS_WGS_CQP_libsize <- PVCA(counts = t(log(as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQP_libsize)+0.5)),
                        meta = metaFiltered_WGS,
                        threshold = pvcaThreshold,
                        inter = FALSE)

save(pvca_WIS_WGS_Raw, 
     pvca_WIS_WGS_VSNM,
     pvca_WIS_WGS_CBS,
     pvcaWIS_WGS_CQ,
     pvca_WIS_WGS_CQP,
     pvca_WIS_WGS_CQ_libsize,
     pvca_WIS_WGS_CQP_libsize,
     file = "Interim_data/pvca_WIS_WGS_raw_vsnm_cbs_cq_21Aug23.RData")

pvca_WIS_WGS_Raw
pvca_WIS_WGS_VSNM
pvca_WIS_WGS_CBS
pvcaWIS_WGS_CQ
pvca_WIS_WGS_CQP
pvca_WIS_WGS_CQ_libsize
pvca_WIS_WGS_CQP_libsize

##----------------------------RNA----------------------------##

# Load CQ data
load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_19Aug23.RData")
load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP_19Aug23.RData")
load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_libsize_19Aug23.RData")
load("Interim_data/tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP_libsize_19Aug23.RData")

metaFiltered_RNA <- metadataSamplesAllQC_HiSeq_RNA[,c("sample_type",
                                                      "disease_type",
                                                      "data_submitting_center_label")]
pvcaThreshold <- 0.7

source("00-functions.R") # for PVCA() function
pvca_WIS_RNA_Raw <- PVCA(counts = t(log(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA+0.5)),
                         meta = metaFiltered_RNA,
                         threshold = pvcaThreshold,
                         inter = FALSE)

pvca_WIS_RNA_VSNM <- PVCA(counts = t(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM),
                          meta = metaFiltered_RNA,
                          threshold = pvcaThreshold,
                          inter = FALSE)

pvca_WIS_RNA_CBS <- PVCA(counts = t(log(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CBS+0.5)),
                         meta = metaFiltered_RNA,
                         threshold = pvcaThreshold,
                         inter = FALSE)

pvcaWIS_RNA_CQ <- PVCA(counts = t(log(as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ)+0.5)),
                       meta = metaFiltered_RNA,
                       threshold = pvcaThreshold,
                       inter = FALSE)

pvca_WIS_RNA_CQP <- PVCA(counts = t(log(as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP)+0.5)),
                         meta = metaFiltered_RNA,
                         threshold = pvcaThreshold,
                         inter = FALSE)

pvca_WIS_RNA_CQ_libsize <- PVCA(counts = t(log(as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ_libsize)+0.5)),
                                meta = metaFiltered_RNA,
                                threshold = pvcaThreshold,
                                inter = FALSE)

pvca_WIS_RNA_CQP_libsize <- PVCA(counts = t(log(as.data.frame(tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQP_libsize)+0.5)),
                                 meta = metaFiltered_RNA,
                                 threshold = pvcaThreshold,
                                 inter = FALSE)

save(pvca_WIS_RNA_Raw, 
     pvca_WIS_RNA_VSNM,
     pvca_WIS_RNA_CBS,
     pvcaWIS_RNA_CQ,
     pvca_WIS_RNA_CQP,
     pvca_WIS_RNA_CQ_libsize,
     pvca_WIS_RNA_CQP_libsize,
     file = "Interim_data/pvca_WIS_RNA_raw_vsnm_cbs_cq_21Aug23.RData")

pvca_WIS_RNA_Raw
pvca_WIS_RNA_VSNM
pvca_WIS_RNA_CBS
pvcaWIS_RNA_CQ
pvca_WIS_RNA_CQP
pvca_WIS_RNA_CQ_libsize
pvca_WIS_RNA_CQP_libsize

##-----------------------------------------------------------------##
##----------------------------Plot PVCA----------------------------##
##-----------------------------------------------------------------##
## NOTE: The above PVCAs finalized use and comparison of raw vs. VSNM vs. CQ data
## These are plotted below

require(ggpubr)
require(ggsci)
# WGS
load("Interim_data/pvca_WIS_WGS_raw_vsnm_cbs_cq_21Aug23.RData", verbose = TRUE)

pvcaDataWGS <- as.data.frame(rbind(pvca_WIS_WGS_Raw,
                                   pvca_WIS_WGS_VSNM,
                                   pvcaWIS_WGS_CQ)) %>%
  rownames_to_column("group")
pvcaDataWGS$group <- c("Raw Counts", "VSNM", "ConQuR")
pvcaDataWGS.melted <- melt(pvcaDataWGS, id.vars = "group")
pvcaDataWGS.melted$group <- factor(pvcaDataWGS.melted$group, 
                                levels = c("Raw Counts", "VSNM", "ConQuR"))

pvcaPlotWGS <- ggplot(pvcaDataWGS.melted, aes(x = variable, y = value, fill = group)) + 
  geom_bar(aes(fill = group), position = "dodge", stat = "identity") + 
  geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(x = "Technical & Biological Effects",
       y = "Weighted average proportion variance",
       title = "WGS samples: PVCA of batch effect correction procedures\n(Only using WIS-overlapping genera)") +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.3, 0.75)) + 
  scale_x_discrete(labels=c("sample_type" = "Sample Type", 
                            "disease_type" = "Disease Type", 
                            "data_submitting_center_label" = "Sequencing Center",
                            "resid" = "Residual \n(not explained by technical variation)")) +
  scale_fill_nejm(name = "Data types", labels = c("Raw Count Data", 
                                                  "VSNM Normalized Data",
                                                  "ConQuR Normalized Data"))
ggsave(plot = pvcaPlotWGS,
       filename = "Figures/PVCA-VSNM-vs-ConQuR-WGS-batch-correction-WIS-feature-subset.jpeg", 
       dpi = "retina",
       width = 10, height = 6, units = "in")

# RNA
load("Interim_data/pvca_WIS_RNA_raw_vsnm_cbs_cq_21Aug23.RData", verbose = TRUE)

pvcaDataRNA <- as.data.frame(rbind(pvca_WIS_RNA_Raw,
                                   pvca_WIS_RNA_VSNM,
                                   pvcaWIS_RNA_CQ)) %>%
  rownames_to_column("group")
pvcaDataRNA$group <- c("Raw Counts", "VSNM", "ConQuR")
pvcaDataRNA.melted <- melt(pvcaDataRNA, id.vars = "group")
pvcaDataRNA.melted$group <- factor(pvcaDataRNA.melted$group, 
                                   levels = c("Raw Counts", "VSNM", "ConQuR"))

pvcaPlotRNA <- ggplot(pvcaDataRNA.melted, aes(x = variable, y = value, fill = group)) + 
  geom_bar(aes(fill = group), position = "dodge", stat = "identity") + 
  geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(x = "Technical & Biological Effects",
       y = "Weighted average proportion variance",
       title = "RNA-Seq samples: PVCA of batch effect correction procedures\n(Only using WIS-overlapping genera)") +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.3, 0.75)) + 
  scale_x_discrete(labels=c("sample_type" = "Sample Type", 
                            "disease_type" = "Disease Type", 
                            "data_submitting_center_label" = "Sequencing Center",
                            "resid" = "Residual \n(not explained by technical variation)")) +
  scale_fill_nejm(name = "Data types", labels = c("Raw Count Data", 
                                                  "VSNM Normalized Data",
                                                  "ConQuR Normalized Data"))
ggsave(plot = pvcaPlotRNA,
       filename = "Figures/PVCA-VSNM-vs-ConQuR-RNA-batch-correction-WIS-feature-subset.jpeg", 
       dpi = "retina",
       width = 10, height = 6, units = "in")
