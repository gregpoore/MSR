#-----------------------------------------------------------------------------
# S09-ML-tcga-multiclass.R
# Copyright (c) 2022--, Greg Poore
# Purpose: Perform multiclass ML on TCGA PT and BDN samples
#-----------------------------------------------------------------------------

# Note: This was run interactively on Barnacle2 using the following interactive instance parameters and conda environment r4
# srun --nodes=1 --cpus-per-task=16 --mem=64gb --time=24:00:00 --pty bash -i

#-------------------------------#
# Load dependencies
require(tidyr)
require(caret) # for model building
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(xgboost) # for machine learning
require(randomForest) # for machine learning
require(tibble) # for df operations
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("../../Interim_data/KU_T2T_Full_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../../Interim_data/KU_T2T_WIS_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../../Interim_data/KU_T2T_BIO_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../../Interim_data/KU_T2T_Filt_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("../../Interim_data/KU_T2T_BIOFilt_data_for_ml_tcga_by_seq_center_14Oct23.RData")

# Load ConQuR-corrected data
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_BCM.RData")
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ_CMS.RData")
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_BCM.RData")
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ_CMS.RData")
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_BCM.RData")
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ_CMS.RData")
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_BCM.RData")
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ_CMS.RData")
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_BCM.RData")
load("../../Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_RNA_CQ_CMS.RData")

mlMulticlass <- function(metaData = metadataSamplesAllQC_HiSeq_WGS,
                         countData = tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ,
                         sampleType = "Primary Tumor",
                         seqCenter = "AllSeqCenters",
                         modelType = "xgbTree",
                         numResampleIter = 1,
                         numKFold = 10){
  
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- gsub('([[:punct:]])|\\s+','',seqCenter)
  numSeed <- 42
  xgbGrid <- data.frame(nrounds = 10,
                        max_depth = 4,
                        eta = .1,
                        gamma = 0,
                        colsample_bytree = .7,
                        min_child_weight = 1,
                        subsample = .8)
  
  metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  countDataFilt <- countData[rownames(metaDataFilt),]
  
  # if(seqCenter %in% c("AllSeqCenters","VSNM","WGS")){
  #   metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  #   metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  #   countDataFilt <- countData[rownames(metaDataFilt),]
  # } else{
  #   print("Filtering by seqCenter...")
  #   metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% filter(data_submitting_center_label == seqCenter) %>% droplevels()
  #   metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  #   countDataFilt <- countData[rownames(metaDataFilt),]
  # }
  
  mlDataY <- metaDataFilt
  mlDataX <- countDataFilt
  
  trainX <- mlDataX
  trainY <- mlDataY[,"predY"]
  
  set.seed(numSeed) # have to restate seed again, as ctrl defines the cross validation sampling during training
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       sampling = "up",
                       summaryFunction = multiClassSummary,
                       classProbs = TRUE,
                       verboseIter = FALSE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  print("Now training model...")
  set.seed(numSeed)
  mlModel <- train(x = trainX,
                   y = trainY,
                   method = modelType,
                   preProcess = c("zv"),
                   nthread = 1,
                   trControl = trainControl(method = "repeatedcv",
                                            number = numKFold,
                                            repeats = numResampleIter,
                                            sampling = "up",
                                            summaryFunction = multiClassSummary,
                                            classProbs = TRUE,
                                            verboseIter = FALSE,
                                            savePredictions = TRUE,
                                            allowParallel=TRUE),
                   # metric = "ROC",
                   tuneGrid = xgbGrid,
                   verbose = FALSE)
  
  resPredFun <- function(mlModel){
    resPred <- mlModel$pred
    
    ## Split folds and calculate perf on each fold
    resPredSplit <- split(resPred, resPred$Resample)
    repX_perf <- list()
    for(zz in seq_along(resPredSplit)){
      resPredSingleRep <- resPredSplit[[zz]]
      predProbs <- resPredSingleRep
      multiClass <- resPredSingleRep
      rep_perfTmp <- multiClassSummary(multiClass, lev = levels(multiClass$obs))
      repX_perf[[zz]] <- data.frame(as.list(rep_perfTmp))
    }
    
    # SUMMARIZE MODEL PERFORMANCES
    rep_perf <- do.call(rbind, repX_perf)
    # print(rep_perf)
    return(rep_perf)
  }
  
  print("Obtaining performance values...")
  resPredAll <- mlModel$pred
  perfCombinedAll <-  multiClassSummary(resPredAll, lev = levels(resPredAll$obs))
  repPerfCombinedAll <- resPredFun(mlModel)
  
  ## Save predictions and perf
  # Preds
  baseFilename <- paste0("multiclassCV_",st,"_",sc,"_k",numKFold,"_modelType_",modelType)
  write.csv(resPredAll, file = paste0("pred",baseFilename,".csv"))
  save(resPredAll, file = paste0("pred",baseFilename,".RData"))
  # Overall perf
  write.csv(perfCombinedAll, file = paste0("perfCombinedAll",baseFilename,".csv"))
  save(perfCombinedAll, file = paste0("perfCombinedAll",baseFilename,".RData"))
  # Rep perf
  write.csv(resPredAll, file = paste0("repPerfCombinedAll",baseFilename,".csv"))
  save(resPredAll, file = paste0("repPerfCombinedAll",baseFilename,".RData"))
  
  res <- list(resPredAll=resPredAll,
              perfCombinedAll=perfCombinedAll,
              repPerfCombinedAll=repPerfCombinedAll)
  return(res)
}

require(splitstackshape)

# NOTE: It may be necessary to run these one at a time (the multi-threading can appear to hang otherwise)

##----------------------WIS----------------------##
# WGS PT
# ku_cqWGS_PT_WIS <- mlMulticlass(metaData = metaKUT2TFinalNonzero_WIS_HiSeq_WGS,
#                             countData = kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ,
#                             sampleType = "Primary Tumor", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_PT_WIS, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_PT_WIS_15Oct23.RData")

# # WGS BDN
# ku_cqWGS_BDN_WIS <- mlMulticlass(metaData = metaKUT2TFinalNonzero_WIS_HiSeq_WGS,
#                              countData = kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ,
#                              sampleType = "Blood Derived Normal", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_BDN_WIS, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_WIS_15Oct23.RData")

# # RNA
# ku_cqRNA_PT_WIS <- mlMulticlass(metaData = metaKUT2TFinalNonzero_WIS_HiSeq_RNA,
#                              countData = kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "RNA_CQ")
# save(ku_cqRNA_PT_WIS, file = "../../Interim_data/multiclass_KU_T2T_CQ_RNA_PT_WIS_15Oct23.RData")

##----------------------BIO----------------------##
# # WGS PT
# ku_cqWGS_PT_BIO <- mlMulticlass(metaData = metaKUT2TFinalNonzero_BIO_HiSeq_WGS,
#                              countData = kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_PT_BIO, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_PT_BIO_15Oct23.RData")

# WGS BDN
# ku_cqWGS_BDN_BIO <- mlMulticlass(metaData = metaKUT2TFinalNonzero_BIO_HiSeq_WGS,
#                                   countData = kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ,
#                                   sampleType = "Blood Derived Normal", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_BDN_BIO, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_BIO_15Oct23.RData")

# # RNA PT
# ku_cqRNA_PT_BIO <- mlMulticlass(metaData = metaKUT2TFinalNonzero_BIO_HiSeq_RNA,
#                              countData = kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "RNA_CQ")
# save(ku_cqRNA_PT_BIO, file = "../../Interim_data/multiclass_KU_T2T_CQ_RNA_PT_BIO_15Oct23.RData")

# # CMS PT
# ku_CMS_PT_BIO <- mlMulticlass(metaData = metaKUT2TFinalNonzero_BIO_HiSeq_CMS,
#                              countData = kuT2TFinalNonzero_BIO_HiSeq_CMS,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "CMS_CQ")
# save(ku_CMS_PT_BIO, file = "../../Interim_data/multiclass_KU_T2T_CMSonly_RNA_PT_BIO_21Oct23.RData")

##----------------------Filt----------------------##
# # WGS PT
# ku_cqWGS_PT_Filt <- mlMulticlass(metaData = metaKUT2TFinalNonzero_Filt_HiSeq_WGS,
#                              countData = kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_PT_Filt, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_PT_Filt_15Oct23.RData")

# # WGS BDN
# ku_cqWGS_BDN_Filt <- mlMulticlass(metaData = metaKUT2TFinalNonzero_Filt_HiSeq_WGS,
#                                   countData = kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ,
#                                   sampleType = "Blood Derived Normal", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_BDN_Filt, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_Filt_15Oct23.RData")

# # RNA PT
# ku_cqRNA_PT_Filt <- mlMulticlass(metaData = metaKUT2TFinalNonzero_Filt_HiSeq_RNA,
#                              countData = kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "RNA_CQ")
# save(ku_cqRNA_PT_Filt, file = "../../Interim_data/multiclass_KU_T2T_CQ_RNA_PT_Filt_15Oct23.RData")

##----------------------BIOFilt----------------------##
# # WGS PT
# ku_cqWGS_PT_BIOFilt <- mlMulticlass(metaData = metaKUT2TFinalNonzero_BIOFilt_HiSeq_WGS,
#                              countData = kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_PT_BIOFilt, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_PT_BIOFilt_15Oct23.RData")

# # WGS BDN
# ku_cqWGS_BDN_BIOFilt <- mlMulticlass(metaData = metaKUT2TFinalNonzero_BIOFilt_HiSeq_WGS,
#                                   countData = kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ,
#                                   sampleType = "Blood Derived Normal", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_BDN_BIOFilt, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_BIOFilt_15Oct23.RData")

# # RNA PT
# ku_cqRNA_PT_BIOFilt <- mlMulticlass(metaData = metaKUT2TFinalNonzero_BIOFilt_HiSeq_RNA,
#                              countData = kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "RNA_CQ")
# save(ku_cqRNA_PT_BIOFilt, file = "../../Interim_data/multiclass_KU_T2T_CQ_RNA_PT_BIOFilt_15Oct23.RData")

##----------------------Full----------------------##

# # WGS PT
# ku_cqWGS_PT_Full <- mlMulticlass(metaData = metaKUT2TFinalNonzero_Full_HiSeq_WGS,
#                              countData = kuT2TFinalNonzero_Full_HiSeq_WGS_CQ,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_PT_Full, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_PT_Full_15Oct23.RData")

# # WGS BDN
# ku_cqWGS_BDN_Full <- mlMulticlass(metaData = metaKUT2TFinalNonzero_Full_HiSeq_WGS,
#                                   countData = kuT2TFinalNonzero_Full_HiSeq_WGS_CQ,
#                                   sampleType = "Blood Derived Normal", numKFold = 10, seqCenter = "WGS_CQ")
# save(ku_cqWGS_BDN_Full, file = "../../Interim_data/multiclass_KU_T2T_CQ_WGS_BDN_Full_15Oct23.RData")

# # RNA PT
# ku_cqRNA_PT_Full <- mlMulticlass(metaData = metaKUT2TFinalNonzero_Full_HiSeq_RNA,
#                              countData = kuT2TFinalNonzero_Full_HiSeq_RNA_CQ,
#                              sampleType = "Primary Tumor", numKFold = 10, seqCenter = "RNA_CQ")
# save(ku_cqRNA_PT_Full, file = "../../Interim_data/multiclass_KU_T2T_CQ_RNA_PT_Full_15Oct23.RData")




