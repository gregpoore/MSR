#-----------------------------------------------------------------------------
# S07-ML-tcga-multiclass.R
# Copyright (c) 2022--, Greg Poore
# Purpose: Perform multiclass ML on TCGA PT and BDN samples using WIS-overlapping bacterial genera
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
load("../../Interim_data/tcga_WIS_hiseq_VSNM_CBS_CQ_31Aug23.RData", verbose = TRUE)
# load("data_for_multiclass_ml_tcga_wgs24July22.RData")

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

## RNA
# Primary Tumor
cqRNA_PT <- mlMulticlass(metaData = metadataSamplesAllQC_HiSeq_RNA,
                      countData = tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_CQ,
                      sampleType = "Primary Tumor", numKFold = 10, seqCenter = "RNA_CQ")
save(cqRNA_PT, file = "../Interim_data/multiclass_PT_results_RNA_CQ_30Sept23.RData")

## WGS
# Primary Tumor
cq_PT <- mlMulticlass(metaData = metadataSamplesAllQC_HiSeq_WGS,
                       countData = tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ,
                       sampleType = "Primary Tumor", numKFold = 10, seqCenter = "WGS_CQ")
save(cq_PT, file = "../Interim_data/multiclass_PT_results_wgs_CQ_28Sept23.RData")

# Blood Derived Normal
cq_BDN <- mlMulticlass(metaData = metadataSamplesAllQC_HiSeq_WGS,
                        countData = tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_CQ,
                        sampleType = "Blood Derived Normal", numKFold = 10, seqCenter = "WGS_CQ")
save(cq_BDN, file = "../Interim_data/multiclass_BDN_results_wgs_CQ_28Sept23.RData")

#------------------------------------------#
# Run for VSNM and all seq center data
#------------------------------------------#

## VSNM WGS
# Primary Tumor
vsnmMulti_PT <- mlMulticlass(metaData = metadataSamplesAllQC_HiSeq_WGS,
                        countData = tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM,
                        sampleType = "Primary Tumor", numKFold = 10, seqCenter = "VSNM")
save(vsnmMulti_PT, file = "../Interim_data/multiclass_PT_results_wgs_VSNM_28Sept23.RData")

# Blood Derived Normal
vsnmMulti_BDN <- mlMulticlass(metaData = metadataSamplesAllQC_HiSeq_WGS,
                         countData = tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_VSNM,
                         sampleType = "Blood Derived Normal", numKFold = 10, seqCenter = "VSNM")
save(vsnmMulti_BDN, file = "../Interim_data/multiclass_BDN_results_wgs_VSNM_28Sept23.RData")

## VSNM RNA
# Primary Tumor
vsnmMultiRNA_PT <- mlMulticlass(metaData = metadataSamplesAllQC_HiSeq_RNA,
                             countData = tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_VSNM,
                             sampleType = "Primary Tumor", numKFold = 10, seqCenter = "VSNM")
save(vsnmMultiRNA_PT, file = "../../Interim_data/multiclass_PT_results_RNA_VSNM_30Sept23.RData")





