#-----------------------------------------------------------------------------
# S02-ML-train-test-tcga-wis-feature-subset-per-seqcenter-and-vsnm.R
# Copyright (c) 2022--, Greg Poore
# Purpose: Run machine learning on all available TCGA cancer types by seq center
#-----------------------------------------------------------------------------

#-------------------------------#
# Load dependencies
require(splitstackshape)
require(reshape2)
require(tidyr)
require(caret) # for model building
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(tibble) # for df operations
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning
require(cvAUC)

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("data_for_ml_tcga_by_seq_center_and_experimental_strategy25July22.RData")

# GBM HYPERPARAMETER SEARCH GRID BELOW (DEFAULT PER CARET PACKAGE)
kfoldGBMGrid <- data.frame(n.trees=150, interaction.depth=3,
                                         shrinkage=0.1,
                                         n.minobsinnode=1)

# Model setup -- MODIFY AS NEEDED; ALL THESE COMPARISONS WILL BE TESTED
sampleTypeList <- c("Primary Tumor vs Solid Tissue Normal",
                    "Blood Derived Normal",
                    "Primary Tumor")
# MODIFY AS NEEDED
datasetList <- list(# Subset raw count data
  tcgaGenusKrakenQCFiltWIS_HMS,
  tcgaGenusKrakenQCFiltWIS_BCM,
  tcgaGenusKrakenQCFiltWIS_MDA,
  tcgaGenusKrakenQCFiltWIS_WashU,
  tcgaGenusKrakenQCFiltWIS_Broad_WGS,
  tcgaGenusKrakenQCFiltWIS_UNC,
  tcgaGenusKrakenQCFiltWIS_CMS,
  tcgaGenusKrakenQCFiltWIS_Broad_RNA,
  # VSNM batch corrected data
  vsnmDataGenusKrakenQCFiltWIS)

datasetListNames <- c(# Subset raw count data
  "tcgaGenusKrakenQCFiltWIS_HMS",
  "tcgaGenusKrakenQCFiltWIS_BCM",
  "tcgaGenusKrakenQCFiltWIS_MDA",
  "tcgaGenusKrakenQCFiltWIS_WashU",
  "tcgaGenusKrakenQCFiltWIS_Broad_WGS",
  "tcgaGenusKrakenQCFiltWIS_UNC",
  "tcgaGenusKrakenQCFiltWIS_CMS",
  "tcgaGenusKrakenQCFiltWIS_Broad_RNA",
  # VSNM batch corrected data
  "vsnmDataGenusKrakenQCFiltWIS")

metadataList <- list(# Subset metadata
  metadataSamplesAllQC_HiSeq_HMS,
  metadataSamplesAllQC_HiSeq_BCM,
  metadataSamplesAllQC_HiSeq_MDA,
  metadataSamplesAllQC_HiSeq_WashU,
  metadataSamplesAllQC_HiSeq_Broad_WGS,
  metadataSamplesAllQC_HiSeq_UNC,
  metadataSamplesAllQC_HiSeq_CMS,
  metadataSamplesAllQC_HiSeq_Broad_RNA,
  # Full metadata
  metadataSamplesAllQC)

metadataListNames <- list(# Subset metadata
  "metadataSamplesAllQC_HiSeq_HMS",
  "metadataSamplesAllQC_HiSeq_BCM",
  "metadataSamplesAllQC_HiSeq_MDA",
  "metadataSamplesAllQC_HiSeq_WashU",
  "metadataSamplesAllQC_HiSeq_Broad_WGS",
  "metadataSamplesAllQC_HiSeq_UNC",
  "metadataSamplesAllQC_HiSeq_CMS",
  "metadataSamplesAllQC_HiSeq_Broad_RNA",
  # Full metadata
  "metadataSamplesAllQC")

baseNamePerDatasetResultsFile <- "perfML_train_test_tcga_wis_feature_subset_seqcenter_vsnm_25July22"
baseNameAllResultsFile <- "perfML_train_test_tcga_wis_feature_subset_seqcenter_vsnm_ALL_25July22"

trainProp <- 0.7
numIter <- 1

# MATCHED METADATA DATA FRAME FOR ALL COUNT DATA:
caretTuneGrid <- kfoldGBMGrid
numKFold <- 4
numResampleIter <- 1
prroc_roc <- list()
prroc_pr <- list()
perf <- list()
rep_perf <- list()
perfTmp <- list()
rep_perfTmp <- list()
perfTmp2 <- list()
rep_perfTmp2 <- list()

for(jj in seq_along(datasetList)){

  dataTmp <- datasetList[[jj]]
  datasetName <- datasetListNames[[jj]]

  metaTmpQC <- metadataList[[jj]]
  metaTmpQC$disease_type <- factor(metaTmpQC$disease_type)
  metadataName <- metadataListNames[[jj]]

  for(kk in seq_along(sampleTypeList)){
    st <- sampleTypeList[kk]
    print(st)
    
    if(st == "Primary Tumor vs Solid Tissue Normal"){
      metaTmp <- metaTmpQC
      metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% c("Primary Tumor",
                                                                "Solid Tissue Normal"),])
    } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
      metaTmp <- metaTmpQC
      metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% st,])
    } else if(st == "Stage I vs IV"){
      metaTmp <- metaTmpPath
      metaTmp2 <- droplevels(metaTmp[(metaTmp$sample_type %in% "Primary Tumor") &
                                       (metaTmp$pathologic_stage_label_binned %in% c("Stage1","Stage4")),])
    }
    
    print(seq_along(levels(metaTmp2$disease_type)))
    
    for(ii in seq_along(levels(metaTmp2$disease_type))){
      
      dt <- levels(metaTmp2$disease_type)[ii]
      print(dt)
      
      if(st == "Primary Tumor vs Solid Tissue Normal"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- factor(gsub('([[:punct:]])|\\s+','',metaTmp3$sample_type))
        positiveClass <- "PrimaryTumor"
        negativeClass <- "SolidTissueNormal"
      } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
        metaTmp3 <- metaTmp2
        metaTmp3$predY <- factor(ifelse(metaTmp2$disease_type == dt, 
                                        yes = dt, 
                                        no = "OtherCancerType"),
                                 levels = c(dt, "OtherCancerType"))
        positiveClass <- gsub('([[:punct:]])|\\s+','',dt)
        negativeClass <- "OtherCancerType"
      } else if(st == "Stage I vs IV"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- metaTmp3$pathologic_stage_label_binned
        positiveClass <- "Stage4"
        negativeClass <- "Stage1"
      }
      
      print(table(metaTmp3$predY))
      
      # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
      if(length(table(metaTmp3$predY)) < 2){next}
      
      # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 20 SAMPLES IN EITHER CLASS
      if(any(table(metaTmp3$predY) < 20)){next}
      
      minorityClassSize <- min(table((metaTmp3$predY)))
      majorityClassSize <- max(table((metaTmp3$predY)))
      
      minorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == min(table(metaTmp3$predY)))])
      majorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == max(table(metaTmp3$predY)))])
      
      mlDataY <- metaTmp3
      mlDataX <- dataTmp[rownames(mlDataY),]

      repX_perf <- list()
      for(pp in 1:numIter){
        # USE 90% OF DATA FOR TRAINING AND 10% FOR TESTING
        set.seed(pp)
        index <- createDataPartition(mlDataY$predY, p = trainProp, list = FALSE)
        trainX <- mlDataX[index,]
        trainY <- mlDataY[index,"predY"]
        refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))

        testX <- mlDataX[-index,]
        testY <- mlDataY[-index,"predY"]
        refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
        
        set.seed(pp)
        ctrl <- trainControl(method = "repeatedcv",
                             number = numKFold,
                             repeats = numResampleIter,
                             sampling = "up",
                             summaryFunction = twoClassSummary,
                             classProbs = TRUE,
                             # verboseIter = TRUE,
                             savePredictions = TRUE,
                             allowParallel=TRUE)
        
        mlModel <- train(x = trainX,
                         y = refactoredTrainY,
                         method = "gbm",
                         # preProcess = c("scale","center"),
                         trControl = ctrl,
                         # verbose = TRUE,
                         metric = "ROC",
                         tuneGrid = caretTuneGrid)

        predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
        fg <- predProbs[refactoredTestY == positiveClass]
        bg <- predProbs[refactoredTestY == negativeClass]
        
        rep_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
        rep_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
        
        multiClass <- data.frame(obs = refactoredTestY,
                                 pred = predict(mlModel, newdata = testX),
                                 predict(mlModel, newdata = testX, type = "prob"))
        repX_perf[[pp]] <- data.frame(auroc=rep_roc$auc,
                               aupr=rep_pr$auc.integral,
                               rep=paste0("Fold",pp), 
                               diseaseType = dt,
                               sampleType = st,
                               datasetName = datasetName,
                               metadataName = metadataName,
                               minorityClassSize = minorityClassSize,
                               majorityClassSize = majorityClassSize,
                               minorityClassName = minorityClassName,
                               majorityClassName = majorityClassName)

        rm(mlModel)
      }
      rep_perf[[ii]] <- do.call(rbind, repX_perf)
      
      #--------------------------------------#
      # Save performance into relevant files #
      #--------------------------------------#
      
      filepathPerfStatsPerFold <- paste0("./perfDataPerFold__",datasetName)
      
      if(!( dir.exists( file.path(filepathPerfStatsPerFold)))){
        dir.create(file.path(filepathPerfStatsPerFold))
      }

      filenamePerFoldCSV <- paste0(filepathPerfStatsPerFold,"/",
                            dt,
                            " -- ",
                            st,
                            " -- PerfPerFold.csv")
      
      write.csv(rep_perf[[ii]], file = filenamePerFoldCSV)
      
    }
    # SUMMARIZE RESULTS
    rep_perfTmp[[kk]] <- do.call(rbind, rep_perf)
  }
  # SUMMARIZE RESULTS
  rep_perfTmp2[[jj]] <- do.call(rbind, rep_perfTmp)

  write.csv(rep_perfTmp2[[jj]], file = paste0("rep_",baseNamePerDatasetResultsFile,"_",datasetName,".csv"))
}

# SUMMARIZE RESULTS
rep_perf1VsAll <- do.call(rbind, rep_perfTmp2)

write.csv(rep_perf1VsAll, file = paste0("rep_",baseNameAllResultsFile,".csv"))

#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------