#-----------------------------------------------------------------------------
# 00-functions.R
# Copyright (c) 2021--, Greg Poore
# Purpose: List main functions for R scripts
#-----------------------------------------------------------------------------

plotControlsRaw <- function(seqCenter="HMS", 
                            inputSampleType="Primary Tumor",
                            qvalSize = 2.5,
                            statSpacingROC=0.5,
                            statSpacingPR=0.25,
                            alphaVal=0.2,
                            tipLength=0.01,
                            qvalAsterisks=FALSE){
  require(rstatix)
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Primary Tumor vs Solid Tissue Normal"){st = "PT_vs_NAT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  ## AUROC
  mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
    filter(sampleType == inputSampleType) %>%
    filter(grepl(seqCenter,datasetNameShort)) %>%
    distinct() %>% droplevels() %>%
    wilcox_test(AUROC ~ datasetNameShort, exact = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    mutate(p.adj = signif(p.adj, digits=3)) %>%
    add_xy_position(x = "datasetNameShort", step.increase = statSpacingROC) -> roc.stat.test
    
  if(all(roc.stat.test$y.position == 1)){
    roc.stat.test$y.position <- seq(from=1.05, by=0.2, length.out = length(roc.stat.test$y.position))
  }
  # print(data.frame(roc.stat.test))
  
  mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
    filter(sampleType == inputSampleType) %>%
    filter(grepl(seqCenter,datasetNameShort)) %>%
    distinct() %>% droplevels() %>%
    ggboxplot(x = "datasetNameShort",
             y = "AUROC", fill = "datasetNameShort",
             legend = "none",
             notch = TRUE,
             xlab = "",
             add = "jitter",
             add.params = list(alpha=alphaVal),
             palette = "nejm") +
    rotate_x_text(90) + 
    stat_pvalue_manual(roc.stat.test, 
                       # label = "q = {p.adj}", 
                       label = ifelse(qvalAsterisks, 
                                      yes = "{p.adj.signif}", 
                                      no = "q = {p.adj}"),
                       tip.length = tipLength,
                       size = qvalSize) +
    geom_hline(yintercept = 0.5, linetype="dotted") + 
    labs(fill = "datasetNameShort") + 
    theme(plot.title = element_text(hjust=0.5)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
                       limits = c(0,1.01*max(roc.stat.test$y.position))) -> plotROC
  
  ## AUPR
  mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
    filter(sampleType == inputSampleType) %>%
    filter(grepl(seqCenter,datasetNameShort)) %>%
    distinct() %>% droplevels() %>%
    wilcox_test(AUPR ~ datasetNameShort, exact = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    mutate(p.adj = signif(p.adj, digits=3)) %>%
    add_xy_position(x = "datasetNameShort", step.increase = statSpacingPR) -> pr.stat.test
  
  if(all(pr.stat.test$y.position == 1)){
    pr.stat.test$y.position <- seq(from=1.05, by=0.2, length.out = length(pr.stat.test$y.position))
  }
  # print(data.frame(pr.stat.test))
  
  mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
    filter(sampleType == inputSampleType) %>%
    filter(grepl(seqCenter,datasetNameShort)) %>%
    distinct() %>% droplevels() %>%
    ggboxplot(x = "datasetNameShort",
             y = "AUPR", fill = "datasetNameShort",
             legend = "none",
             notch = TRUE,
             add = "jitter",
             add.params = list(alpha=alphaVal),
             xlab = "",
             palette = "nejm") +
    rotate_x_text(90) + 
    stat_pvalue_manual(pr.stat.test, 
                       # label = "q = {p.adj}", 
                       label = "{p.adj.signif}", 
                       tip.length = tipLength,
                       size = qvalSize) +
    labs(fill = "datasetNameShort") + 
    theme(plot.title = element_text(hjust=0.5)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
                       limits = c(0,1.01*max(pr.stat.test$y.position))) -> plotPR
  
  combinedPlot <- ggarrange(plotROC, plotPR, ncol = 2) 
  combinedPlotAnnotated <- annotate_figure(combinedPlot, top = text_grob(paste0(seqCenter," | ",inputSampleType,"\nActual vs. Controls"), 
                                        color = "black", face = "bold", size = 14))
  print(combinedPlotAnnotated)
  baseName <- ifelse(qvalAsterisks, yes = "control_v2_signif_scrambled_overlay_raw_data",
                     no = "control_v2_scrambled_overlay_raw_data_")
  ggsave(filename = paste0("Figures/",baseName,"_",seqCenter,"_",st,".svg"),
           plot = combinedPlotAnnotated,
           dpi = "retina", units = "in", width = 6, height = 7)
}


cmPredHeatmapTCGA <- function(predsPT, predsBDN, plotString, midPT = 80, midBDN = 80, midAuto = FALSE,
                          ptOnlyFlag = FALSE, dimPT = 8, dimBDN = 8){
  require(caret)
  require(tidyr)
  # Generate heatmaps of confusion matrices
  # Source: https://stackoverflow.com/questions/7421503/how-to-plot-a-confusion-matrix-using-heatmaps-in-r
  # RColorBrewer::display.brewer.all()

  cm_PT <- confusionMatrix(predsPT$pred, predsPT$obs)

  if(midAuto){
    midPT <- round(0.70*max(cm_PT$table))
    print(sprintf("Auto PT midpoint: %d", midPT))
  }
  
  ## Primary Tumor
  hm_pt_df <- as.data.frame.matrix(cm_PT$table) %>% rownames_to_column("y")
  hm_pt <- tibble(hm_pt_df)
  # print(hm_pt)
  hm_pt <- hm_pt %>% gather(x, value, -y)
  # hm <- cm_mlMulticlassPredPT_df
  hm_pt <- hm_pt %>%
    mutate(x = factor(x), # alphabetical order by default
           y = factor(y, levels = rev(unique(y)))) %>% # force reverse alphabetical order
    mutate(diagFill = ifelse(y==x,max(hm_pt$value),value)) # force the diagonal to have the highest fill value (ie red/orange)

  # Plot PT
  p_pt <- ggplot(hm_pt, aes(x=x, y=y, fill=diagFill)) +
    geom_tile(color="black") + theme_bw() + coord_equal() +
    theme_pubr() + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_fill_gradient2(low = "white",mid = "orange", high = "red", midpoint = midPT) +
    # scale_fill_gradient2(low = "white",mid = "yellow", high = "orange", midpoint = midPT) +
    # scale_fill_distiller(palette="RdYlBu", direction=-1) + # or YlOrRd
    guides(fill="none") + # removing legend for `fill`
    rotate_x_text(30) +
    # scale_x_discrete(position = "top") +
    labs(title = paste0(plotString," | Primary Tumor | TCGA Bins Full | Multiclass ML"),
         x = "Reference", y = "Predicted") + # using a title instead
    geom_text(aes(label=value), color="black") # printing values

  if(ptOnlyFlag){

    print("PT stats:")
    print(cm_PT$overall)

    print(p_pt)
    ggsave(paste0("Figures/mlMulticlass_tcga_",plotString,"_xgboost_PT.jpeg"), dpi = "retina",
           width = dimPT, height = dimPT, units = "in")
    
    res <- list(hm_pt=hm_pt,
                p_pt=p_pt)
    
  } else{
    cm_BDN <- confusionMatrix(predsBDN$pred, predsBDN$obs)
    if(midAuto){
      midBDN <- round(0.70*max(cm_BDN$table))
      print(sprintf("Auto BDN midpoint: %d", midBDN))
    }

    ## Blood Derived Normal
    hm_bdn_df <- as.data.frame.matrix(cm_BDN$table) %>% rownames_to_column("y")
    hm_bdn <- tibble(hm_bdn_df)
    hm_bdn <- hm_bdn %>% gather(x, value, -y)
    # hm <- cm_mlMulticlassPredPT_df
    hm_bdn <- hm_bdn %>%
      mutate(x = factor(x), # alphabetical order by default
             y = factor(y, levels = rev(unique(y)))) %>% # force reverse alphabetical order
      mutate(diagFill = ifelse(y==x,max(hm_bdn$value),value)) # force the diagonal to have the highest fill value (ie red/orange)

    # Plot BDN
    p_bdn <- ggplot(hm_bdn, aes(x=x, y=y, fill=diagFill)) +
      geom_tile(color="black") + theme_bw() + coord_equal() +
      theme_pubr() + theme(plot.title = element_text(hjust = 0.5)) + 
      scale_fill_gradient2(low = "white",mid = "orange", high = "red", midpoint = midBDN) +
      # scale_fill_distiller(palette="RdYlBu", direction=-1) + # or YlOrRd
      guides(fill="none") + # removing legend for `fill`
      rotate_x_text(30) +
      # scale_x_discrete(position = "top") +
      labs(title = paste0(plotString," | Blood Derived Normal | TCGA Bins Full | Multiclass ML"),
           x = "Reference", y = "Predicted") + # using a title instead
      geom_text(aes(label=value), color="black") # printing values

    # Print stats
    print("BDN stats:")
    print(cm_BDN$overall)
    cat('\n')
    print("PT stats:")
    print(cm_PT$overall)
    
    print(p_bdn)
    ggsave(paste0("Figures/tcga_diag_mlMulticlass_",plotString,"_xgboost_BDN.jpeg"), dpi = "retina",
           width = dimBDN, height = dimBDN, units = "in")
    print(p_pt)
    ggsave(paste0("Figures/tcga_diag_mlMulticlass_",plotString,"_xgboost_PT.jpeg"), dpi = "retina",
           width = dimPT, height = dimPT, units = "in")
    
    res <- list(hm_pt=hm_pt,
                hm_bdn=hm_bdn,
                cm_BDN=cm_BDN,
                cm_PT=cm_PT,
                p_pt=p_pt,
                p_bdn=p_bdn)
  }
  
  return(res)
}


PVCA <- function(counts, meta, threshold, inter){
  # Not my function. Copied from NIEHS:
  # https://www.niehs.nih.gov/research/resources/software/biostatistics/pvca/index.cfm
  
  counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
  cor.counts <- cor(counts.center)
  dim(cor.counts)
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val/eigen.val.sum
  meta <- as.data.frame(meta)
  
  all <- 0
  npc.in <- 0
  for(i in 1:n.eigen){
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if(all > threshold){break}
  }
  if (npc.in < 3) {npc <- 3}
  
  pred.list <- colnames(meta)
  meta <- droplevels(meta)
  
  n.preds <- ncol(meta) + 1
  if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}
  
  ran.pred.list <- c()
  for(i in 1:ncol(meta)){
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
  }
  ##interactions
  if(inter){
    for(i in 1:(ncol(meta)-1)){
      for(j in (i+1):ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep=" ~ ")
  ran.var.mat <- NULL
  for(i in 1:npc.in){
    dat <- cbind(eigen.mat[,i],meta)
    colnames(dat) <- c("pc",colnames(meta))
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit)
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
  wgt.vec <- eigen.val/eigen.val.sum
  prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
  std.prop.var <- prop.var/sum(prop.var)
  std.prop.var
}


barplotPerf <- function(inputData=mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter, 
                        sampleTypeInput, seqCenterAbbrev, 
                        plotWidthSingle = 6,
                        plotWidthCombined = 10,
                        intFlag = FALSE,
                        ciWidth = 0.99,
                        singlePlotHeight = 4,
                        combPlotHeight = 4,
                        fileNameString = "vsnmVsRaw",
                        factorCQ = FALSE){
  inputDataFilt <- inputData %>%
    filter(sampleType == sampleTypeInput) %>%
    filter(grepl(seqCenterAbbrev,datasetName)) %>%
    reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName",
                              "metadataName","minorityClassSize","majorityClassSize",
                              "minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
    summarySE(measurevar = "value", 
              groupvars = c("datasetName","metadataName","variable",
                            "abbrev","minorityClassSize","majorityClassSize",
                            "minorityClassName","majorityClassName",
                            "nullAUPR","nullAUROC"),
              conf.interval = ciWidth) %>%
    mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
    mutate(datasetName = gsub("vb.+","Raw",datasetName)) %>%
    mutate(datasetName = gsub("snm.+","VSNM",datasetName)) %>%
    mutate(datasetName = gsub("tcgaGenusKrakenAllFiltWIS_HiSeq_[WGS|RNA]_","",datasetName))

    if(factorCQ){
      inputDataFilt <- inputDataFilt %>%
      mutate(datasetName = factor(case_when(
        grepl("CQ",datasetName) ~ "ConQuR",
        grepl("VSNM",datasetName) ~ "VSNM",
        grepl(seqCenterAbbrev,datasetName) ~ "Raw"
        ), levels = c("Raw","VSNM","ConQuR")))
    }
  
  inputDataFilt %>%
    filter(variable == "AUROC") %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_manual(values = c("#ADB6B6FF","#925E9FFF","#E18727FF"), name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUROC",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUROC
  
  if(intFlag){
    fileNameAUROC <- paste0("Figures/barplot_AUROC_tcga_INT_",fileNameString,"_",seqCenterAbbrev,"_",
                            gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                            ".jpeg")
  } else{
    fileNameAUROC <- paste0("Figures/barplot_AUROC_tcga_",fileNameString,"_",seqCenterAbbrev,"_",
                            gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                            ".jpeg")
  }
  ggsave(filename = fileNameAUROC,
         plot = barPlotAUROC,
         dpi = "retina", units = "in", height = singlePlotHeight, width = plotWidthSingle)
  
  inputDataFilt %>%
    filter(variable == "AUPR") %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_manual(values = c("#ADB6B6FF","#925E9FFF","#E18727FF"), name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUPR",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUPR
  
  if(intFlag){
    fileNameAUPR <- paste0("Figures/barplot_AUPR_tcga_INT_",fileNameString,"_",seqCenterAbbrev,"_",
                           gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                           ".jpeg")
  } else{
    fileNameAUPR <- paste0("Figures/barplot_AUPR_tcga_",fileNameString,"_",seqCenterAbbrev,"_",
                           gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                           ".jpeg")
  }
  ggsave(filename = fileNameAUPR,
         plot = barPlotAUPR,
         dpi = "retina", units = "in", height = singlePlotHeight, width = plotWidthSingle)
  
  inputDataFilt %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    facet_wrap(vars(variable)) +
    scale_fill_manual(values = c("#ADB6B6FF","#925E9FFF","#E18727FF"), name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "Area Under Curve",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotCombined
  
  print(barPlotCombined)
  
  if(intFlag){
    fileName <- paste0("Figures/barplot_combined_tcga_INT_",fileNameString,"_",seqCenterAbbrev,"_",
                       gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                       ".jpeg")
  } else{
    fileName <- paste0("Figures/barplot_combined_tcga_",fileNameString,"_",seqCenterAbbrev,"_",
                       gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                       ".jpeg")
  }
  ggsave(filename = fileName,
         plot = barPlotCombined,
         dpi = "retina", units = "in", height = combPlotHeight, width = plotWidthCombined)
}

# Make one summary barplot 
barplotSummaryPerf <- function(inputData=mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int, 
                        sampleTypeInput, 
                        seqCenterAbbrev="All", 
                        plotWidthSingle = 6,
                        plotWidthCombined = 10,
                        intFlag = FALSE,
                        ciWidth = 0.99,
                        fileNameString = "vsnmVsRaw",
                        factorCQ = FALSE){
  inputDataFilt <- inputData %>%
    filter(sampleType == sampleTypeInput) %>%
    mutate(metadataName = "All") %>%
    select(AUROC, AUPR, abbrev,diseaseType,
           sampleType,datasetName,metadataName) %>%
    mutate(datasetName = gsub("vb.+","Raw",datasetName)) %>%
    mutate(datasetName = gsub("snm.+","VSNM",datasetName)) %>%
    mutate(datasetName = gsub("tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_|tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_","",datasetName)) %>%
    mutate(datasetName = case_when(
      datasetName == "Raw" ~ "Raw",
      grepl("CQ",datasetName) ~ "ConQuR",
      grepl("VSNM",datasetName) ~ "VSNM",
      datasetName %in% c("HMS","BCM","MDA","WashU","Broad_WGS","UNC","CMS","Broad_RNA") ~ "Raw"
    )) %>%
    reshape2::melt(id.vars = c("abbrev","diseaseType","sampleType","datasetName",
                               "metadataName")) %>%
    summarySE(measurevar = "value", 
              groupvars = c("datasetName", "metadataName","variable","abbrev"),
              conf.interval = ciWidth)

  if(factorCQ){
      inputDataFilt <- inputDataFilt %>%
        mutate(datasetName = factor(datasetName, levels = c("Raw","VSNM","ConQuR")))
    }
  
  inputDataFilt %>%
    filter(variable == "AUROC") %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_manual(values = c("#ADB6B6FF","#925E9FFF","#E18727FF"), name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUROC",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUROC
  
  if(intFlag){
    fileNameAUROC <- paste0("Figures/barplot_AUROC_tcga_INT_",fileNameString,"_",seqCenterAbbrev,"_",
                            gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                            ".jpeg")
  } else{
    fileNameAUROC <- paste0("Figures/barplot_AUROC_tcga_",fileNameString,"_",seqCenterAbbrev,"_",
                            gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                            ".jpeg")
  }
  ggsave(filename = fileNameAUROC,
         plot = barPlotAUROC,
         dpi = "retina", units = "in", height = 3, width = plotWidthSingle)
  
  inputDataFilt %>%
    filter(variable == "AUPR") %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_manual(values = c("#ADB6B6FF","#925E9FFF","#E18727FF"), name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUPR",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUPR
  
  if(intFlag){
    fileNameAUPR <- paste0("Figures/barplot_AUPR_tcga_INT_",fileNameString,"_",seqCenterAbbrev,"_",
                           gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                           ".jpeg")
  } else{
    fileNameAUPR <- paste0("Figures/barplot_AUPR_tcga_",fileNameString,"_",seqCenterAbbrev,"_",
                           gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                           ".jpeg")
  }
  ggsave(filename = fileNameAUPR,
         plot = barPlotAUPR,
         dpi = "retina", units = "in", height = 3, width = plotWidthSingle)
  
  inputDataFilt %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    facet_wrap(vars(variable), nrow = 2) +
    scale_fill_manual(values = c("#ADB6B6FF","#925E9FFF","#E18727FF"), name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "Area Under Curve",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotCombined
  
  print(barPlotCombined)
  
  if(intFlag){
    fileName <- paste0("Figures/barplot_combined_tcga_INT_",fileNameString,"_",seqCenterAbbrev,"_",
                       gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                       ".jpeg")
  } else{
    fileName <- paste0("Figures/barplot_combined_tcga_",fileNameString,"_",seqCenterAbbrev,"_",
                       gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                       ".jpeg")
  }
  ggsave(filename = fileName,
         plot = barPlotCombined,
         dpi = "retina", units = "in", height = 5, width = plotWidthCombined)
  
  # return(inputDataFilt)
}


# Make one summary barplot 
barplotSummaryPerfNoHu <- function(inputData=mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter_Int, 
                        sampleTypeInput, 
                        seqCenterAbbrev="All", 
                        plotWidthSingle = 6,
                        plotWidthCombined = 10,
                        intFlag = FALSE,
                        ciWidth = 0.99,
                        prefix2Remove = "vbDataBarnDFReconciledQC_Genus_",
                        fileNameString = "vsnmVsRaw",
                        factorCQ = FALSE,
                        outputData = FALSE){
  inputDataFilt <- inputData %>%
    filter(sampleType == sampleTypeInput) %>%
    mutate(metadataName = "All") %>%
    select(AUROC, AUPR, abbrev,diseaseType,
           sampleType,datasetName,metadataName) %>%
    mutate(datasetName = gsub(prefix2Remove,"",datasetName)) %>%
    mutate(datasetName = gsub("HiSeq_","",datasetName)) %>%
    # mutate(datasetName = gsub("snm.+","VSNM",datasetName)) %>%
    mutate(datasetName = gsub("tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_|tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_","",datasetName)) %>%
    mutate(datasetName = case_when(
      grepl("NoHu",datasetName) ~ "NoHu",
      grepl("CQ",datasetName) ~ "ConQuR",
      grepl("VSNM",datasetName) ~ "VSNM",
      datasetName %in% c("HMS","BCM","MDA","WashU","Broad_WGS","UNC","CMS","Broad_RNA") ~ "Full"
    )) %>%
    reshape2::melt(id.vars = c("abbrev","diseaseType","sampleType","datasetName",
                               "metadataName")) %>%
    summarySE(measurevar = "value", 
              groupvars = c("datasetName", "metadataName","variable","abbrev"),
              conf.interval = ciWidth)

  if(factorCQ){
      inputDataFilt <- inputDataFilt %>%
        mutate(datasetName = factor(datasetName, levels = c("Raw","VSNM","ConQuR")))
    }
  
  inputDataFilt %>%
    filter(variable == "AUROC") %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_manual(values = c("#ADB6B6FF","#00468BFF"), name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUROC",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUROC
  
  if(intFlag){
    fileNameAUROC <- paste0("Figures/barplot_AUROC_tcga_INT_",fileNameString,"_",seqCenterAbbrev,"_",
                            gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                            ".jpeg")
  } else{
    fileNameAUROC <- paste0("Figures/barplot_AUROC_tcga_",fileNameString,"_",seqCenterAbbrev,"_",
                            gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                            ".jpeg")
  }
  ggsave(filename = fileNameAUROC,
         plot = barPlotAUROC,
         dpi = "retina", units = "in", height = 3, width = plotWidthSingle)
  
  inputDataFilt %>%
    filter(variable == "AUPR") %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_manual(values = c("#ADB6B6FF","#00468BFF"), name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUPR",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUPR
  
  if(intFlag){
    fileNameAUPR <- paste0("Figures/barplot_AUPR_tcga_INT_",fileNameString,"_",seqCenterAbbrev,"_",
                           gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                           ".jpeg")
  } else{
    fileNameAUPR <- paste0("Figures/barplot_AUPR_tcga_",fileNameString,"_",seqCenterAbbrev,"_",
                           gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                           ".jpeg")
  }
  ggsave(filename = fileNameAUPR,
         plot = barPlotAUPR,
         dpi = "retina", units = "in", height = 3, width = plotWidthSingle)
  
  inputDataFilt %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    facet_wrap(vars(variable), nrow = 2) +
    scale_fill_manual(values = c("#ADB6B6FF","#00468BFF"), name = "Data type") +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "Area Under Curve",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotCombined
  
  print(barPlotCombined)
  
  if(intFlag){
    fileName <- paste0("Figures/barplot_combined_tcga_INT_",fileNameString,"_",seqCenterAbbrev,"_",
                       gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                       ".jpeg")
  } else{
    fileName <- paste0("Figures/barplot_combined_tcga_",fileNameString,"_",seqCenterAbbrev,"_",
                       gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                       ".jpeg")
  }
  ggsave(filename = fileName,
         plot = barPlotCombined,
         dpi = "retina", units = "in", height = 3, width = plotWidthCombined)

  ## Identify cancer types within non-overlapping CIs
  # AUROC
  cmpAUROC <- inputDataFilt %>% filter(variable == "AUROC") %>%
    mutate(low = ifelse(value-ci<0,0,value-ci), 
           high = ifelse(value+ci>1,1,value+ci)) %>%
    select(-metadataName, -variable) %>%
    select(datasetName, abbrev, low, high) %>%
    pivot_wider(names_from = c("datasetName"),
                values_from = c("low","high")) %>%
    mutate(overlapNoHu = (between(low_NoHu,low_Full,high_Full) |
                            between(high_NoHu,low_Full,high_Full) ),
           overlapHu = (between(low_Full,low_NoHu,high_NoHu) |
                          between(high_Full,low_NoHu,high_NoHu) )) %>%
    mutate(ovComb = overlapNoHu | overlapHu) %>%
    select(-overlapNoHu, -overlapHu) %>%
    mutate(factorHigh = ifelse( ((low_Full+high_Full)/2) > ((low_NoHu+high_NoHu)/2), "Full","NoHu")) %>%
    filter(!ovComb)

  # AUROC
  cmpAUPR <- inputDataFilt %>% filter(variable == "AUPR") %>%
    mutate(low = ifelse(value-ci<0,0,value-ci), 
           high = ifelse(value+ci>1,1,value+ci)) %>%
    select(-metadataName, -variable) %>%
    select(datasetName, abbrev, low, high) %>%
    pivot_wider(names_from = c("datasetName"),
                values_from = c("low","high")) %>%
    mutate(overlapNoHu = (between(low_NoHu,low_Full,high_Full) |
                            between(high_NoHu,low_Full,high_Full) ),
           overlapHu = (between(low_Full,low_NoHu,high_NoHu) |
                          between(high_Full,low_NoHu,high_NoHu) )) %>%
    mutate(ovComb = overlapNoHu | overlapHu) %>%
    select(-overlapNoHu, -overlapHu) %>%
    mutate(factorHigh = ifelse( ((low_Full+high_Full)/2) > ((low_NoHu+high_NoHu)/2), "Full","NoHu")) %>%
    filter(!ovComb)

  if(nrow(cmpAUROC)==0){
    print("No differing AUROCs")
  } else{
    print("Differing AUROCs:")
    print(cmpAUROC)
  }
  if(nrow(cmpAUPR)==0){
    print("No differing AUPRs")
  } else{
    print("Differing AUPRs:")
    print(cmpAUPR)
  }
  
  if(outputData){
    return(inputDataFilt)
  }
}

#--------Expanded enrichment test function--------#

enrichmentFxn2 <- function(myPath = "Supporting_scripts/S10-ML-10k-tcga-raw-vsnm-subsets/features__",
                          totalFeatures = colnames(snmDataSampleTypeWithExpStrategy),
                          seqCenter = "HMS",
                          sampleType = "Primary Tumor",
                          dataSetRaw = "vbDataBarnDFReconciledQC_Int_HiSeq",
                          dataSetVSNM = "snmDataSampleTypeWithExpStrategy_HiSeq",
                          cancerAbbrevs = abbreviationsTCGA_Allcancer,
                          plotWidth = 5,
                          fileNameString = "vsnmVsRaw",
                          kendallOnlyFlag = FALSE,
                          showCM = FALSE){

  if(seqCenter == ""){
    dataSetRaw_SeqCenter <- dataSetRaw
    dataSetVSNM_SeqCenter <- dataSetVSNM

  } else{
    dataSetRaw_SeqCenter <- paste0(dataSetRaw,"_",seqCenter)
    dataSetVSNM_SeqCenter <- paste0(dataSetVSNM,"_",seqCenter)
  }
  
  
  
  # Get cancer types
  cancerTypes <- gsub(" --.+","",
                      grep(paste0(" -- ",sampleType," -- "),
                           list.files(paste0(myPath,dataSetRaw_SeqCenter)),value = TRUE))
  
  chi2List <- list()
  fisherList <- list()
  kendallList <- list()
  for(ii in 1:length(cancerTypes)){
    cancerType <- cancerTypes[ii]
    
    fileName <- paste0(cancerType," -- ",sampleType," -- Features.csv")
    
    featuresRaw <- read.csv(paste0(myPath,dataSetRaw_SeqCenter,"/",fileName)) %>%
      mutate(rank = row_number())
    featuresVSNM <- read.csv(paste0(myPath,dataSetVSNM_SeqCenter,"/",fileName)) %>%
      mutate(rank = row_number())
    
    print(sprintf("Length raw list: %d",length(featuresRaw$Taxa)))
    print(sprintf("Length VSNM list: %d",length(featuresVSNM$Taxa)))
    
    totalFeatures_NOT_VSNM <- totalFeatures[!(totalFeatures %in% featuresVSNM$Taxa)]
    totalFeatures_NOT_Raw <- totalFeatures[!(totalFeatures %in% featuresRaw$Taxa)]
    
    chi2df <- matrix(c(length(intersect(featuresRaw$Taxa,featuresVSNM$Taxa)),
                       sum( totalFeatures_NOT_VSNM %in% featuresRaw$Taxa),
                       sum( totalFeatures_NOT_Raw %in% featuresVSNM$Taxa),
                       sum(totalFeatures_NOT_VSNM %in% totalFeatures_NOT_Raw) ),
                     nrow = 2, ncol = 2)
    colnames(chi2df) <- c("Raw+","Raw-")
    rownames(chi2df) <- c("VSNM+","VSNM-")
    # Chi2 test
    chi2Test <- chisq.test(chi2df)
    chi2List[[ii]] <- chi2List
    if(showCM){
      print(cancerType)
      print(chi2df)
    }
    
    # Fisher exact test
    fisherTest <- fisher.test(chi2df)
    # print(fisherTest$estimate)
    # print(fisherTest$null.value)
    fisherList[[ii]] <- data.frame(OR = as.numeric(fisherTest$estimate),
                                   OR_low = fisherTest$conf.int[1],
                                   OR_high = fisherTest$conf.int[2],
                                   p = fisherTest$p.value,
                                   CT = cancerType,
                                   SC = seqCenter,
                                   ST = sampleType)
    
    #--------Kendall tau--------#
    featuresRawMissing <- data.frame(Taxa = totalFeatures_NOT_Raw,
                                     varImp = 0,
                                     rank = max(featuresRaw$rank)+1)
    featuresVSNMMissing <- data.frame(Taxa = totalFeatures_NOT_VSNM,
                                      varImp = 0,
                                      rank = max(featuresVSNM$rank)+1)
    featuresRawFull <- rbind(featuresRaw, featuresRawMissing)
    featuresVSNMFull <- rbind(featuresVSNM, featuresVSNMMissing)
    
    featuresRawVSNMFull <- featuresRawFull %>%
      left_join(featuresVSNMFull, by = "Taxa")
    
    kendallTest <- cor.test(x = featuresRawVSNMFull$varImp.x,
                            y = featuresRawVSNMFull$varImp.y,
                            method = "kendall")
    
    kendallList[[ii]] <- data.frame(tau = as.numeric(kendallTest$estimate),
                                   p = kendallTest$p.value,
                                   CT = cancerType,
                                   SC = seqCenter,
                                   ST = sampleType)
    
  }

  cat('\n')
  
  fisherListDf <- do.call(rbind,fisherList) %>%
    rstatix::adjust_pvalue(p.col = "p",method = "BH") %>%
    rstatix::add_significance(p.col = "p.adj") %>%
    mutate(abbrev = cancerAbbrevs[CT,"abbrev"])
  
  kendallListDf <- do.call(rbind,kendallList) %>%
    rstatix::adjust_pvalue(p.col = "p",method = "BH") %>%
    rstatix::add_significance(p.col = "p.adj") %>%
    mutate(abbrev = cancerAbbrevs[CT,"abbrev"])
  
  fisherKendallCombinedDf <- fisherListDf %>%
    rename(p.fisher = p, p.adj.fisher = p.adj,
           p.adj.signif.fisher = p.adj.signif) %>%
    left_join(kendallListDf[,c("tau","p","p.adj","p.adj.signif","abbrev")],
              by = "abbrev") %>%
    rename(p.kendall = p, p.adj.kendall = p.adj,
           p.adj.signif.kendall = p.adj.signif)
  
  p.adj.signif.fisher.kendall <- c(fisherKendallCombinedDf$p.adj.signif.fisher,
                                   fisherKendallCombinedDf$p.adj.signif.kendall)
  p.adj.signif.kendall.fisher <- c(fisherKendallCombinedDf$p.adj.signif.kendall,
                                   fisherKendallCombinedDf$p.adj.signif.fisher)

  if(kendallOnlyFlag){

    # Plot Kendall tau
    kendallListDf %>%
      ggplot(aes(x = reorder(abbrev, -tau,median), y = tau)) +
      geom_bar(stat="identity", color="black", position=position_dodge()) +
      theme_pubr() +
      rotate_x_text(30) +
      labs(x = "TCGA Cancer Type", 
           y = "Kendall tau",
           title = paste(seqCenter,sampleType,sep = " | ")) +
      geom_text(aes(label = p.adj.signif, y = tau), vjust = -0.4) +
      ylim(c(0,1.1*max(kendallListDf$tau))) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") -> barPlotKendallTau

    print(barPlotKendallTau)
    fileNameTau <- paste0("Figures/kendall_tau_barplot_tcga_",fileNameString,"_",seqCenter,"_",
                       gsub('([[:punct:]])|\\s+','',sampleType),
                       ".jpeg")
    ggsave(filename = fileNameTau,
           plot = barPlotKendallTau,
           dpi = "retina", units = "in", height = 3.5, width = plotWidth)

    # Plot Fisher and Kendall p-values
    fisherKendallCombinedDf %>%
      mutate(log.p.adj.fisher = -log10(p.adj.fisher),
             log.p.adj.kendall = -log10(p.adj.kendall)) %>%
      select(abbrev, log.p.adj.kendall, p.adj.signif.kendall) %>%
      reshape2::melt(id.vars = c("abbrev","p.adj.signif.kendall")) %>%
      mutate(variable = factor(gsub("log\\.p\\.adj\\.","",variable), levels = c("kendall"))) %>%
      mutate(p.adj.signif.combined = p.adj.signif.kendall) %>%
      # mutate(p.adj.signif.combined = `if`(head(as.character(variable),1)=="fisher",
      #                                p.adj.signif.fisher.kendall,
      #                                p.adj.signif.kendall.fisher)) %>%
      ggplot(aes(x = reorder(abbrev,-value,median), y = value, fill=variable)) +
      geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
      geom_text(aes(label = p.adj.signif.combined, y = value),
                position = position_dodge(0.9), 
                hjust = -0.2,
                vjust = 0.8,
                size = 3,
                angle = 90) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
      ylim(c(0,1.1*max( c(-log10(fisherKendallCombinedDf$p.adj.fisher),
                          -log10(fisherKendallCombinedDf$p.adj.kendall)) ))) +
      scale_fill_nejm(name = "Test type") +
      theme_pubr() +
      rotate_x_text(30) +
      labs(x = "TCGA Cancer Type", 
           y = "-Log(p-adjust)",
           title = paste(seqCenter,sampleType,sep = " | ")) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") -> pvalPlot

    print(pvalPlot)
    fileNamePval <- paste0("Figures/pvalue_Kendall_barplot_tcga_",fileNameString,"_",seqCenter,"_",
                       gsub('([[:punct:]])|\\s+','',sampleType),
                       ".jpeg")
    ggsave(filename = fileNamePval,
           plot = pvalPlot,
           dpi = "retina", units = "in", height = 3.5, width = plotWidth)

    res <- list(fisherListDf=fisherListDf,
              pvalPlot=pvalPlot,
              kendallListDf=kendallListDf,
              fisherKendallCombinedDf=fisherKendallCombinedDf)


  } else{

    # Plot Fisher exact ORs
    fisherListDf %>%
      ggplot(aes(x = reorder(abbrev,-OR,median), y = OR)) +
      geom_bar(stat="identity", color="black", position=position_dodge()) +
      geom_errorbar(aes(ymin=OR_low, ymax=OR_high), width=.2,
                    position=position_dodge(.9)) +
      theme_pubr() +
      rotate_x_text(30) +
      labs(x = "TCGA Cancer Type", 
           y = "Odds ratio feature enrichment",
           title = paste(seqCenter,sampleType,sep = " | ")) +
      geom_text(aes(label = p.adj.signif, y = OR_high), vjust = -0.4) +
      ylim(c(0,1.1*max(fisherListDf$OR_high))) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") -> barPlotOR
    
    print(barPlotOR)
    fileNameOR <- paste0("Figures/odds_ratio_barplot_tcga_",fileNameString,"_",seqCenter,"_",
                       gsub('([[:punct:]])|\\s+','',sampleType),
                       ".jpeg")
    ggsave(filename = fileNameOR,
           plot = barPlotOR,
           dpi = "retina", units = "in", height = 3.5, width = plotWidth)
    
    # Plot Kendall tau
    kendallListDf %>%
      ggplot(aes(x = reorder(abbrev, -tau,median), y = tau)) +
      geom_bar(stat="identity", color="black", position=position_dodge()) +
      theme_pubr() +
      rotate_x_text(30) +
      labs(x = "TCGA Cancer Type", 
           y = "Kendall tau",
           title = paste(seqCenter,sampleType,sep = " | ")) +
      geom_text(aes(label = p.adj.signif, y = tau), vjust = -0.4) +
      ylim(c(0,1.1*max(kendallListDf$tau))) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") -> barPlotKendallTau
    
    print(barPlotKendallTau)
    fileNameTau <- paste0("Figures/kendall_tau_barplot_tcga_",fileNameString,"_",seqCenter,"_",
                       gsub('([[:punct:]])|\\s+','',sampleType),
                       ".jpeg")
    ggsave(filename = fileNameTau,
           plot = barPlotKendallTau,
           dpi = "retina", units = "in", height = 3.5, width = plotWidth)
    
    # Plot Fisher and Kendall p-values
    fisherKendallCombinedDf %>%
      mutate(log.p.adj.fisher = -log10(p.adj.fisher),
             log.p.adj.kendall = -log10(p.adj.kendall)) %>%
      select(abbrev, log.p.adj.fisher, log.p.adj.kendall,
             p.adj.signif.fisher, p.adj.signif.kendall) %>%
      reshape2::melt(id.vars = c("abbrev","p.adj.signif.fisher","p.adj.signif.kendall")) %>%
      mutate(variable = factor(gsub("log\\.p\\.adj\\.","",variable), levels = c("kendall","fisher"))) %>%
      mutate(p.adj.signif.combined = `if`(head(as.character(variable),1)=="fisher",
                                     p.adj.signif.fisher.kendall,
                                     p.adj.signif.kendall.fisher)) %>%
      ggplot(aes(x = reorder(abbrev,-value,median), y = value, fill=variable)) +
      geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
      geom_text(aes(label = p.adj.signif.combined, y = value),
                position = position_dodge(0.9), 
                hjust = -0.2,
                vjust = 0.8,
                size = 3,
                angle = 90) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
      ylim(c(0,1.1*max( c(-log10(fisherKendallCombinedDf$p.adj.fisher),
                          -log10(fisherKendallCombinedDf$p.adj.kendall)) ))) +
      scale_fill_nejm(name = "Test type") +
      theme_pubr() +
      rotate_x_text(30) +
      labs(x = "TCGA Cancer Type", 
           y = "-Log(p-adjust)",
           title = paste(seqCenter,sampleType,sep = " | ")) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") -> pvalPlot
    
    print(pvalPlot)
    fileNamePval <- paste0("Figures/pvalue_Fisher_Kendall_barplot_tcga_",fileNameString,"_",seqCenter,"_",
                       gsub('([[:punct:]])|\\s+','',sampleType),
                       ".jpeg")
    ggsave(filename = fileNamePval,
           plot = pvalPlot,
           dpi = "retina", units = "in", height = 3.5, width = plotWidth)

    res <- list(fisherListDf=fisherListDf,
              barPlotOR=barPlotOR,
              pvalPlot=pvalPlot,
              kendallListDf=kendallListDf,
              fisherKendallCombinedDf=fisherKendallCombinedDf)

  }
  
  return(res)
}


runSeqCenterFXN <- function(inputDataDf,
                            intFlagInput=TRUE,
                            fileString = "vsnmVsRaw",
                            factorCQInput = FALSE){
  
  source("00-functions.R") # for barplotSummaryPerf() and barplotPerf() functions
  
  #----------------All----------------#
  barplotSummaryPerf(inputData = inputDataDf,
                     seqCenterAbbrev="All",
                     sampleTypeInput = "Blood Derived Normal",
                     intFlag = intFlagInput,
                     plotWidthSingle = 8,
                     plotWidthCombined = 8,
                     fileNameString = fileString,
                     factorCQ = factorCQInput)
  
  barplotSummaryPerf(inputData = inputDataDf,
                     seqCenterAbbrev="All",
                     sampleTypeInput = "Primary Tumor",
                     intFlag = intFlagInput,
                     plotWidthSingle = 12,
                     plotWidthCombined = 12,
                     fileNameString = fileString,
                     factorCQ = factorCQInput)
  
  barplotSummaryPerf(inputData = inputDataDf,
                     seqCenterAbbrev="All",
                     sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                     intFlag = intFlagInput,
                     plotWidthSingle = 5,
                     plotWidthCombined = 5,
                     fileNameString = fileString,
                     factorCQ = factorCQInput)
  
  #----------------WGS----------------#
  # HMS
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="HMS",
              sampleTypeInput = "Blood Derived Normal",
              intFlag = FALSE,
              plotWidthSingle = 6,
              plotWidthCombined = 8,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="HMS",
              sampleTypeInput = "Primary Tumor",
              intFlag = intFlagInput,
              plotWidthSingle = 6,
              plotWidthCombined = 8,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="HMS",
              sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
              intFlag = intFlagInput,
              plotWidthSingle = 3,
              plotWidthCombined = 5,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  # BCM
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="BCM",
              sampleTypeInput = "Blood Derived Normal",
              intFlag = intFlagInput,
              plotWidthSingle = 4,
              plotWidthCombined = 6,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="BCM",
              sampleTypeInput = "Primary Tumor",
              intFlag = intFlagInput,
              plotWidthSingle = 4,
              plotWidthCombined = 6,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="BCM",
              sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
              intFlag = intFlagInput,
              plotWidthSingle = 3,
              plotWidthCombined = 5,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  # MDA
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="MDA",
              sampleTypeInput = "Blood Derived Normal",
              intFlag = intFlagInput,
              plotWidthSingle = 4,
              plotWidthCombined = 6,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="MDA",
              sampleTypeInput = "Primary Tumor",
              intFlag = intFlagInput,
              plotWidthSingle = 4,
              plotWidthCombined = 6,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  # WashU
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="WashU",
              sampleTypeInput = "Blood Derived Normal",
              intFlag = intFlagInput,
              plotWidthSingle = 3,
              plotWidthCombined = 5,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="WashU",
              sampleTypeInput = "Primary Tumor",
              intFlag = intFlagInput,
              plotWidthSingle = 3,
              plotWidthCombined = 5,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  # Broad_WGS
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="Broad_WGS",
              sampleTypeInput = "Blood Derived Normal",
              intFlag = intFlagInput,
              plotWidthSingle = 4,
              plotWidthCombined = 6,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="Broad_WGS",
              sampleTypeInput = "Primary Tumor",
              intFlag = intFlagInput,
              plotWidthSingle = 6,
              plotWidthCombined = 8,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="Broad_WGS",
              sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
              intFlag = intFlagInput,
              plotWidthSingle = 3,
              plotWidthCombined = 5,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  #----------------RNA----------------#
  
  # UNC
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="UNC",
              sampleTypeInput = "Primary Tumor",
              intFlag = intFlagInput,
              plotWidthSingle = 8,
              plotWidthCombined = 18,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="UNC",
              sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
              intFlag = intFlagInput,
              plotWidthSingle = 6,
              plotWidthCombined = 10,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  # CMS
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="CMS",
              sampleTypeInput = "Primary Tumor",
              intFlag = intFlagInput,
              plotWidthSingle = 3,
              plotWidthCombined = 5,
              fileNameString = fileString,
              factorCQ = factorCQInput)
  
  barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="CMS",
              sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
              intFlag = intFlagInput,
              plotWidthSingle = 2,
              plotWidthCombined = 3,
              fileNameString = fileString,
              factorCQ = factorCQInput)
}






runSeqCenterEnrichmentFxn2 <- function(pathInput = "Supporting_scripts/S10-ML-10k-tcga-raw-vsnm-subsets/features__",
                            totalFeaturesInput = colnames(snmDataSampleTypeWithExpStrategy),
                            datasetRawInput = "vbDataBarnDFReconciledQC_Int_HiSeq",
                            datasetVSNMInput = "snmDataSampleTypeWithExpStrategy_HiSeq",
                            fileNameStringInput = "vsnmVsRaw",
                            kendallOnlyFlagInput = FALSE,
                            wgsRun = TRUE,
                            rnaRun = TRUE){
  
  source("00-functions.R") # for enrichmentFxn2() function

  if(wgsRun){

    #----------------WGS----------------#
    # HMS
    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="HMS",
                  sampleType = "Blood Derived Normal",
                  plotWidth = 5)

    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="HMS",
                  sampleType = "Primary Tumor",
                  plotWidth = 5)

    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="HMS",
                  sampleType = "Primary Tumor vs Solid Tissue Normal",
                  plotWidth = 2)

    # BCM
    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="BCM",
                  sampleType = "Blood Derived Normal",
                  plotWidth = 4)

    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="BCM",
                  sampleType = "Primary Tumor",
                  plotWidth = 4)

    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="BCM",
                  sampleType = "Primary Tumor vs Solid Tissue Normal",
                  plotWidth = 2)

    # MDA
    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="MDA",
                  sampleType = "Blood Derived Normal",
                  plotWidth = 4)

    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="MDA",
                  sampleType = "Primary Tumor",
                  plotWidth = 4)

    # WashU
    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="WashU",
                  sampleType = "Blood Derived Normal",
                  plotWidth = 3)

    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="WashU",
                  sampleType = "Primary Tumor",
                  plotWidth = 3)

    # Broad_WGS
    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="Broad_WGS",
                  sampleType = "Blood Derived Normal",
                  plotWidth = 4)

    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="Broad_WGS",
                  sampleType = "Primary Tumor vs Solid Tissue Normal",
                  plotWidth = 2)

  } else{
    print("WGS **not** run")
  }
  
  if(rnaRun){

    #----------------RNA----------------#
    # UNC
    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="UNC",
                  sampleType = "Primary Tumor",
                  plotWidth = 10)

    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="UNC",
                  sampleType = "Primary Tumor vs Solid Tissue Normal",
                  plotWidth = 6)

    # CMS
    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="CMS",
                  sampleType = "Primary Tumor",
                  plotWidth = 3)

    enrichmentFxn2(myPath = pathInput,
                  totalFeatures = totalFeaturesInput,
                  dataSetRaw = datasetRawInput,
                  dataSetVSNM = datasetVSNMInput,
                  fileNameString = fileNameStringInput,
                  kendallOnlyFlag = kendallOnlyFlagInput,
                  seqCenter="CMS",
                  sampleType = "Primary Tumor vs Solid Tissue Normal",
                  plotWidth = 2)

  } else{
    print("RNA **not** run")

  }
  
}


runAncomBC_1VsAll_OGUs <- function(metaString = "metaRSFinal9090_Nonzero_HiSeq",
                                   countString = "rs210PanFinal9090_Nonzero_HiSeq",
                                   dataString = "rs210Pan_Filt9090",
                                   taxTable = taxRS210_ff_combTaxaSpeciesZebraWithViruses,
                                   makeTaxTable = FALSE,
                                   qvalCutoff = 0.05,
                                   showTopX = 3,
                                   showTopXFlag=FALSE,
                                   ancombcLibCut = 0,
                                   sampleTypes = c("Primary Tumor","Blood Derived Normal"),
                                   SeqCenters = c("HMS","BCM","MDA","WashU","Broad_WGS"),
                                   taxaPlotLabel = "genus"){
  
  require(ANCOMBC)
  require(ggrepel)
  require(ggsci)
  
  if(makeTaxTable){
    print("Making synthetic taxa table")
  } else{
    # Note that taxa and count tables have to be a matrix
    taxTableFormatted <- taxTable %>%
      select(-UNITN,-UHGG,-PATH,-WIS,-cov_All) %>%
      column_to_rownames("OGU")
      
    psTaxTable <- taxTableFormatted %>% 
      rename(Domain=domain,Phylum=phylum,Class=class,Order=order,
             Family=family,Genus=genus,Species=species) %>%
      as.matrix()
  }
  
  
  plotPrefix <- paste0("Figures/ancombc_",dataString,"/")
  # Create folder for plots if doesn't exist
  plotFolder <- paste0("Figures/ancombc_",dataString)
  if(!( dir.exists( file.path(plotFolder)))){
    dir.create(file.path(plotFolder))
  }
  
  for(zz in seq_along(sampleTypes)){
    
    sampleType <- sampleTypes[zz]
    sampleTypeFormatted <- ifelse(sampleType=="Primary Tumor","PT","BDN")
    
    for(jj in seq_along(SeqCenters)){
      
      SeqCenterFormatted <- SeqCenters[jj]
      SeqCenter <- SeqCenterFormatted
      # SeqCenterFormatted <- gsub('([[:punct:]])|\\s+','',SeqCenter)
      metaDataFilt <- eval(as.name(paste0(metaString,"_",SeqCenterFormatted))) %>%
        mutate(investigation_formatted = gsub("^TCGA-","",investigation)) %>%
        filter(sample_type == sampleType) %>%
        droplevels()
      countDataFilt <- eval(as.name(paste0(countString,"_",SeqCenterFormatted)))[rownames(metaDataFilt),]
      
      cancerTypes <- as.character(unique(metaDataFilt$disease_type))
      
      for(ii in seq_along(cancerTypes)){
        Dz <- cancerTypes[ii]
        DzFormatted <- gsub('([[:punct:]])|\\s+','',Dz)
        
        metaDataFilt$predY <- factor(ifelse(metaDataFilt$disease_type == Dz,
                                            yes = Dz, no = "Other"),
                                     levels = c("Other",Dz))
        countDataFilt_Dz <- countDataFilt[which(metaDataFilt$disease_type == Dz),]
        countDataFilt_Other <- countDataFilt[which(metaDataFilt$disease_type != Dz),]
        investigationText <- metaDataFilt$investigation_formatted[which(metaDataFilt$disease_type == Dz)[1]]
        
        # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
        # in case a seq center only has 1 cancer type
        if(length(table(metaDataFilt$predY)) < 2){next}
        
        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 10 SAMPLES IN EITHER CLASS
        if(any(table(metaDataFilt$predY) < 10)){next}
        
        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF ancombcLibCut READS/SAMPLE IN EITHER CLASS
        if(all(rowSums(countDataFilt_Dz) < ancombcLibCut)){next}
        if(all(rowSums(countDataFilt_Other) < ancombcLibCut)){next}
        
        # If sufficient samples, then:
        print(SeqCenter)
        print(Dz)
        print(sprintf("Number of samples (Dz|Other): %d | %d", 
                      unname(table(metaDataFilt$predY)[2]),
                      unname(table(metaDataFilt$predY)[1])))

        if(makeTaxTable){
          taxTableFormatted <- data.frame(domain = colnames(countDataFilt),
                                         phylum = colnames(countDataFilt),
                                         class = colnames(countDataFilt),
                                         order = colnames(countDataFilt),
                                         family = colnames(countDataFilt),
                                         genus = colnames(countDataFilt),
                                         species = colnames(countDataFilt),
                                         row.names = colnames(countDataFilt))

          psTaxTable <- taxTableFormatted %>% 
            rename(Domain=domain,Phylum=phylum,Class=class,Order=order,
                   Family=family,Genus=genus,Species=species) %>%
            as.matrix()
        }
        
        ps_1vsAll <- phyloseq(otu_table(countDataFilt, taxa_are_rows = FALSE), 
                                               tax_table(as.matrix(psTaxTable)), 
                                               sample_data(metaDataFilt))
        # print(ps_1vsAll)
        print(sprintf("Read count distribution:"))
        print(summary(rowSums(otu_table(ps_1vsAll)))) # Sample read distribution
        print(sprintf("Now running ANCOM-BC..."))
        ancombc_OGU_1vsAll_X <- ancombc(phyloseq = ps_1vsAll, 
                                       formula = "predY", p_adj_method = "BH", zero_cut = 0.999, 
                                       lib_cut = ancombcLibCut, 
                                       # group = "sample_type", struc_zero = FALSE, neg_lb = FALSE,
                                       tol = 1e-5, max_iter = 100, conserve = FALSE, alpha = 0.05, global = FALSE)
        
        print(sprintf("Plotting and saving data..."))
        ancom_res_df_OGU_1vsAll_X <- data.frame(
          beta = unlist(ancombc_OGU_1vsAll_X$res$beta),
          se = unlist(ancombc_OGU_1vsAll_X$res$se),
          W = unlist(ancombc_OGU_1vsAll_X$res$W),
          p_val = unlist(ancombc_OGU_1vsAll_X$res$p_val),
          q_val = unlist(ancombc_OGU_1vsAll_X$res$q_val),
          diff_abn = ifelse(unlist(ancombc_OGU_1vsAll_X$res$q_val)<=qvalCutoff, yes = TRUE, no = FALSE),
          genus = taxTableFormatted[row.names(ancombc_OGU_1vsAll_X$res$beta),"genus"],
          species = taxTableFormatted[row.names(ancombc_OGU_1vsAll_X$res$beta),"species"],
          OGUs = row.names(ancombc_OGU_1vsAll_X$res$beta),
          row.names = row.names(ancombc_OGU_1vsAll_X$res$beta))
        ancom_res_df_OGU_1vsAll_X$diff_name_flag <- ancom_res_df_OGU_1vsAll_X$diff_abn + 0
        ancom_res_df_OGU_1vsAll_X_sorted <- ancom_res_df_OGU_1vsAll_X[order(ancom_res_df_OGU_1vsAll_X$q_val),]
        
        ancom_res_df_OGU_1vsAll_X_sorted$diff_label <- ""
        ancom_res_df_OGU_1vsAll_X_sorted$diff_label[1:showTopX] <- ifelse(ancom_res_df_OGU_1vsAll_X_sorted$diff_abn, 
                                                                         yes = paste0(ancom_res_df_OGU_1vsAll_X_sorted$OGUs,"\n(",ancom_res_df_OGU_1vsAll_X_sorted[,taxaPlotLabel],")"),
                                                                         no = "")[1:showTopX]
        
        print(sprintf("Number of differentially abundant OGUs: %d", 
                      sum(ancom_res_df_OGU_1vsAll_X_sorted$diff_name_flag)))
        cat("\n")
        
        ## Set up sub-axis CT (right) vs. Others (left) labels
        library(grid)
        text_high <- textGrob(investigationText, gp=gpar(fontsize=11, fontface="bold.italic"))
        text_low <- textGrob("Other", gp=gpar(fontsize=11, fontface="bold.italic"))
        text_high_xpos <- ifelse(max(ancom_res_df_OGU_1vsAll_X_sorted$beta)<=0.2, yes = 0.2, no = max(ancom_res_df_OGU_1vsAll_X_sorted$beta))
        text_low_xpos <- ifelse(min(ancom_res_df_OGU_1vsAll_X_sorted$beta)>=-0.2, yes = -0.2, no = min(ancom_res_df_OGU_1vsAll_X_sorted$beta))
        yval <- -log10(ancom_res_df_OGU_1vsAll_X_sorted$p_val)
        text_ypos <- ifelse(1.1*max(yval)<=-log10(0.05), 
                            yes = 1.1*-log10(0.05), no = 1.1*max(yval))
        text_ypos <- ifelse(is.finite(text_ypos), yes = text_ypos, no = 1.1*max(yval[is.finite(yval)]))
        
        if(showTopXFlag){
          ancom_res_df_OGU_1vsAll_X_sorted %>%
          ggplot(aes(x = beta, y = -log10(p_val), color = diff_abn, label = diff_label)) + geom_point(size = 2) +
          theme_bw() + geom_hline(yintercept=-log10(0.05), col="black", linetype='dashed') +
          geom_vline(xintercept=c(-0.2, 0.2), col="black", linetype='dashed') + 
          theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + coord_fixed(ylim = c(0, NA), clip = "off") +
          scale_color_aaas() +
          labs(x = "log-fold change\n(ANCOM-BC beta)", y = "-Log10(p-value)", 
               color = paste0("Differentially\nabundant\nOGUs\n(q<=",qvalCutoff,")"), 
               title = paste(paste0(SeqCenter," (",sampleType,")"), Dz, "", sep = "\n")) +
          geom_label_repel(force = 20, size = 2, box.padding = 2, point.padding = 1e-06, label.size = 0.2, show.legend = FALSE, color = "black", max.overlaps = 10) +
          theme(plot.margin = unit(c(1,1,2,1), "lines")) +
          annotation_custom(text_high,xmin=text_high_xpos,xmax=text_high_xpos,ymin=text_ypos, ymax=text_ypos) + 
          annotation_custom(text_low,xmin=text_low_xpos,xmax=text_low_xpos,ymin=text_ypos, ymax=text_ypos) -> p
          } else{
            ancom_res_df_OGU_1vsAll_X_sorted %>%
          ggplot(aes(x = beta, y = -log10(p_val), color = diff_abn)) + geom_point(size = 2) +
          theme_bw() + geom_hline(yintercept=-log10(0.05), col="black", linetype='dashed') +
          geom_vline(xintercept=c(-0.2, 0.2), col="black", linetype='dashed') + 
          theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + coord_fixed(ylim = c(0, NA), clip = "off") +
          scale_color_aaas() +
          labs(x = "log-fold change\n(ANCOM-BC beta)", y = "-Log10(p-value)", 
               color = paste0("Differentially\nabundant\nOGUs\n(q<=",qvalCutoff,")"), 
               title = paste(paste0(SeqCenter," (",sampleType,")"), Dz, "", sep = "\n")) +
          # geom_label_repel(force = 20, size = 2, box.padding = 2, point.padding = 1e-06, label.size = 0.2, show.legend = FALSE, color = "black", max.overlaps = 10) +
          theme(plot.margin = unit(c(1,1,2,1), "lines")) +
          annotation_custom(text_high,xmin=text_high_xpos,xmax=text_high_xpos,ymin=text_ypos, ymax=text_ypos) + 
          annotation_custom(text_low,xmin=text_low_xpos,xmax=text_low_xpos,ymin=text_ypos, ymax=text_ypos) -> p
          }
        
        ggsave(plot=p, filename = paste0(plotPrefix, sampleTypeFormatted, "_", 
                                         SeqCenterFormatted,"_",DzFormatted,".jpeg"), 
               dpi = "retina", width = 8, height = 6, units = "in")
        # Write data to file
        ancom_res_df_OGU_1vsAll_X_sorted %>% 
          write.csv(file = paste0(plotPrefix, "data_",  sampleTypeFormatted, "_", 
                                  SeqCenterFormatted,"_",DzFormatted,".csv"))
      }
    }
  }
}


ancomReplot <- function(seqCenter = "BCM",
                        sampleType = "Primary Tumor",
                        dataString = "kuT2T_BIO",
                        numticks = 3,
                        fontSize = 8,
                        pointSize = 1,
                        plotWidth = 14,
                        plotHeight = 2,
                        nRowFlag = FALSE,
                        nRow = 4){
  require(ggpubr)
  require(ggsci)
  require(gridExtra)
  require(scales)

  # Load abbreviations
  abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", 
                                          stringsAsFactors = FALSE, row.names = 1)
  abbreviationsTCGA_AllcancerMod <- abbreviationsTCGA_Allcancer %>%
    rownames_to_column("CT") %>%
    mutate(CT = gsub('([[:punct:]])|\\s+','',CT)) %>%
    column_to_rownames("CT")
  
  ancomPlots <- list()
  sampleTypeFormatted <- ifelse(sampleType=="Primary Tumor","PT","BDN")
  fileList <- list.files(paste0("Figures/ancombc_",dataString,"/"))
  fileListFilt <- fileList[grepl(paste0("^data_",sampleTypeFormatted,"_",seqCenter),fileList)]
  CTs <- gsub("\\.csv$","",gsub(paste0("^data_",sampleTypeFormatted,"_",seqCenter,"_"),"",fileListFilt))
  
  for(ii in 1:length(CTs)){
    CT <- CTs[ii]
    abbrev <- abbreviationsTCGA_AllcancerMod[CT,"abbrev"]
    plotPrefix <- paste0("Figures/ancombc_",dataString,"/")
    sampleTypeFormatted <- ifelse(sampleType=="Primary Tumor","PT","BDN")
    ancombcData <- read.csv(file = paste0(plotPrefix, "data_",  sampleTypeFormatted, "_", 
                                          seqCenter,"_",CT,".csv"))
    ancombcData %>%
      ggplot(aes(x = beta, y = -log10(p_val), color = diff_abn)) + 
      geom_point(size = pointSize) +
      theme_pubr() + 
      geom_hline(yintercept=-log10(0.05), col="black", linetype='solid', linewidth=0.2) +
      geom_vline(xintercept=c(-0.2, 0.2), col="black", linetype='solid', linewidth=0.2) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_aaas() +
      labs(x="Log-fold change",y="-Log10(p-value)",title = abbrev) +
      # rotate_x_text(30) +
      scale_x_continuous(breaks = trans_breaks(identity, identity, n = numticks)) +
      theme(legend.position = "none",
            plot.title = element_text(face="bold"),
            text = element_text(size=fontSize)) -> ancomPlots[[ii]]
  }
  # Create folder for plots if doesn't exist
  outputPrefix <- paste0("Figures/ancombcMerged_",dataString,"/")
  outputFolder <- paste0("Figures/ancombcMerged_",dataString)
  if(!( dir.exists( file.path(outputFolder)))){
    dir.create(file.path(outputFolder))
  }
  if(nRowFlag){
      pGrid <- do.call("grid.arrange", c(ancomPlots, nrow=nRow,clip=TRUE))
    } else{
      pGrid <- do.call("grid.arrange", c(ancomPlots, ncol=length(CTs),clip=TRUE))
    }
  ggsave(plot = pGrid,
         filename = paste0(outputPrefix,seqCenter,"_",sampleTypeFormatted,".jpeg"),
         dpi = 900, units = "in", width = plotWidth, height = plotHeight)
}