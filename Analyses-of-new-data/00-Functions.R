#-----------------------------------------------------------------------------
# 00-Functions.R
# Copyright (c) 2021--, Greg Poore
# Purpose: List main functions for R scripts
#-----------------------------------------------------------------------------

barplotPerf <- function(inputData=mlPerfAll10k_Allcancer_Raw_VSNM_SeqCenter, 
                        sampleTypeInput, 
                        seqCenterAbbrev, 
                        plotWidthSingle = 6,
                        plotWidthCombined = 10,
                        prefixRaw = "rs210PanFinal7575_Nonzero_HiSeq_",
                        intFlag = FALSE,
                        ciWidth = 0.95,
                        singlePlotHeight = 3.5,
                        combPlotHeight = 3.5,
                        fileNameString = "rs210Pan_WIS",
                        errorbarPlotFlag = FALSE,
                        factorCQ = FALSE,
                        factorHD = FALSE # HD for host depletion
                        ){

  plotPrefix <- paste0("Figures/mlBarplots_",fileNameString,"/")
  # Create folder for plots if doesn't exist
  plotFolder <- paste0("Figures/mlBarplots_",fileNameString)
  if(!( dir.exists( file.path(plotFolder)))){
    dir.create(file.path(plotFolder))
  }

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
    mutate(datasetName = gsub(prefixRaw,"",datasetName)) %>%
    mutate(datasetName = gsub("vsnm.+","VSNM",datasetName)) %>%
    mutate(datasetName = gsub("WGS_|RNA_","",datasetName))

    if(factorCQ){
      inputDataFilt <- inputDataFilt %>%
      mutate(datasetName = factor(case_when(
        grepl("CQ",datasetName) ~ "ConQuR",
        grepl("VSNM",datasetName) ~ "VSNM",
        grepl(seqCenterAbbrev,datasetName) ~ "Raw"
        ), levels = c("Raw","VSNM","ConQuR")))
    }

    if(factorHD){ # HD for host depletion
      inputDataFilt <- inputDataFilt %>%
      mutate(datasetName = factor(case_when(
        grepl("HG38",datasetName) ~ "hg38",
        grepl("T2T",datasetName) ~ "T2T",
        grepl("Pan",datasetName) ~ "Pangenome"
        ), levels = c("hg38","T2T","Pangenome")))
    }
  
  inputDataFilt %>%
    filter(variable == "AUROC") %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),
                  lty="11",position = position_dodge(0.9)) +
    scale_fill_manual(values = c("#ADB6B6FF","#925E9FFF","#E18727FF"), name = "Data type") +
    # scale_fill_nejm(name = "Data type") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUROC",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUROC

  fileNameAUROC <- paste0(plotPrefix,"auroc_",seqCenterAbbrev,"_",
                          gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                          ".jpeg")
  
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
    geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),
                  lty="11",position = position_dodge(0.9)) +
    scale_fill_manual(values = c("#ADB6B6FF","#925E9FFF","#E18727FF"), name = "Data type") +
    # scale_fill_nejm(name = "Data type") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUPR",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUPR

  fileNameAUPR <- paste0(plotPrefix,"aupr_",seqCenterAbbrev,"_",
                           gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                           ".jpeg")
  
  ggsave(filename = fileNameAUPR,
         plot = barPlotAUPR,
         dpi = "retina", units = "in", height = singlePlotHeight, width = plotWidthSingle)
  
  inputDataFilt %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),
                  lty="11",position = position_dodge(0.9)) +
    geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),
                  lty="11",position = position_dodge(0.9)) +
    facet_wrap(vars(variable)) +
    scale_fill_manual(values = c("#ADB6B6FF","#925E9FFF","#E18727FF"), name = "Data type") +
    # scale_fill_nejm(name = "Data type") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "Area Under Curve",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotCombined
  
  print(barPlotCombined)

  fileName <- paste0(plotPrefix,"combined_",seqCenterAbbrev,"_",
                       gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                       ".jpeg")
  
  ggsave(filename = fileName,
         plot = barPlotCombined,
         dpi = "retina", units = "in", height = combPlotHeight, width = plotWidthCombined)

  if(errorbarPlotFlag){
    inputDataFilt %>%
    ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),
                  width=0.4,size=0.6,position = position_dodge(0.9)) +
    geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="11",position = position_dodge(0.9)) +
    geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="11",position = position_dodge(0.9)) +
    geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    ggtitle(paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    rotate_x_text(30) + 
    scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") -> ebPlot
  
  print(ebPlot)

  fileNameEB <- paste0(plotPrefix,"ebPlot_",seqCenterAbbrev,"_",
                       gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                       ".jpeg")
  
  ggsave(filename = fileNameEB,
         plot = ebPlot,
         dpi = "retina", units = "in", height = singlePlotHeight, width = plotWidthSingle)

  }

  # return(inputDataFilt)
}

# Make one summary barplot 
barplotSummaryPerf <- function(inputData, 
                        sampleTypeInput, 
                        prefixRaw = "rs210PanFinal7575_Nonzero_HiSeq_",
                        seqCenterAbbrev="All", 
                        plotWidthSingle = 6,
                        plotWidthCombined = 10,
                        intFlag = FALSE,
                        ciWidth = 0.95,
                        fileNameString = "rs210Pan_WIS",
                        factorCQ = FALSE,
                        factorHD = FALSE # HD for host depletion
                        ){

  plotPrefix <- paste0("Figures/mlBarplots_",fileNameString,"/")
  # Create folder for plots if doesn't exist
  plotFolder <- paste0("Figures/mlBarplots_",fileNameString)
  if(!( dir.exists( file.path(plotFolder)))){
    dir.create(file.path(plotFolder))
  }

  inputDataFilt <- inputData %>%
    filter(sampleType == sampleTypeInput) %>%
    mutate(metadataName = "All") %>%
    select(AUROC, AUPR, abbrev,diseaseType,
           sampleType,datasetName,metadataName) %>%
    mutate(datasetName = gsub(prefixRaw,"",datasetName)) %>%
    mutate(datasetName = gsub("vsnm.+","VSNM",datasetName)) %>%
    mutate(datasetName = gsub("WGS_|RNA_","",datasetName)) %>%
    # mutate(datasetName = gsub("tcgaGenusKrakenAllFiltWIS_HiSeq_WGS_|tcgaGenusKrakenAllFiltWIS_HiSeq_RNA_","",datasetName)) %>%
    mutate(datasetName = case_when(
      datasetName == "Raw" ~ "Raw",
      grepl("CQ",datasetName) ~ "ConQuR",
      grepl("VSNM",datasetName) ~ "VSNM",
      grepl("HG38",datasetName) ~ "hg38",
      grepl("T2T",datasetName) ~ "T2T",
      grepl("Pan",datasetName) ~ "Pangenome",
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

  if(factorHD){ # HD for host depletion
    inputDataFilt <- inputDataFilt %>%
        mutate(datasetName = factor(datasetName, levels = c("hg38","T2T","Pangenome")))
    }

  # plotColors <- c("#ADB6B6FF","#925E9FFF","#E18727FF") # original
  plotColors <- c("#ADB6B6FF","#E18727FF","#925E9FFF")
  
  inputDataFilt %>%
    filter(variable == "AUROC") %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_manual(values = plotColors, name = "Data type") +
    # scale_fill_nejm(name = "Data type") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUROC",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUROC

  fileNameAUROC <- paste0(plotPrefix,"auroc_",seqCenterAbbrev,"_",
                          gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                          ".jpeg")
  
  ggsave(filename = fileNameAUROC,
         plot = barPlotAUROC,
         dpi = "retina", units = "in", height = 3.5, width = plotWidthSingle)
  
  inputDataFilt %>%
    filter(variable == "AUPR") %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_manual(values = plotColors, name = "Data type") +
    # scale_fill_nejm(name = "Data type") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "AUPR",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotAUPR

  fileNameAUPR <- paste0(plotPrefix,"aupr_",seqCenterAbbrev,"_",
                         gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                         ".jpeg")
  
  ggsave(filename = fileNameAUPR,
         plot = barPlotAUPR,
         dpi = "retina", units = "in", height = 3.5, width = plotWidthSingle)
  
  inputDataFilt %>%
    ggplot(aes(x = reorder(abbrev,value,median), y = value, fill = datasetName)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), 
                      ymax=ifelse(value+ci>1,1,value+ci)), width=.2,
                  position=position_dodge(.9)) +
    facet_wrap(vars(variable), nrow = 2) +
    scale_fill_manual(values = plotColors, name = "Data type") +
    # scale_fill_nejm(name = "Data type") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    theme_pubr() +
    rotate_x_text(30) +
    labs(x = "TCGA Cancer Type", 
         y = "Area Under Curve",
         title = paste(seqCenterAbbrev,sampleTypeInput,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotCombined
  
  print(barPlotCombined)

  fileName <- paste0(plotPrefix,"combined_",seqCenterAbbrev,"_",
                       gsub('([[:punct:]])|\\s+','',sampleTypeInput),
                       ".jpeg")
  
  ggsave(filename = fileName,
         plot = barPlotCombined,
         dpi = "retina", units = "in", height = 3.5, width = plotWidthCombined)
  
  # return(inputDataFilt)
}

#--------Expanded enrichment test function--------#

enrichmentFxn2 <- function(myPath = "Supporting_scripts/S24-RS210-Pan-RawVsCQ-Filt7575-seqcenter/features__",
                          totalFeatures = colnames(rs210PanFinal7575_Nonzero_HiSeq_WGS),
                          seqCenter = "HMS",
                          sampleType = "Primary Tumor",
                          dataSetRaw = "rs210PanFinal7575_Nonzero_HiSeq",
                          dataSetVSNM = "rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ",
                          cancerAbbrevs = abbreviationsTCGA_Allcancer,
                          plotWidth = 5,
                          fileNameString = "rs210Pan_Filt7575",
                          kendallOnlyFlag = FALSE,
                          showCM = FALSE){

  plotPrefix <- paste0("Figures/mlEnrich_",fileNameString,"/")
  # Create folder for plots if doesn't exist
  plotFolder <- paste0("Figures/mlEnrich_",fileNameString)
  if(!( dir.exists( file.path(plotFolder)))){
    dir.create(file.path(plotFolder))
  }

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
  
  # chi2List <- list()
  fisherList <- list()
  kendallList <- list()
  for(ii in 1:length(cancerTypes)){
    cancerType <- cancerTypes[ii]
    
    fileName <- paste0(cancerType," -- ",sampleType," -- Features.csv")
    
    featuresRaw <- read.csv(paste0(myPath,dataSetRaw_SeqCenter,"/",fileName)) %>%
      mutate(rank = row_number())
    featuresVSNM <- read.csv(paste0(myPath,dataSetVSNM_SeqCenter,"/",fileName)) %>%
      mutate(rank = row_number())
    
    # print(sprintf("Length raw list: %d",length(featuresRaw$Taxa)))
    # print(sprintf("Length VSNM list: %d",length(featuresVSNM$Taxa)))
    
    totalFeatures_NOT_VSNM <- totalFeatures[!(totalFeatures %in% featuresVSNM$Taxa)]
    totalFeatures_NOT_Raw <- totalFeatures[!(totalFeatures %in% featuresRaw$Taxa)]
    
    chi2df <- matrix(c(length(intersect(featuresRaw$Taxa,featuresVSNM$Taxa)),
                       sum( totalFeatures_NOT_VSNM %in% featuresRaw$Taxa),
                       sum( totalFeatures_NOT_Raw %in% featuresVSNM$Taxa),
                       sum(totalFeatures_NOT_VSNM %in% totalFeatures_NOT_Raw) ),
                     nrow = 2, ncol = 2)
    colnames(chi2df) <- c("Raw+","Raw-")
    rownames(chi2df) <- c("VSNM+","VSNM-")
    # # Chi2 test
    # chi2Test <- chisq.test(chi2df)
    # chi2List[[ii]] <- chi2List
    if(showCM){
      print(cancerType)
      print(chi2df)
    }
    
    # Fisher exact test
    fisherTest <- fisher.test(chi2df)
    # print(fisherTest$estimate)
    # print(fisherTest$null.value)
    fisherList[[ii]] <- data.frame(OR = ifelse(as.numeric(fisherTest$estimate)==Inf,
                                                yes = 1e6, # substitute high value for Inf OR
                                                no = as.numeric(fisherTest$estimate)),
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

  # cat('\n')
  
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
    fileNameTau <- paste0(plotPrefix,"kendall_tau_",seqCenter,"_",
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
    fileNamePval <- paste0(plotPrefix,"Figures/pvalue_Kendall_",seqCenter,"_",
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
      # geom_errorbar(aes(ymin=OR_low, ymax=OR_high), width=.2,
      #               position=position_dodge(.9)) +
      theme_pubr() +
      rotate_x_text(30) +
      labs(x = "TCGA Cancer Type", 
           y = "log10(odds ratio)",
           title = paste(seqCenter,sampleType,sep = " | ")) +
      # geom_text(aes(label = p.adj.signif, y = OR_high), vjust = -0.4) +
      geom_text(aes(label = p.adj.signif, y = OR), vjust = -0.4) +
      # ylim(c(0,1.1*max(fisherListDf$OR))) +
      scale_y_log10(limits = c(1,5*max(fisherListDf$OR))) +
      # ylim(c(0,1.1*max(fisherListDf$OR_high))) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") -> barPlotOR
    
    print(barPlotOR)
    fileNameOR <- paste0(plotPrefix,"OR_",seqCenter,"_",
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
    fileNameTau <- paste0(plotPrefix,"kendall_tau_",seqCenter,"_",
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
           y = "-log10(p-adjust)",
           title = paste(seqCenter,sampleType,sep = " | ")) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") -> pvalPlot
    
    print(pvalPlot)
    fileNamePval <- paste0(plotPrefix,"pvalue_Fisher_Kendall_",seqCenter,"_",
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

runBarplotPerfCQFXN <- function(inputDataDf,
                                intFlagInput=FALSE,
                                prefixRawInput = "rs210PanFinal7575_Nonzero_HiSeq_",
                                fileString = "rs210Pan_Filt7575",
                                errorbarPlotFlagInput = FALSE,
                                factorCQInput = FALSE,
                                factorHDInput = FALSE # HD for host depletion
                                ){

  source("00-Functions.R") # for barplotPerf() function
  
  #----------------BDN----------------#
  if(sum(grepl("CQ",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Blood Derived Normal")>0){
    barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="CQ",
              sampleTypeInput = "Blood Derived Normal",
              prefixRaw = prefixRawInput,
              intFlag = FALSE,
              plotWidthSingle = 8,
              plotWidthCombined = 12,
              errorbarPlotFlag = errorbarPlotFlagInput,
              fileNameString = fileString,
              factorCQ = factorCQInput,
              factorHD = factorHDInput)
    } else{
      print("BDN not run!")
    }

  #----------------PT----------------#
  if(sum(grepl("CQ",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor")>0){
    barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="CQ",
              sampleTypeInput = "Primary Tumor",
              prefixRaw = prefixRawInput,
              intFlag = FALSE,
              plotWidthSingle = 8,
              plotWidthCombined = 14,
              errorbarPlotFlag = errorbarPlotFlagInput,
              fileNameString = fileString,
              factorCQ = factorCQInput,
              factorHD = factorHDInput)
    } else{
      print("PT not run!")
    }

  #----------------STN----------------#
  if(sum(grepl("CQ",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor vs Solid Tissue Normal")>0){
    barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="CQ",
              sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
              prefixRaw = prefixRawInput,
              intFlag = FALSE,
              plotWidthSingle = 5,
              plotWidthCombined = 8,
              errorbarPlotFlag = errorbarPlotFlagInput,
              fileNameString = fileString,
              factorCQ = factorCQInput,
              factorHD = factorHDInput)
    } else{
      print("PT vs STN not run!")
    }

}


runSeqCenterFXN <- function(inputDataDf,
                            intFlagInput=FALSE,
                            prefixRawInput = "rs210PanFinal7575_Nonzero_HiSeq_",
                            fileString = "rs210Pan_Filt7575",
                            errorbarPlotFlagInput = FALSE,
                            factorCQInput = FALSE,
                            factorHDInput = FALSE # HD for host depletion
                            ){
  
  source("00-Functions.R") # for barplotSummaryPerf() and barplotPerf() functions
  
  #----------------All----------------#
  barplotSummaryPerf(inputData = inputDataDf,
                     seqCenterAbbrev="All",
                     sampleTypeInput = "Blood Derived Normal",
                     intFlag = intFlagInput,
                     prefixRaw = prefixRawInput,
                     plotWidthSingle = 8,
                     plotWidthCombined = 8,
                     fileNameString = fileString,
                     factorCQ = factorCQInput,
                     factorHD = factorHDInput)
  
  barplotSummaryPerf(inputData = inputDataDf,
                     seqCenterAbbrev="All",
                     sampleTypeInput = "Primary Tumor",
                     intFlag = intFlagInput,
                     prefixRaw = prefixRawInput,
                     plotWidthSingle = 12,
                     plotWidthCombined = 12,
                     fileNameString = fileString,
                     factorCQ = factorCQInput,
                     factorHD = factorHDInput)
  
  barplotSummaryPerf(inputData = inputDataDf,
                     seqCenterAbbrev="All",
                     sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                     intFlag = intFlagInput,
                     prefixRaw = prefixRawInput,
                     plotWidthSingle = 6,
                     plotWidthCombined = 6,
                     fileNameString = fileString,
                     factorCQ = factorCQInput,
                     factorHD = factorHDInput)
  
  #----------------WGS----------------#
  # HMS
  if(sum(grepl("HMS",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Blood Derived Normal")>0){
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="HMS",
                sampleTypeInput = "Blood Derived Normal",
                prefixRaw = prefixRawInput,
                intFlag = FALSE,
                plotWidthSingle = 5, # 6
                plotWidthCombined = 8,
                # singlePlotHeight = 5,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("HMS BDN not run!")
  }

  if(sum(grepl("HMS",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor")>0){  
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="HMS",
                sampleTypeInput = "Primary Tumor",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 5, # 6
                plotWidthCombined = 8,
                # singlePlotHeight = 5,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("HMS PT not run!")
  }
   
  if(sum(grepl("HMS",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor vs Solid Tissue Normal")>0){  
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="HMS",
                sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 2, # 3
                plotWidthCombined = 5,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("HMS PT vs STN not run!")
  }
  
  
  # BCM
  if(sum(grepl("BCM",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Blood Derived Normal")>0){
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="BCM",
                sampleTypeInput = "Blood Derived Normal",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 3, # 4
                plotWidthCombined = 6,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("BCM BDN not run!")
  }

  if(sum(grepl("BCM",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor")>0){  
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="BCM",
                sampleTypeInput = "Primary Tumor",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 3.5, # 4
                plotWidthCombined = 6,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("BCM PT not run!")
  }

  if(sum(grepl("BCM",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor vs Solid Tissue Normal")>0){    
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="BCM",
                sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 2.5, # 3
                plotWidthCombined = 5,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("BCM PT vs STN not run!")
  }
  
  
  # MDA
  if(sum(grepl("MDA",inputDataDf$datasetName) & 
         inputDataDf$sampleType=="Blood Derived Normal")>0){
    barplotPerf(inputData = inputDataDf,
              seqCenterAbbrev="MDA",
              sampleTypeInput = "Blood Derived Normal",
                prefixRaw = prefixRawInput,
              intFlag = intFlagInput,
              plotWidthSingle = 3.5, # 4
              plotWidthCombined = 6,
              errorbarPlotFlag = errorbarPlotFlagInput,
              fileNameString = fileString,
              factorCQ = factorCQInput,
              factorHD = factorHDInput)
  } else{
    print("MDA BDN not run!")
  }
  
  if(sum(grepl("MDA",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor")>0){
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="MDA",
                sampleTypeInput = "Primary Tumor",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 3.5, # 4
                plotWidthCombined = 6,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("MDA PT not run!")
  }
  
  
  # WashU
  if(sum(grepl("WashU",inputDataDf$datasetName) & 
         inputDataDf$sampleType=="Blood Derived Normal")>0){
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="WashU",
                sampleTypeInput = "Blood Derived Normal",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 2.5, # 3
                plotWidthCombined = 5,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("WashU BDN not run!")
  }
   
  if(sum(grepl("WashU",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor")>0){ 
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="WashU",
                sampleTypeInput = "Primary Tumor",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 2.5, # 3
                plotWidthCombined = 5,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("WashU PT not run!")
  }
  
  # Broad_WGS
  if(sum(grepl("Broad_WGS",inputDataDf$datasetName) & 
         inputDataDf$sampleType=="Blood Derived Normal")>0){
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="Broad_WGS",
                sampleTypeInput = "Blood Derived Normal",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 3.5, # 4
                plotWidthCombined = 6,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("Broad_WGS BDN not run!")
  }

  if(sum(grepl("Broad_WGS",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor")>0){   
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="Broad_WGS",
                sampleTypeInput = "Primary Tumor",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 5, # 6
                plotWidthCombined = 8,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("Broad_WGS PT not run!")
  }
    
  if(sum(grepl("Broad_WGS",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor vs Solid Tissue Normal")>0){
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="Broad_WGS",
                sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 2, # 3
                plotWidthCombined = 5,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("Broad_WGS PT vs STN not run!")
  }
  
  #----------------RNA----------------#
  
  # UNC
  if(sum(grepl("UNC",inputDataDf$datasetName) & 
         inputDataDf$sampleType=="Primary Tumor")>0){
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="UNC",
                sampleTypeInput = "Primary Tumor",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 8,
                plotWidthCombined = 18,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("UNC PT not run!")
  }

  if(sum(grepl("UNC",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor vs Solid Tissue Normal")>0){  
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="UNC",
                sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 6,
                plotWidthCombined = 10,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("UNC PT vs STN not run!")
  }
  
  # CMS
  if(sum(grepl("CMS",inputDataDf$datasetName) & 
           inputDataDf$sampleType=="Primary Tumor")>0){
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="CMS",
                sampleTypeInput = "Primary Tumor",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 2.5, # 3
                plotWidthCombined = 5,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("CMS PT not run!")
  }

  if(sum(grepl("CMS",inputDataDf$datasetName) & 
       inputDataDf$sampleType=="Primary Tumor vs Solid Tissue Normal")>0){  
    barplotPerf(inputData = inputDataDf,
                seqCenterAbbrev="CMS",
                sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                prefixRaw = prefixRawInput,
                intFlag = intFlagInput,
                plotWidthSingle = 2, # 2
                plotWidthCombined = 3,
                errorbarPlotFlag = errorbarPlotFlagInput,
                fileNameString = fileString,
                factorCQ = factorCQInput,
                factorHD = factorHDInput)
  } else{
    print("CMS PT vs STN not run!")
  }
  
}



runSeqCenterEnrichmentFxn2 <- function(pathInput = "Supporting_scripts/S24-RS210-Pan-RawVsCQ-Filt7575-seqcenter/features__",
                            totalFeaturesInput = colnames(rs210PanFinal7575_Nonzero_HiSeq_WGS),
                            datasetRawInput = "rs210PanFinal7575_Nonzero_HiSeq",
                            datasetVSNMInput = "rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ",
                            fileNameStringInput = "rs210Pan_Filt7575",
                            kendallOnlyFlagInput = FALSE,
                            combineFlag = TRUE,
                            combinedWGSonly = FALSE){

  datasetVSNMInput_RNA <- gsub("WGS","RNA",datasetVSNMInput)
  
  source("00-Functions.R") # for enrichmentFxn2() function
  #----------------WGS----------------#
  # HMS
  enrich_HMS_BDN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="HMS",
                sampleType = "Blood Derived Normal",
                plotWidth = 5),
  error=function(w) print("HMS BDN not run!")
  )

  enrich_HMS_PT <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="HMS",
                sampleType = "Primary Tumor",
                plotWidth = 5),
  error=function(w) print("HMS PT not run!")
  )

  enrich_HMS_STN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="HMS",
                sampleType = "Primary Tumor vs Solid Tissue Normal",
                plotWidth = 2.5),
  error=function(w) print("HMS PT vs STN not run!")
  )

  # BCM
  enrich_BCM_BDN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="BCM",
                sampleType = "Blood Derived Normal",
                plotWidth = 4),
  error=function(w) print("BCM BDN not run!")
  )

  enrich_BCM_PT <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="BCM",
                sampleType = "Primary Tumor",
                plotWidth = 4),
  error=function(w) print("BCM PT not run!")
  )

  enrich_BCM_STN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="BCM",
                sampleType = "Primary Tumor vs Solid Tissue Normal",
                plotWidth = 3),
  error=function(w) print("BCM PT vs STN not run!")
  )

  # MDA
  enrich_MDA_BDN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="MDA",
                sampleType = "Blood Derived Normal",
                plotWidth = 4),
  error=function(w) print("MDA BDN not run!")
  )

  enrich_MDA_PT <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="MDA",
                sampleType = "Primary Tumor",
                plotWidth = 4),
  error=function(w) print("MDA PT not run!")
  )

  # WashU
  enrich_WashU_BDN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="WashU",
                sampleType = "Blood Derived Normal",
                plotWidth = 3),
  error=function(w) print("WashU BDN not run!")
  )

  enrich_WashU_PT <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="WashU",
                sampleType = "Primary Tumor",
                plotWidth = 3),
  error=function(w) print("WashU PT not run!")
  )

  # Broad_WGS
  enrich_Broad_WGS_BDN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="Broad_WGS",
                sampleType = "Blood Derived Normal",
                plotWidth = 4),
  error=function(w) print("Broad_WGS BDN not run!")
  )

  enrich_Broad_WGS_PT <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="Broad_WGS",
                sampleType = "Primary Tumor",
                plotWidth = 5),
  error=function(w) print("Broad_WGS PT not run!")
  )

  enrich_Broad_WGS_STN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="Broad_WGS",
                sampleType = "Primary Tumor vs Solid Tissue Normal",
                plotWidth = 2.5),
  error=function(w) print("Broad_WGS PT vs STN not run!")
  )

  #----------------RNA----------------#
  # UNC
  enrich_UNC_PT <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput_RNA,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="UNC",
                sampleType = "Primary Tumor",
                plotWidth = 10),
  error=function(w) print("UNC PT not run!")
  )

  enrich_UNC_STN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput_RNA,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="UNC",
                sampleType = "Primary Tumor vs Solid Tissue Normal",
                plotWidth = 6),
  error=function(w) print("UNC PT vs STN not run!")
  )

  # CMS
  enrich_CMS_PT <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput_RNA,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="CMS",
                sampleType = "Primary Tumor",
                plotWidth = 3),
  error=function(w) print("CMS PT not run!")
  )

  enrich_CMS_STN <- tryCatch(enrichmentFxn2(myPath = pathInput,
                totalFeatures = totalFeaturesInput,
                dataSetRaw = datasetRawInput,
                dataSetVSNM = datasetVSNMInput_RNA,
                fileNameString = fileNameStringInput,
                kendallOnlyFlag = kendallOnlyFlagInput,
                seqCenter="CMS",
                sampleType = "Primary Tumor vs Solid Tissue Normal",
                plotWidth = 2.5),
  error=function(w) print("CMS PT vs STN not run!")
  )

  # print(enrich_UNC_PT)
  # print(length(enrich_CMS_STN))

  if(combineFlag){

    plotPrefix <- paste0("Figures/mlEnrich_",fileNameStringInput,"/")
    # Create folder for plots if doesn't exist
    plotFolder <- paste0("Figures/mlEnrich_",fileNameStringInput)
    if(!( dir.exists( file.path(plotFolder)))){
      dir.create(file.path(plotFolder))
    }

    if(combinedWGSonly){
      enrichList <- list(# HMS
                       enrich_HMS_BDN,
                       enrich_HMS_PT,
                       enrich_HMS_STN,
                       # BCM
                       enrich_BCM_BDN,
                       enrich_BCM_PT,
                       enrich_BCM_STN,
                       # MDA
                       enrich_MDA_BDN,
                       enrich_MDA_PT,
                       # WashU
                       enrich_WashU_BDN,
                       enrich_WashU_PT,
                       # Broad_WGS
                       enrich_Broad_WGS_BDN,
                       enrich_Broad_WGS_PT,
                       enrich_Broad_WGS_STN)
    } else{
      enrichList <- list(# HMS
                       enrich_HMS_BDN,
                       enrich_HMS_PT,
                       enrich_HMS_STN,
                       # BCM
                       enrich_BCM_BDN,
                       enrich_BCM_PT,
                       enrich_BCM_STN,
                       # MDA
                       enrich_MDA_BDN,
                       enrich_MDA_PT,
                       # WashU
                       enrich_WashU_BDN,
                       enrich_WashU_PT,
                       # Broad_WGS
                       enrich_Broad_WGS_BDN,
                       enrich_Broad_WGS_PT,
                       enrich_Broad_WGS_STN,
                       # # UNC
                       # enrich_UNC_PT,
                       # enrich_UNC_STN,
                       # CMS
                       enrich_CMS_PT,
                       enrich_CMS_STN)
    }
    
    # The following line removes all list entries that lack an associated data frame (length > 1)
    # NOTE: If the tryCatch() above fails, it returns a string of length 1
    enrichListFilt <- Filter(length,lapply(enrichList, Filter, f = function(x) length(x) > 1))

    # The following line extracts the Fisher and Kendall test results from every list in the filtered list of lists (enrichListFilt)
    enrichListFiltDf <- do.call(rbind,lapply(enrichListFilt, function(m) m[["fisherKendallCombinedDf"]]))

    #----------------Combine outputs----------------#
    kf_Comb <- enrichListFiltDf

    kf_CombP <- kf_Comb %>%
      group_by(abbrev, ST) %>%
      # NOTE: If combined p-values are smaller than double.xmin, they will be 0 (Inf log), so replace to estimate
      summarise(p.fisher.comb = ifelse(survcomp::combine.test(p.fisher)==0,
                                        yes = .Machine$double.xmin,
                                        no = survcomp::combine.test(p.fisher)),
                p.kendall.comb = ifelse(survcomp::combine.test(p.kendall)==0,
                                        yes = .Machine$double.xmin,
                                        no = survcomp::combine.test(p.kendall)),
                tau.comb = median(tau),
                tau.se = ifelse(is.na(sd(tau)/n()),0,sd(tau)/n()),
                OR.comb = median(OR),
                OR.se = ifelse(is.na(sd(OR)/n()),0,sd(OR)/n()),
                count = n()) %>%
      rstatix::adjust_pvalue("p.fisher.comb", "p.fisher.comb.adj") %>%
      rstatix::adjust_pvalue("p.kendall.comb", "p.kendall.comb.adj") %>%
      rstatix::add_significance(p.col = "p.fisher.comb.adj") %>%
      rstatix::add_significance(p.col = "p.kendall.comb.adj")

    # res <- list(kf_CombP=kf_CombP,
    #             kf_Comb=kf_Comb)
    # return(res)

    fileNameString <- fileNameStringInput
    seqCenter <- "All"
    plotWidthInput <- c(8,8,5)
    sampleTypeInput <- c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal")
    for(ii in 1:length(sampleTypeInput)){
      sampleType <- sampleTypeInput[ii]
      plotWidth <- plotWidthInput[ii]
      kf_CombP_Filt <- kf_CombP %>%
        filter(ST == sampleType)
      
      p.adj.signif.fisher.kendall <- c(kf_CombP_Filt$p.fisher.comb.adj.signif,
                                       kf_CombP_Filt$p.kendall.comb.adj.signif)
      p.adj.signif.kendall.fisher <- c(kf_CombP_Filt$p.kendall.comb.adj.signif,
                                       kf_CombP_Filt$p.fisher.comb.adj.signif)
      
      kf_CombP_Filt %>%
        mutate(log.p.fisher.comb.adj = -log10(p.fisher.comb.adj),
               log.p.kendall.comb.adj = -log10(p.kendall.comb.adj)) %>%
        select(abbrev, log.p.fisher.comb.adj, log.p.kendall.comb.adj,
               p.fisher.comb.adj.signif, p.kendall.comb.adj.signif, count) %>%
        reshape2::melt(id.vars = c("abbrev","p.fisher.comb.adj.signif","p.kendall.comb.adj.signif","count")) %>%
        mutate(variable = factor(case_when(
          variable == "log.p.fisher.comb.adj" ~ "Fisher",
          variable == "log.p.kendall.comb.adj" ~ "Kendall",
        ), levels = c("Kendall","Fisher"))) %>%
        mutate(p.adj.signif.combined = `if`(head(as.character(variable),1)=="Fisher",
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
        geom_text(aes(label = count, y = value),
                  position = position_dodge(0.9), 
                  hjust = 0.5,
                  vjust = 1.5,
                  size = 2,
                  angle = 0,
                  color = "white") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
        ylim(c(0,1.1*max( c(-log10(kf_CombP_Filt$p.fisher.comb.adj),
                            -log10(kf_CombP_Filt$p.kendall.comb.adj)) ))) +
        scale_fill_nejm(name = "Test type") +
        theme_pubr() +
        rotate_x_text(30) +
        theme(axis.text.x = element_text(size=10)) +
        labs(x = "TCGA Cancer Type", 
             y = "-log10(p-adjust)",
             title = paste(seqCenter,sampleType,sep = " | ")) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "right") -> pvalPlot
      
      print(pvalPlot)
      fileNamePval <- paste0(plotPrefix,"pvalue_combined_Fisher_Kendall_",seqCenter,"_",
                             gsub('([[:punct:]])|\\s+','',sampleType),
                             ".jpeg")
      ggsave(filename = fileNamePval,
             plot = pvalPlot,
             dpi = "retina", units = "in", height = 3, width = plotWidth)
      
      kf_CombP_Filt %>%
        ggplot(aes(x = reorder(abbrev, -tau.comb,median), y = tau.comb)) +
        geom_bar(stat="identity", color="black", position=position_dodge()) +
        # geom_errorbar(aes(ymin=tau.comb-tau.se, ymax=tau.comb+tau.se), width=.2) +
        theme_pubr() +
        geom_text(aes(label = count, y = tau.comb),
                  position = position_dodge(0.9), 
                  hjust = 0.5,
                  vjust = 1.5,
                  size = 2,
                  angle = 0,
                  color = "white") +
        rotate_x_text(30) +
        labs(x = "TCGA Cancer Type", 
             y = "Kendall tau",
             title = paste(seqCenter,sampleType,sep = " | ")) +
        geom_text(aes(label = p.kendall.comb.adj.signif, y = tau.comb), vjust = -0.4) +
        ylim(c(0,1.1*max( (kf_CombP_Filt$tau.comb) ))) +
        # ylim(c(0,1.1*max( (kf_CombP_Filt$tau.comb+kf_CombP_Filt$tau.se) ))) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "right") -> barPlotKendallTau
      
      print(barPlotKendallTau)
      fileNameTau <- paste0(plotPrefix,"kendall_combined_tau_",seqCenter,"_",
                            gsub('([[:punct:]])|\\s+','',sampleType),
                            ".jpeg")
      ggsave(filename = fileNameTau,
             plot = barPlotKendallTau,
             dpi = "retina", units = "in", height = 3.5, width = plotWidth)
      
      kf_CombP_Filt %>%
        ggplot(aes(x = reorder(abbrev, -OR.comb,median), y = OR.comb)) +
        geom_bar(stat="identity", color="black", position=position_dodge()) +
        # geom_errorbar(aes(ymin=OR.comb-OR.se, ymax=OR.comb+OR.se), width=.2) +
        theme_pubr() +
        geom_text(aes(label = count, y = OR.comb),
                  position = position_dodge(0.9), 
                  hjust = 0.5,
                  vjust = 1.5,
                  size = 2,
                  angle = 0,
                  color = "white") +
        rotate_x_text(30) +
        labs(x = "TCGA Cancer Type", 
             y = "log10(odds ratio)",
             title = paste(seqCenter,sampleType,sep = " | ")) +
        geom_text(aes(label = p.fisher.comb.adj.signif, y = OR.comb), vjust = -0.4) +
        # ylim(c(0,1.1*max( (kf_CombP_Filt$OR.comb) ))) +
        scale_y_log10(limits = c(1,5*max( (kf_CombP_Filt$OR.comb) ))) +
        # ylim(c(0,1.1*max( (kf_CombP_Filt$OR.comb+kf_CombP_Filt$OR.se) ))) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "right") -> barPlotOR
      
      print(barPlotOR)
      fileNameOR <- paste0(plotPrefix,"odds_ratio_combined_",seqCenter,"_",
                            gsub('([[:punct:]])|\\s+','',sampleType),
                            ".jpeg")
      ggsave(filename = fileNameOR,
             plot = barPlotOR,
             dpi = "retina", units = "in", height = 3.5, width = plotWidth)
    }
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

export2Qiime <- function(metaString = "metaRSFinal5050_Nonzero_HiSeq",
                         countString = "rs210PanFinal5050_Nonzero_HiSeq",
                         dataString = "rs210Pan_Filt5050"){
  
  filePrefix <- paste0("Qiime_data_and_scripts/Qiime_input_data/",dataString,"/")
  # Create folder for plots if doesn't exist
  fileFolder <- paste0("Qiime_data_and_scripts/Qiime_input_data/",dataString)
  if(!( dir.exists( file.path(fileFolder)))){
    dir.create(file.path(fileFolder))
  }
  
  #----------------PT subsets----------------#
  
  # Subset to primary tumor data
  # WGS
  metaData_HMS_PT <- eval(as.name(paste0(metaString,"_HMS"))) %>%
    filter(sample_type == "Primary Tumor") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_BCM_PT <- eval(as.name(paste0(metaString,"_BCM"))) %>%
    filter(sample_type == "Primary Tumor") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_MDA_PT <- eval(as.name(paste0(metaString,"_MDA"))) %>%
    filter(sample_type == "Primary Tumor") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_WashU_PT <- eval(as.name(paste0(metaString,"_WashU"))) %>%
    filter(sample_type == "Primary Tumor") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_Broad_WGS_PT <- eval(as.name(paste0(metaString,"_Broad_WGS"))) %>%
    filter(sample_type == "Primary Tumor") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  # RNA
  metaData_UNC_PT <- eval(as.name(paste0(metaString,"_UNC"))) %>%
    filter(sample_type == "Primary Tumor") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_CMS_PT <- eval(as.name(paste0(metaString,"_CMS"))) %>%
    filter(sample_type == "Primary Tumor") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  
  ## Write metadata to text files
  # WGS
  write.table(metaData_HMS_PT, 
              file = paste0(filePrefix,"metaData_HMS_PT.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_BCM_PT, 
              file = paste0(filePrefix,"metaData_BCM_PT.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_MDA_PT, 
              file = paste0(filePrefix,"metaData_MDA_PT.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_WashU_PT, 
              file = paste0(filePrefix,"metaData_WashU_PT.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_Broad_WGS_PT, 
              file = paste0(filePrefix,"metaData_Broad_WGS_PT.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  # RNA
  write.table(metaData_UNC_PT, 
              file = paste0(filePrefix,"metaData_UNC_PT.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_CMS_PT, 
              file = paste0(filePrefix,"metaData_CMS_PT.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  
  ## Subset count data
  # WGS
  countData_HMS_PT <- eval(as.name(paste0(countString,"_HMS")))[metaData_HMS_PT$sampleid,]
  countData_BCM_PT <- eval(as.name(paste0(countString,"_BCM")))[metaData_BCM_PT$sampleid,]
  countData_MDA_PT <- eval(as.name(paste0(countString,"_MDA")))[metaData_MDA_PT$sampleid,]
  countData_WashU_PT <- eval(as.name(paste0(countString,"_WashU")))[metaData_WashU_PT$sampleid,]
  countData_Broad_WGS_PT <- eval(as.name(paste0(countString,"_Broad_WGS")))[metaData_Broad_WGS_PT$sampleid,]
  # RNA
  countData_UNC_PT <- eval(as.name(paste0(countString,"_UNC")))[metaData_UNC_PT$sampleid,]
  countData_CMS_PT <- eval(as.name(paste0(countString,"_CMS")))[metaData_CMS_PT$sampleid,]
  
  ## Save count data as biom tables
  # WGS
  countData_HMS_PT_BIOM <- make_biom(t(countData_HMS_PT))
  write_biom(countData_HMS_PT_BIOM, biom_file = paste0(filePrefix,"countData_HMS_PT.biom"))
  
  countData_BCM_PT_BIOM <- make_biom(t(countData_BCM_PT))
  write_biom(countData_BCM_PT_BIOM, biom_file = paste0(filePrefix,"countData_BCM_PT.biom"))
  
  countData_MDA_PT_BIOM <- make_biom(t(countData_MDA_PT))
  write_biom(countData_MDA_PT_BIOM, biom_file = paste0(filePrefix,"countData_MDA_PT.biom"))
  
  countData_WashU_PT_BIOM <- make_biom(t(countData_WashU_PT))
  write_biom(countData_WashU_PT_BIOM, biom_file = paste0(filePrefix,"countData_WashU_PT.biom"))
  
  countData_Broad_WGS_PT_BIOM <- make_biom(t(countData_Broad_WGS_PT))
  write_biom(countData_Broad_WGS_PT_BIOM, biom_file = paste0(filePrefix,"countData_Broad_WGS_PT.biom"))
  # RNA
  countData_UNC_PT_BIOM <- make_biom(t(countData_UNC_PT))
  write_biom(countData_UNC_PT_BIOM, biom_file = paste0(filePrefix,"countData_UNC_PT.biom"))
  
  countData_CMS_PT_BIOM <- make_biom(t(countData_CMS_PT))
  write_biom(countData_CMS_PT_BIOM, biom_file = paste0(filePrefix,"countData_CMS_PT.biom"))
  
  
  #----------------BDN subsets----------------#
  
  # Subset to primary tumor data
  # WGS
  metaData_HMS_BDN <- eval(as.name(paste0(metaString,"_HMS"))) %>%
    filter(sample_type == "Blood Derived Normal") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_BCM_BDN <- eval(as.name(paste0(metaString,"_BCM"))) %>%
    filter(sample_type == "Blood Derived Normal") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_MDA_BDN <- eval(as.name(paste0(metaString,"_MDA"))) %>%
    filter(sample_type == "Blood Derived Normal") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_WashU_BDN <- eval(as.name(paste0(metaString,"_WashU"))) %>%
    filter(sample_type == "Blood Derived Normal") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_Broad_WGS_BDN <- eval(as.name(paste0(metaString,"_Broad_WGS"))) %>%
    filter(sample_type == "Blood Derived Normal") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  
  ## Write metadata to text files
  # WGS
  write.table(metaData_HMS_BDN, 
              file = paste0(filePrefix,"metaData_HMS_BDN.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_BCM_BDN, 
              file = paste0(filePrefix,"metaData_BCM_BDN.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_MDA_BDN, 
              file = paste0(filePrefix,"metaData_MDA_BDN.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_WashU_BDN, 
              file = paste0(filePrefix,"metaData_WashU_BDN.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_Broad_WGS_BDN, 
              file = paste0(filePrefix,"metaData_Broad_WGS_BDN.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  
  ## Subset count data
  # WGS
  countData_HMS_BDN <- eval(as.name(paste0(countString,"_HMS")))[metaData_HMS_BDN$sampleid,]
  countData_BCM_BDN <- eval(as.name(paste0(countString,"_BCM")))[metaData_BCM_BDN$sampleid,]
  countData_MDA_BDN <- eval(as.name(paste0(countString,"_MDA")))[metaData_MDA_BDN$sampleid,]
  countData_WashU_BDN <- eval(as.name(paste0(countString,"_WashU")))[metaData_WashU_BDN$sampleid,]
  countData_Broad_WGS_BDN <- eval(as.name(paste0(countString,"_Broad_WGS")))[metaData_Broad_WGS_BDN$sampleid,]
  
  ## Save count data as biom tables
  # WGS
  countData_HMS_BDN_BIOM <- make_biom(t(countData_HMS_BDN))
  write_biom(countData_HMS_BDN_BIOM, biom_file = paste0(filePrefix,"countData_HMS_BDN.biom"))
  
  countData_BCM_BDN_BIOM <- make_biom(t(countData_BCM_BDN))
  write_biom(countData_BCM_BDN_BIOM, biom_file = paste0(filePrefix,"countData_BCM_BDN.biom"))
  
  countData_MDA_BDN_BIOM <- make_biom(t(countData_MDA_BDN))
  write_biom(countData_MDA_BDN_BIOM, biom_file = paste0(filePrefix,"countData_MDA_BDN.biom"))
  
  countData_WashU_BDN_BIOM <- make_biom(t(countData_WashU_BDN))
  write_biom(countData_WashU_BDN_BIOM, biom_file = paste0(filePrefix,"countData_WashU_BDN.biom"))
  
  countData_Broad_WGS_BDN_BIOM <- make_biom(t(countData_Broad_WGS_BDN))
  write_biom(countData_Broad_WGS_BDN_BIOM, biom_file = paste0(filePrefix,"countData_Broad_WGS_BDN.biom"))
  
  print("#-------PT counts-------#")
  print("HMS")
  print(summary(rowSums(countData_HMS_PT)))
  print("BCM")
  print(summary(rowSums(countData_BCM_PT)))
  print("MDA")
  print(summary(rowSums(countData_MDA_PT)))
  print("WashU")
  print(summary(rowSums(countData_WashU_PT)))
  print("Broad_WGS")
  print(summary(rowSums(countData_Broad_WGS_PT)))
  print("UNC")
  print(summary(rowSums(countData_UNC_PT)))
  print("CMS")
  print(summary(rowSums(countData_CMS_PT)))
  print("#-------BDN counts-------#")
  print("HMS")
  print(summary(rowSums(countData_HMS_BDN)))
  print("BCM")
  print(summary(rowSums(countData_BCM_BDN)))
  print("MDA")
  print(summary(rowSums(countData_MDA_BDN)))
  print("WashU")
  print(summary(rowSums(countData_WashU_BDN)))
  print("Broad_WGS")
  print(summary(rowSums(countData_Broad_WGS_BDN)))
  
}

export2QiimeCQ <- function(metaString = "metaRSFinal5050_Nonzero_HiSeq",
                         countString = "rs210PanFinal5050_Nonzero_HiSeq",
                         dataString = "rs210Pan_Filt5050"){
  
  filePrefix <- paste0("Qiime_data_and_scripts/Qiime_input_data/",dataString,"/")
  # Create folder for plots if doesn't exist
  fileFolder <- paste0("Qiime_data_and_scripts/Qiime_input_data/",dataString)
  if(!( dir.exists( file.path(fileFolder)))){
    dir.create(file.path(fileFolder))
  }
  
  #----------------PT subsets----------------#
  
  # Subset to primary tumor data
  metaData_WGS_PT <- eval(as.name(paste0(metaString,"_WGS"))) %>%
    filter(sample_type == "Primary Tumor") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  metaData_RNA_PT <- eval(as.name(paste0(metaString,"_RNA"))) %>%
    filter(sample_type == "Primary Tumor") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  
  ## Write metadata to text files
  write.table(metaData_WGS_PT, 
              file = paste0(filePrefix,"metaData_WGS_PT.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  write.table(metaData_RNA_PT, 
              file = paste0(filePrefix,"metaData_RNA_PT.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  
  ## Subset count data
  # WGS
  countData_WGS_PT <- eval(as.name(paste0(countString,"_WGS_CQ")))[metaData_WGS_PT$sampleid,]
  countData_RNA_PT <- eval(as.name(paste0(countString,"_RNA_CQ")))[metaData_RNA_PT$sampleid,]
  
  ## Save count data as biom tables
  countData_WGS_PT_BIOM <- make_biom(t(countData_WGS_PT))
  write_biom(countData_WGS_PT_BIOM, biom_file = paste0(filePrefix,"countData_WGS_PT.biom"))
  
  countData_RNA_PT_BIOM <- make_biom(t(countData_RNA_PT))
  write_biom(countData_RNA_PT_BIOM, biom_file = paste0(filePrefix,"countData_RNA_PT.biom"))
  
  #----------------BDN subsets----------------#
  
  # Subset to primary tumor data
  # WGS
  metaData_WGS_BDN <- eval(as.name(paste0(metaString,"_WGS"))) %>%
    filter(sample_type == "Blood Derived Normal") %>%
    select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
    rownames_to_column("sampleid") %>% droplevels()
  
  ## Write metadata to text files
  write.table(metaData_WGS_BDN, 
              file = paste0(filePrefix,"metaData_WGS_BDN.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  
  ## Subset count data
  # WGS
  countData_WGS_BDN <- eval(as.name(paste0(countString,"_WGS_CQ")))[metaData_WGS_BDN$sampleid,]
  
  ## Save count data as biom tables
  # WGS
  countData_WGS_BDN_BIOM <- make_biom(t(countData_WGS_BDN))
  write_biom(countData_WGS_BDN_BIOM, biom_file = paste0(filePrefix,"countData_WGS_BDN.biom"))
  
  print("#-------PT counts-------#")
  print("WGS")
  print(summary(rowSums(countData_WGS_PT)))
  print("RNA")
  print(summary(rowSums(countData_RNA_PT)))
  print("#-------BDN counts-------#")
  print("WGS")
  print(summary(rowSums(countData_WGS_BDN)))
}

alphaBetaFXN <- function(metaData = metaRSFinal9090_Nonzero_HiSeq_HMS,
                        countData = rs210PanFinal9090_Nonzero_HiSeq_HMS,
                        alphaDivType = "Observed",
                        dataString = "rs210Pan_Filt9090_HMS",
                        useTaxTable = TRUE,
                        alphaPlotWidth = 8,
                        alphaPlotHeight = 5,
                        betaPlotWidth = 5,
                        betaPlotHeight = 5,
                        ptOnlyFlag = FALSE,
                        alphaOnlyFlag = FALSE)
{
  require(tibble)
  require(ggpubr)
  require(ggsci)
  
  
  plotPrefix <- paste0("Figures/alphaBeta_",dataString,"/")
  
  # Create folder for plots if doesn't exist
  plotFolder <- paste0("Figures/alphaBeta_",dataString)
  if(!( dir.exists( file.path(plotFolder)))){
    dir.create(file.path(plotFolder))
  }
  
  # Note that taxa and count tables have to be a matrix
  if(useTaxTable){
    psTaxTable <- taxRS210_ff_combTaxaSpeciesZebraWithViruses %>%
    select(-UNITN,-UHGG,-PATH,-WIS,-cov_All) %>%
    column_to_rownames("OGU") %>%
    rename(Domain=domain,Phylum=phylum,Class=class,Order=order,
           Family=family,Genus=genus,Species=species) %>%
    as.matrix()
  }
  
  if(ptOnlyFlag){
    
    metaDataPT <- metaData %>% filter(sample_type == "Primary Tumor") %>% droplevels()
    countDataPT <- as.matrix(countData[rownames(metaDataPT),])
    
    if(useTaxTable){
      psObjPT <- phyloseq(otu_table(countDataPT, taxa_are_rows = FALSE), 
                        sample_data(metaDataPT), tax_table(psTaxTable))
    } else{
      psObjPT <- phyloseq(otu_table(countDataPT, taxa_are_rows = FALSE), 
                      sample_data(metaDataPT))
    }
    
    
    ## Rarefy
    print("Primary tumor read distribution")
    print(summary(sample_sums(psObjPT)))
    rarefyLevelPT <- readline(prompt="Enter rarefaction amount for PT: ")
    rarefyLevelPTInt <- as.integer(rarefyLevelPT)
    
    psObjPT_rare <- rarefy_even_depth(psObjPT, sample.size = rarefyLevelPTInt,
                                      rngseed = 42, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
    if(length(table(sample_data(psObjPT_rare)$investigation))==1){
      print("Rarefied data only has 1 cancer type")
      return(NA)
    }
    
    ## Calculate alpha diversity
    alphaDivPT <- data.frame(estimate_richness(psObjPT_rare, 
                                               measures = c("Observed","Chao1", 
                                                            "ACE", "Shannon", 
                                                            "Simpson", "InvSimpson"))) %>%
      mutate(rarefactionLevel = rarefyLevelPTInt) %>%
      rownames_to_column("sampleid") %>% mutate(sampleid=gsub("^X","",sampleid)) %>%
      mutate(investigation = factor(gsub("^TCGA-","",metaDataPT[sampleid,c("investigation")]))) %>%
      column_to_rownames("sampleid")
    
    ## Plot
    alphaDivPT %>%
      mutate(investigation = reorder(investigation, !!sym(alphaDivType), median)) %>%
      ggboxplot(x = "investigation",
                y = alphaDivType,
                fill = "investigation",
                palette = "igv",
                legend = "none",
                xlab = "") +
      stat_compare_means(label.x.npc = 0.05) -> plotAlphaPT
    ggsave(plot = plotAlphaPT,
           filename = paste0(plotPrefix,"alphaDiv_PT_R",rarefyLevelPTInt,"_",alphaDivType,".jpeg"),
           dpi = "retina", units = "in", width = alphaPlotWidth, height = alphaPlotHeight)
    
    ## Calculate beta diversity
    # Bray-Curtis PT
    brayDistPT = phyloseq::distance(psObjPT_rare, method="bray") # need to rarefy
    brayOrdinationPT = ordinate(psObjPT_rare, method="PCoA", distance=brayDistPT)
    print("Calculating PT Bray-Curtis adonis...")
    brayAdonisPT <- vegan::adonis2(brayDistPT ~ sample_data(psObjPT_rare)$investigation)
    brayAdonisPTString <- sprintf("F=%.1f, p=%0.3f",brayAdonisPT$F[1], brayAdonisPT$`Pr(>F)`[1])
    
    brayPlot2D_PT <- plot_ordination(psObjPT_rare, brayOrdinationPT, color="investigation", axes = c(1,2)) +
      theme(aspect.ratio=1) + theme_pubr() + geom_point(size = 3) + scale_color_igv() + coord_fixed() + theme(legend.position = "right") +
      labs(color="Disease type") + 
      annotate(geom="text", x=Inf,y=Inf, hjust=1, vjust=1, label = brayAdonisPTString)
    ggsave(plot = brayPlot2D_PT,
           filename = paste0(plotPrefix,"brayCurtis_PT_R",rarefyLevelPTInt,".jpeg"),
           dpi = "retina", units = "in", width = betaPlotWidth, height = betaPlotHeight)
    
    # Jaccard PT
    jaccardDistPT = phyloseq::distance(psObjPT_rare, method="jaccard") # need to rarefy
    jaccardOrdinationPT = ordinate(psObjPT_rare, method="PCoA", distance=jaccardDistPT)
    print("Calculating PT Jaccard adonis...")
    jaccardAdonisPT <- vegan::adonis2(jaccardDistPT ~ sample_data(psObjPT_rare)$investigation)
    jaccardAdonisPTString <- sprintf("F=%.1f, p=%0.3f",jaccardAdonisPT$F[1], jaccardAdonisPT$`Pr(>F)`[1])
    
    jaccardPlot2D_PT <- plot_ordination(psObjPT_rare, jaccardOrdinationPT, color="investigation", axes = c(1,2)) +
      theme(aspect.ratio=1) + theme_pubr() + geom_point(size = 3) + scale_color_igv() + coord_fixed() + theme(legend.position = "right") +
      labs(color="Disease type") + 
      annotate(geom="text", x=Inf,y=Inf, hjust=1, vjust=1, label = jaccardAdonisPTString)
    ggsave(plot = jaccardPlot2D_PT,
           filename = paste0(plotPrefix,"jaccard_PT_R",rarefyLevelPTInt,".jpeg"),
           dpi = "retina", units = "in", width = betaPlotWidth, height = betaPlotHeight)
    
    res <- list(alphaDivPT=alphaDivPT,
                # Phyloseq objs
                psObjPT_rare=psObjPT_rare,
                # Bray curtis
                brayDistPT=brayDistPT,
                brayAdonisPT=brayAdonisPT,
                # Jaccard
                jaccardDistPT=jaccardDistPT,
                jaccardAdonisPT=jaccardAdonisPT)
    save(res, file = paste0(plotPrefix,"alphaBeta_Objects_PT_R",rarefyLevelPTInt,".RData"))
    return(res)
    
    
  } else{
    
    metaDataPT <- metaData %>% filter(sample_type == "Primary Tumor") %>% droplevels()
    metaDataBDN <- metaData %>% filter(sample_type == "Blood Derived Normal") %>% droplevels()
    countDataPT <- as.matrix(countData[rownames(metaDataPT),])
    countDataBDN <- as.matrix(countData[rownames(metaDataBDN),])
    
    if(useTaxTable){
      psObjPT <- phyloseq(otu_table(countDataPT, taxa_are_rows = FALSE), 
                        sample_data(metaDataPT), tax_table(psTaxTable))
      psObjBDN <- phyloseq(otu_table(countDataBDN, taxa_are_rows = FALSE), 
                           sample_data(metaDataBDN), tax_table(psTaxTable))
      } else{
        psObjPT <- phyloseq(otu_table(countDataPT, taxa_are_rows = FALSE), 
                        sample_data(metaDataPT))
        psObjBDN <- phyloseq(otu_table(countDataBDN, taxa_are_rows = FALSE), 
                             sample_data(metaDataBDN))
      }
    
    
    ## Rarefy
    print("Primary tumor read distribution")
    print(summary(sample_sums(psObjPT)))
    rarefyLevelPT <- readline(prompt="Enter rarefaction amount for PT: ")
    rarefyLevelPTInt <- as.integer(rarefyLevelPT)
    
    print("Blood read distribution")
    print(summary(sample_sums(psObjBDN)))
    rarefyLevelBDN <- readline(prompt="Enter rarefaction amount for BDN: ")
    rarefyLevelBDNInt <- as.integer(rarefyLevelBDN)
    
    psObjPT_rare <- rarefy_even_depth(psObjPT, sample.size = rarefyLevelPTInt,
                                      rngseed = 42, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
    psObjBDN_rare <- rarefy_even_depth(psObjBDN, sample.size = rarefyLevelBDNInt,
                                       rngseed = 42, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
    
    if(length(table(sample_data(psObjPT_rare)$investigation))==1){
      print("Error: Rarefied PT data only has 1 cancer type")
      return(NA)
    }
    if(length(table(sample_data(psObjBDN_rare)$investigation))==1){
      print("Error: Rarefied PT data only has 1 cancer type")
      return(NA)
    }
    
    ## Calculate alpha diversity
    alphaDivPT <- data.frame(estimate_richness(psObjPT_rare, 
                                               measures = c("Observed","Chao1", 
                                                            "ACE", "Shannon", 
                                                            "Simpson", "InvSimpson"))) %>%
      mutate(rarefactionLevel = rarefyLevelPTInt) %>%
      rownames_to_column("sampleid") %>% mutate(sampleid=gsub("^X","",sampleid)) %>%
      mutate(investigation = factor(gsub("^TCGA-","",metaDataPT[sampleid,c("investigation")]))) %>%
      column_to_rownames("sampleid")
    alphaDivBDN <- data.frame(estimate_richness(psObjBDN_rare, 
                                                measures = c("Observed","Chao1", 
                                                             "ACE", "Shannon", 
                                                             "Simpson", "InvSimpson"))) %>%
      mutate(rarefactionLevel = rarefyLevelBDNInt) %>%
      rownames_to_column("sampleid") %>% mutate(sampleid=gsub("^X","",sampleid)) %>% 
      mutate(investigation = factor(gsub("^TCGA-","",metaDataBDN[sampleid,c("investigation")]))) %>%
      column_to_rownames("sampleid")
    
    ## Plot
    alphaDivPT %>%
      mutate(investigation = reorder(investigation, !!sym(alphaDivType), median)) %>%
      ggboxplot(x = "investigation",
                y = alphaDivType,
                fill = "investigation",
                palette = "igv",
                legend = "none",
                xlab = "") +
      stat_compare_means(label.x.npc = 0.05) -> plotAlphaPT
    ggsave(plot = plotAlphaPT,
           filename = paste0(plotPrefix,"alphaDiv_PT_R",rarefyLevelPTInt,"_",alphaDivType,".jpeg"),
           dpi = "retina", units = "in", width = alphaPlotWidth, height = alphaPlotHeight)
    
    alphaDivBDN %>%
      mutate(investigation = reorder(investigation, !!sym(alphaDivType), median)) %>%
      ggboxplot(x = "investigation",
                y = alphaDivType,
                fill = "investigation",
                palette = "igv",
                legend = "none",
                xlab = "") +
      stat_compare_means(label.x.npc = 0.05) -> plotAlphaBDN
    ggsave(plot = plotAlphaBDN,
           filename = paste0(plotPrefix,"alphaDiv_BDN_R",rarefyLevelPTInt,"_",alphaDivType,".jpeg"),
           dpi = "retina", units = "in", width = alphaPlotWidth, height = alphaPlotHeight)
    
    ## Calculate beta diversity
    # Bray-Curtis PT
    brayDistPT = phyloseq::distance(psObjPT_rare, method="bray") # need to rarefy
    brayOrdinationPT = ordinate(psObjPT_rare, method="PCoA", distance=brayDistPT)
    print("Calculating PT Bray-Curtis adonis...")
    brayAdonisPT <- vegan::adonis2(brayDistPT ~ sample_data(psObjPT_rare)$investigation)
    brayAdonisPTString <- sprintf("F=%.1f, p=%0.3f",brayAdonisPT$F[1], brayAdonisPT$`Pr(>F)`[1])
    
    brayPlot2D_PT <- plot_ordination(psObjPT_rare, brayOrdinationPT, color="investigation", axes = c(1,2)) +
      theme(aspect.ratio=1) + theme_pubr() + geom_point(size = 3) + scale_color_igv() + coord_fixed() + theme(legend.position = "right") +
      labs(color="Disease type") + 
      annotate(geom="text", x=Inf,y=Inf, hjust=1, vjust=1, label = brayAdonisPTString)
    ggsave(plot = brayPlot2D_PT,
           filename = paste0(plotPrefix,"brayCurtis_PT_R",rarefyLevelPTInt,".jpeg"),
           dpi = "retina", units = "in", width = betaPlotWidth, height = betaPlotHeight)
    
    # Bray-Curtis BDN
    brayDistBDN = phyloseq::distance(psObjBDN_rare, method="bray") # need to rarefy
    brayOrdinationBDN = ordinate(psObjBDN_rare, method="PCoA", distance=brayDistBDN)
    print("Calculating BDN Bray-Curtis adonis...")
    brayAdonisBDN <- vegan::adonis2(brayDistBDN ~ sample_data(psObjBDN_rare)$investigation)
    brayAdonisBDNString <- sprintf("F=%.1f, p=%0.3f",brayAdonisBDN$F[1], brayAdonisBDN$`Pr(>F)`[1])
    
    brayPlot2D_BDN <- plot_ordination(psObjBDN_rare, brayOrdinationBDN, color="investigation", axes = c(1,2)) +
      theme(aspect.ratio=1) + theme_pubr() + geom_point(size = 3) + scale_color_igv() + coord_fixed() + theme(legend.position = "right") +
      labs(color="Disease type") + 
      annotate(geom="text", x=Inf,y=Inf, hjust=1, vjust=1, label = brayAdonisBDNString)
    ggsave(plot = brayPlot2D_BDN,
           filename = paste0(plotPrefix,"brayCurtis_BDN_R",rarefyLevelBDNInt,".jpeg"),
           dpi = "retina", units = "in", width = betaPlotWidth, height = betaPlotHeight)
    
    # Jaccard PT
    jaccardDistPT = phyloseq::distance(psObjPT_rare, method="jaccard") # need to rarefy
    jaccardOrdinationPT = ordinate(psObjPT_rare, method="PCoA", distance=jaccardDistPT)
    print("Calculating PT Jaccard adonis...")
    jaccardAdonisPT <- vegan::adonis2(jaccardDistPT ~ sample_data(psObjPT_rare)$investigation)
    jaccardAdonisPTString <- sprintf("F=%.1f, p=%0.3f",jaccardAdonisPT$F[1], jaccardAdonisPT$`Pr(>F)`[1])
    
    jaccardPlot2D_PT <- plot_ordination(psObjPT_rare, jaccardOrdinationPT, color="investigation", axes = c(1,2)) +
      theme(aspect.ratio=1) + theme_pubr() + geom_point(size = 3) + scale_color_igv() + coord_fixed() + theme(legend.position = "right") +
      labs(color="Disease type") + 
      annotate(geom="text", x=Inf,y=Inf, hjust=1, vjust=1, label = jaccardAdonisPTString)
    ggsave(plot = jaccardPlot2D_PT,
           filename = paste0(plotPrefix,"jaccard_PT_R",rarefyLevelPTInt,".jpeg"),
           dpi = "retina", units = "in", width = betaPlotWidth, height = betaPlotHeight)
    
    # Jaccard BDN
    jaccardDistBDN = phyloseq::distance(psObjBDN_rare, method="jaccard") # need to rarefy
    jaccardOrdinationBDN = ordinate(psObjBDN_rare, method="PCoA", distance=jaccardDistBDN)
    print("Calculating BDN Jaccard adonis...")
    jaccardAdonisBDN <- vegan::adonis2(jaccardDistBDN ~ sample_data(psObjBDN_rare)$investigation)
    jaccardAdonisBDNString <- sprintf("F=%.1f, p=%0.3f",jaccardAdonisBDN$F[1], jaccardAdonisBDN$`Pr(>F)`[1])
    
    jaccardPlot2D_BDN <- plot_ordination(psObjBDN_rare, jaccardOrdinationBDN, color="investigation", axes = c(1,2)) +
      theme(aspect.ratio=1) + theme_pubr() + geom_point(size = 3) + scale_color_igv() + coord_fixed() + theme(legend.position = "right") +
      labs(color="Disease type") + 
      annotate(geom="text", x=Inf,y=Inf, hjust=1, vjust=1, label = jaccardAdonisBDNString)
    ggsave(plot = jaccardPlot2D_BDN,
           filename = paste0(plotPrefix,"jaccard_BDN_R",rarefyLevelBDNInt,".jpeg"),
           dpi = "retina", units = "in", width = betaPlotWidth, height = betaPlotHeight)
    
    res <- list(alphaDivPT=alphaDivPT,
                alphaDivBDN=alphaDivBDN,
                # Phyloseq objs
                psObjPT_rare=psObjPT_rare,
                psObjBDN_rare=psObjBDN_rare,
                # Bray curtis
                brayDistPT=brayDistPT,
                brayDistBDN=brayDistBDN,
                brayAdonisPT=brayAdonisPT,
                brayAdonisBDN=brayAdonisBDN,
                # Jaccard
                jaccardDistPT=jaccardDistPT,
                jaccardDistBDN=jaccardDistBDN,
                jaccardAdonisPT=jaccardAdonisPT,
                jaccardAdonisBDN=jaccardAdonisBDN)
    save(res, file = paste0(plotPrefix,"alphaBeta_Objects_PT_R",
                            rarefyLevelPTInt,"BDN_R",rarefyLevelBDNInt,".RData"))
    return(res)
    
  }
  
}

runAlphaBetaSeqCenter <- function(metaString = "metaRSFinal9090_Nonzero_HiSeq",
                                  countString = "rs210PanFinal9090_Nonzero_HiSeq",
                                  dataStringInput = "rs210Pan_Filt9090",
                                  useTaxTableInput = TRUE){
  print("Working on HMS...")
  ab_HMS <- alphaBetaFXN(metaData = eval(as.name(paste0(metaString,"_HMS"))),
                         countData = eval(as.name(paste0(countString,"_HMS"))),
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataStringInput,"_HMS"),
                         alphaPlotWidth = 8)
  print("Working on BCM...")
  ab_BCM <- alphaBetaFXN(metaData = eval(as.name(paste0(metaString,"_BCM"))),
                         countData = eval(as.name(paste0(countString,"_BCM"))),
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataStringInput,"_BCM"),
                         alphaPlotWidth = 8)
  print("Working on MDA...")
  ab_MDA <- alphaBetaFXN(metaData = eval(as.name(paste0(metaString,"_MDA"))),
                         countData = eval(as.name(paste0(countString,"_MDA"))),
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataStringInput,"_MDA"),
                         alphaPlotWidth = 8)
  print("Working on WashU...")
  ab_WashU <- alphaBetaFXN(metaData = eval(as.name(paste0(metaString,"_WashU"))),
                         countData = eval(as.name(paste0(countString,"_WashU"))),
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataStringInput,"_WashU"),
                         alphaPlotWidth = 8)
  print("Working on Broad_WGS...")
  ab_Broad_WGS <- alphaBetaFXN(metaData = eval(as.name(paste0(metaString,"_Broad_WGS"))),
                         countData = eval(as.name(paste0(countString,"_Broad_WGS"))),
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataStringInput,"_Broad_WGS"),
                         alphaPlotWidth = 8)
  print("Working on UNC...")
  ab_UNC <- alphaBetaFXN(metaData = eval(as.name(paste0(metaString,"_UNC"))),
                         countData = eval(as.name(paste0(countString,"_UNC"))),
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataStringInput,"_UNC"),
                         alphaPlotWidth = 8,
                         ptOnlyFlag = TRUE)
  print("Working on CMS...")
  ab_CMS <- alphaBetaFXN(metaData = eval(as.name(paste0(metaString,"_CMS"))),
                         countData = eval(as.name(paste0(countString,"_CMS"))),
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataStringInput,"_CMS"),
                         alphaPlotWidth = 8,
                         ptOnlyFlag = TRUE)
}

runAlphaBetaSeqCenter_TaxaLevel <- function(metaData = metaRSFinal5050_Nonzero_HiSeq,
                                            countData = rs210PanFinal5050_Nonzero_HiSeq,
                                            dataStringInput = "rs210Pan_Filt5050",
                                            taxaLevel = "Species",
                                            useTaxTableInput = FALSE){
  require(phyloseq)
  require(microbiome)
  dataString = paste0(dataStringInput,"_",taxaLevel)
  rsTaxTable <- taxRS210_ff_combTaxaSpeciesZebraWithViruses %>%
    select(-UNITN,-UHGG,-PATH,-WIS,-cov_All) %>%
    column_to_rownames("OGU") %>%
    rename(Domain=domain,Phylum=phylum,Class=class,Order=order,
           Family=family,Genus=genus,Species=species) %>%
    as.matrix()
  
  psObj <- phyloseq(otu_table(countData, taxa_are_rows = FALSE), 
                    sample_data(metaData), tax_table(rsTaxTable))
  psObj_aggr = aggregate_taxa(psObj, taxaLevel)
  countData_aggr <- data.frame(t(otu_table(psObj_aggr)))
  
  # Subset metadata
  metaData_HMS <- metaData %>% 
    filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
  metaData_BCM <- metaData %>% 
    filter(data_submitting_center_label == "Baylor College of Medicine") %>% droplevels()
  metaData_MDA <- metaData %>% 
    filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% droplevels()
  metaData_WashU <- metaData %>% 
    filter(data_submitting_center_label == "Washington University School of Medicine") %>% droplevels()
  metaData_Broad_WGS <- metaData %>% 
    filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
    filter(experimental_strategy == "WGS") %>% droplevels()
  metaData_UNC <- metaData %>% 
    filter(data_submitting_center_label == "University of North Carolina") %>% droplevels()
  metaData_CMS <- metaData %>% 
    filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% droplevels()
  
  # Subset count data
  countData_aggr_HMS <- countData_aggr[rownames(metaData_HMS),]
  countData_aggr_BCM <- countData_aggr[rownames(metaData_BCM),]
  countData_aggr_MDA <- countData_aggr[rownames(metaData_MDA),]
  countData_aggr_WashU <- countData_aggr[rownames(metaData_WashU),]
  countData_aggr_Broad_WGS <- countData_aggr[rownames(metaData_Broad_WGS),]
  countData_aggr_UNC <- countData_aggr[rownames(metaData_UNC),]
  countData_aggr_CMS <- countData_aggr[rownames(metaData_CMS),]
  
  print("Working on HMS...")
  ab_HMS <- alphaBetaFXN(metaData = metaData_HMS,
                         countData = countData_aggr_HMS,
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataString,"_HMS"),
                         alphaPlotWidth = 8)
  print("Working on BCM...")
  ab_BCM <- alphaBetaFXN(metaData = metaData_BCM,
                         countData = countData_aggr_BCM,
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataString,"_BCM"),
                         alphaPlotWidth = 8)
  print("Working on MDA...")
  ab_MDA <- alphaBetaFXN(metaData = metaData_MDA,
                         countData = countData_aggr_MDA,
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataString,"_MDA"),
                         alphaPlotWidth = 8)
  print("Working on WashU...")
  ab_WashU <- alphaBetaFXN(metaData = metaData_WashU,
                           countData = countData_aggr_WashU,
                           alphaDivType = "Observed",
                           useTaxTable = useTaxTableInput,
                           dataString = paste0(dataString,"_WashU"),
                           alphaPlotWidth = 8)
  print("Working on Broad_WGS...")
  ab_Broad_WGS <- alphaBetaFXN(metaData = metaData_Broad_WGS,
                               countData = countData_aggr_Broad_WGS,
                               alphaDivType = "Observed",
                               useTaxTable = useTaxTableInput,
                               dataString = paste0(dataString,"_Broad_WGS"),
                               alphaPlotWidth = 8)
  print("Working on UNC...")
  ab_UNC <- alphaBetaFXN(metaData = metaData_UNC,
                         countData = countData_aggr_UNC,
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataString,"_UNC"),
                         alphaPlotWidth = 8,
                         ptOnlyFlag = TRUE)
  print("Working on CMS...")
  ab_CMS <- alphaBetaFXN(metaData = metaData_CMS,
                         countData = countData_aggr_CMS,
                         alphaDivType = "Observed",
                         useTaxTable = useTaxTableInput,
                         dataString = paste0(dataString,"_CMS"),
                         alphaPlotWidth = 8,
                         ptOnlyFlag = TRUE)
}

cmPredHeatmapTCGA <- function(predsPT, 
                              predsBDN, 
                              plotString, 
                              midPT = 80, 
                              midBDN = 80, 
                              midAuto = FALSE,
                              ptOnlyFlag = FALSE, 
                              dimPT = 8, 
                              dimBDN = 8){

  plotPrefix <- paste0("Figures/multiclassML_",plotString,"/")
  
  # Create folder for plots if doesn't exist
  plotFolder <- paste0("Figures/multiclassML_",plotString)
  if(!( dir.exists( file.path(plotFolder)))){
    dir.create(file.path(plotFolder))
  }

  require(caret)
  require(tidyr)
  require(tibble)
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
    ggsave(paste0(plotPrefix,"PT_xgboost.jpeg"), dpi = "retina",
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
    ggsave(paste0(plotPrefix,"BDN_xgboost.jpeg"), dpi = "retina",
           width = dimBDN, height = dimBDN, units = "in")
    print(p_pt)
    ggsave(paste0(plotPrefix,"PT_xgboost_PT.jpeg"), dpi = "retina",
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

plotReadCounts <- function(metaData,
                           countData,
                           dataString,
                           plotWidths = c(4,4.5,2,6,4,5,4)){
  require(ggsci)
  
  plotPrefix <- paste0("Figures/readCounts_",dataString,"/")
  # Create folder for plots if doesn't exist
  plotFolder <- paste0("Figures/readCounts_",dataString)
  if(!( dir.exists( file.path(plotFolder)))){
    dir.create(file.path(plotFolder))
  }
  
  print("Metadata and count data are in the same order:")
  print(all(rownames(metaData) == rownames(countData)))
  metaData_Reads <- metaData %>%
    mutate(readCounts = rowSums(countData))
  
  seqCenters <- c("Baylor College of Medicine",
                  "Broad Institute of MIT and Harvard",
                  "Canada's Michael Smith Genome Sciences Centre",
                  "Harvard Medical School",
                  "MD Anderson - Institute for Applied Cancer Science",
                  "University of North Carolina",
                  "Washington University School of Medicine")
  # plotWidths <- c(4,4.5,2,6,4,5,4)
  for(ii in 1:length(seqCenters)){
    seqCenter <- seqCenters[ii]
    metaData_Reads %>%
      filter(data_submitting_center_label == seqCenter) %>%
      mutate(investigation = gsub("^TCGA-","",investigation)) %>%
      ggplot(aes(x=reorder(investigation,readCounts,median),
                 y=readCounts,fill=investigation)) +
      geom_boxplot() +
      scale_y_log10() +
      scale_fill_igv() +
      labs(x = "", y = "Filtered microbial read counts") +
      theme_pubr() +
      rotate_x_text(30) +
      stat_compare_means(label.x.npc = 0.1, size=2) +
      theme(legend.position = "none")
    ggsave(filename = paste0(plotPrefix,gsub('([[:punct:]])|\\s+','',seqCenter),
                             "_allSamples.jpeg"),
           dpi = "retina", units = "in", width = plotWidths[ii], height = 3.5)
  }
}

ancomReplot <- function(seqCenter = "BCM",
                        sampleType = "Primary Tumor",
                        dataString = "kuT2T_BIO",
                        numticks = 3,
                        fontSize = 8,
                        pointSize = 1,
                        plotWidth = 14,
                        plotHeight = 2){
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
  pGrid <- do.call("grid.arrange", c(ancomPlots, ncol=length(CTs),clip=TRUE))
  ggsave(plot = pGrid,
         filename = paste0(outputPrefix,seqCenter,"_",sampleTypeFormatted,".jpeg"),
         dpi = 900, units = "in", width = plotWidth, height = plotHeight)
}

alphaReplot <- function(seqCenter = "BCM",
                        sampleType = "Primary Tumor",
                        dataString = "kuT2T_BIO",
                        alphaDivType = "Observed",
                        plotWidth = 8,
                        plotHeight = 3,
                        fontSize = 8,
                        ptOnlyFlag = FALSE){
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
    column_to_rownames("CT") %>%
    mutate(plotColors = pal_igv("default")(length(abbrev)))
  
  plotColors <- as.character(abbreviationsTCGA_AllcancerMod$plotColors)
  names(plotColors) <- as.character(abbreviationsTCGA_AllcancerMod$abbrev)
  
  # ancomPlots <- list()
  sampleTypeFormatted <- ifelse(sampleType=="Primary Tumor","PT","BDN")
  filePath <- paste0("Figures/alphaBeta_",dataString,"_",seqCenter,"/")
  fileList <- list.files(filePath)
  fileListFilt <- fileList[grepl("\\.RData$",fileList)]
  if(length(fileListFilt)>1){
    print(fileListFilt)
    idx <- as.integer(readline("Which file to use? (select #): "))
    fileListFilt <- fileListFilt[idx]
  }
  load(paste0(filePath,fileListFilt), verbose = TRUE) # loads res object
  
  if(ptOnlyFlag){
    print("PT only!")
    alphaDivPT <- res$alphaDivPT
    
    fileName <- gsub("\\.RData","",gsub("alphaBeta_Objects_","",fileListFilt))
    
    alphaDivPT %>%
      ggplot(aes(x=reorder(investigation, !!sym(alphaDivType), median),
                 y=Observed,fill=investigation)) +
      geom_boxplot(outlier.size = 1) +
      scale_fill_manual(values = plotColors) +
      theme_pubr() +
      rotate_x_text(30) +
      labs(x="Cancer types",y=alphaDivType, title = paste(seqCenter, "PT")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position = "none",
            plot.title = element_text(face="bold"),
            text = element_text(size=fontSize)) +
      stat_compare_means(label.x.npc = 0.05,
                         label.y.npc = 0.95,
                         size = 3) -> alphaPlotPT
    ggsave(plot = alphaPlotPT,
           filename = paste0(filePath,"mergedAlphaDivPlot","_",fileName,".jpeg"),
           dpi = 900, units = "in", width = plotWidth, height = plotHeight)
  } else{
    alphaDivPT <- res$alphaDivPT
    alphaDivBDN <- res$alphaDivBDN
    
    fileName <- gsub("\\.RData","",gsub("alphaBeta_Objects_","",fileListFilt))
    
    alphaDivPT %>%
      ggplot(aes(x=reorder(investigation, !!sym(alphaDivType), median),
                 y=Observed,fill=investigation)) +
      geom_boxplot(outlier.size = 1) +
      scale_fill_manual(values = plotColors) +
      theme_pubr() +
      rotate_x_text(30) +
      labs(x="Cancer types",y=alphaDivType, title = paste(seqCenter, "PT")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position = "none",
            plot.title = element_text(face="bold"),
            text = element_text(size=fontSize)) +
      stat_compare_means(label.x.npc = 0.05,
                         label.y.npc = 0.95,
                         size = 3) -> alphaPlotPT
    
    alphaDivBDN %>%
      ggplot(aes(x=reorder(investigation, !!sym(alphaDivType), median),
                 y=Observed,fill=investigation)) +
      geom_boxplot(outlier.size = 1) +
      scale_fill_manual(values = plotColors) +
      theme_pubr() +
      rotate_x_text(30) +
      labs(x="Cancer types",y=alphaDivType, title = paste(seqCenter,"BDN")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position = "none",
            plot.title = element_text(face="bold"),
            text = element_text(size=fontSize)) +
      stat_compare_means(label.x.npc = 0.05,
                         label.y.npc = 0.95,
                         size = 3) -> alphaPlotBDN
    
    pGrid <- grid.arrange(alphaPlotPT,
                          alphaPlotBDN, ncol=2, clip=TRUE)
    ggsave(plot = pGrid,
           filename = paste0(filePath,"mergedAlphaDivPlot","_",fileName,".jpeg"),
           dpi = 900, units = "in", width = plotWidth, height = plotHeight)
  }
}