# 02-KU-T2T-processing.R
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
# Import TCGA data, subset to BIO features, and subset to nonzero samples
# Loaded in 04-read-depth-and-zebra-plots.R
#----------------------------------------------------------#

load("Interim_data/ku_hg38_T2T_Pan_NoHu_14Oct23.RData", verbose = TRUE)
# metaKUHG38
# metaKUT2T
# metaKUPan
# kuHG38Final
# kuT2TFinal
# kuPanFinal

load("Interim_data/KU_T2T_BIO_data_for_ml_tcga_by_seq_center_14Oct23.RData", verbose=TRUE)

#-------------------HG38-------------------#
kuHG38Final_BIO <- kuHG38Final[,intersect(colnames(kuHG38Final),
                                          colnames(kuT2TFinalNonzero_BIO_HiSeq_WGS))]

metaKUHG38FinalNonzero_BIO <- droplevels(metaKUHG38[which(rowSums(kuHG38Final_BIO)!=0),])
kuHG38FinalNonzero_BIO <- kuHG38Final_BIO[rownames(metaKUHG38FinalNonzero_BIO),]
dim(kuHG38FinalNonzero_BIO) # 8044   294

#-------------------T2T-------------------#
kuT2TFinal_BIO <- kuT2TFinal[,intersect(colnames(kuT2TFinal),
                                          colnames(kuT2TFinalNonzero_BIO_HiSeq_WGS))]

metaKUT2TFinalNonzero_BIO <- droplevels(metaKUT2T[which(rowSums(kuT2TFinal_BIO)!=0),])
kuT2TFinalNonzero_BIO <- kuT2TFinal_BIO[rownames(metaKUT2TFinalNonzero_BIO),]
dim(kuT2TFinalNonzero_BIO) # 7827   294

#-------------------Pan-------------------#
kuPanFinal_BIO <- kuPanFinal[,intersect(colnames(kuPanFinal),
                                          colnames(kuT2TFinalNonzero_BIO_HiSeq_WGS))]

metaKUPanFinalNonzero_BIO <- droplevels(metaKUPan[which(rowSums(kuPanFinal_BIO)!=0),])
kuPanFinalNonzero_BIO <- kuPanFinal_BIO[rownames(metaKUPanFinalNonzero_BIO),]
dim(kuPanFinalNonzero_BIO) # 7259   294

#-------------------Find common samples-------------------#

sampleInt <- Reduce(intersect,list(rownames(metaKUHG38FinalNonzero_BIO),
                                   rownames(metaKUT2TFinalNonzero_BIO),
                                   rownames(metaKUPanFinalNonzero_BIO)))
length(sampleInt) # 7259

featureInt <- Reduce(intersect,list(colnames(kuHG38Final_BIO),
                                    colnames(kuT2TFinal_BIO),
                                    colnames(kuPanFinal_BIO)))

# Final metadata
metaKUAllFinalNonzero_BIO <- metaKUPanFinalNonzero_BIO[sampleInt,]

# Final count data
kuAllFinalNonzero_BIO_HG38 <- kuHG38FinalNonzero_BIO[rownames(metaKUAllFinalNonzero_BIO),featureInt]
kuAllFinalNonzero_BIO_T2T <- kuT2TFinalNonzero_BIO[rownames(metaKUAllFinalNonzero_BIO),featureInt]
kuAllFinalNonzero_BIO_Pan <- kuPanFinalNonzero_BIO[rownames(metaKUAllFinalNonzero_BIO),featureInt]

# Save for supp tables
kuAllFinalNonzero_BIO_HG38 %>%
  rownames_to_column("sampleid") %>%
  write.csv("Supp_tables/kuAllFinalNonzero_BIO_HG38.csv",
            row.names = FALSE)
kuAllFinalNonzero_BIO_T2T %>%
  rownames_to_column("sampleid") %>%
  write.csv("Supp_tables/kuAllFinalNonzero_BIO_T2T.csv",
            row.names = FALSE)
kuAllFinalNonzero_BIO_Pan %>%
  rownames_to_column("sampleid") %>%
  write.csv("Supp_tables/kuAllFinalNonzero_BIO_Pan.csv",
            row.names = FALSE)
metaKUAllFinalNonzero_BIO %>%
  rownames_to_column("sampleid") %>%
  write.csv("Supp_tables/metaKUAllFinalNonzero_BIO_hg38_T2T_Pan.csv",
            row.names = FALSE)

#----------------------------------------------------------#
# Subset to Illumina HiSeq
#----------------------------------------------------------#

metaKUAllFinalNonzero_BIO %>% count(data_submitting_center_label, experimental_strategy)
metaKUAllFinalNonzero_BIO %>% count(cgc_platform) # HiSeq accounts for 94.2% of samples (6836 / 7259)
metaKUAllFinalNonzero_BIO %>% count(cgc_platform, data_submitting_center_label)
metaKUAllFinalNonzero_BIO %>% count(portion_is_ffpe, data_submitting_center_label)

# Subset metadata to Illumina HiSeq samples and remove the only 2 Hopkins samples
metaKUAllFinalNonzero_BIO_HiSeq <- metaKUAllFinalNonzero_BIO %>%
  filter(!is.na(data_submitting_center_label)) %>%
  filter(portion_is_ffpe == "NO") %>%
  filter(!(data_submitting_center_label %in% c("Johns Hopkins / University of Southern California",
                                               "MD Anderson - RPPA Core Facility (Proteomics)"))) %>%
  filter(cgc_platform == "Illumina HiSeq") %>% droplevels()

# Subset count data to Illumina HiSeq samples
kuAllFinalNonzero_BIO_HG38_HiSeq <- kuAllFinalNonzero_BIO_HG38[rownames(metaKUAllFinalNonzero_BIO_HiSeq),]
kuAllFinalNonzero_BIO_T2T_HiSeq <- kuAllFinalNonzero_BIO_T2T[rownames(metaKUAllFinalNonzero_BIO_HiSeq),]
kuAllFinalNonzero_BIO_Pan_HiSeq <- kuAllFinalNonzero_BIO_Pan[rownames(metaKUAllFinalNonzero_BIO_HiSeq),]

dim(kuAllFinalNonzero_BIO_HG38_HiSeq) # 6795  294
dim(kuAllFinalNonzero_BIO_T2T_HiSeq) # 6795  294
dim(kuAllFinalNonzero_BIO_Pan_HiSeq) # 6795  294

#----------------------------------------------------------#
# Subset data to individual seq centers
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUAllFinalNonzero_BIO_HiSeq_HMS <- metaKUAllFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaKUAllFinalNonzero_BIO_HiSeq_BCM <- metaKUAllFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaKUAllFinalNonzero_BIO_HiSeq_MDA <- metaKUAllFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaKUAllFinalNonzero_BIO_HiSeq_WashU <- metaKUAllFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaKUAllFinalNonzero_BIO_HiSeq_Broad_WGS <- metaKUAllFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metaKUAllFinalNonzero_BIO_HiSeq_UNC <- metaKUAllFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaKUAllFinalNonzero_BIO_HiSeq_CMS <- metaKUAllFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset HG38 count data by seq center--------------------#
# WGS per seqcenter (note that Broad has both WGS and RNA, so separate objects are made for both)
kuAllFinalNonzero_BIO_HG38_HiSeq_HMS <- kuAllFinalNonzero_BIO_HG38_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_HMS),]
kuAllFinalNonzero_BIO_HG38_HiSeq_BCM <- kuAllFinalNonzero_BIO_HG38_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_BCM),]
kuAllFinalNonzero_BIO_HG38_HiSeq_MDA <- kuAllFinalNonzero_BIO_HG38_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_MDA),]
kuAllFinalNonzero_BIO_HG38_HiSeq_WashU <- kuAllFinalNonzero_BIO_HG38_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_WashU),]
kuAllFinalNonzero_BIO_HG38_HiSeq_Broad_WGS <- kuAllFinalNonzero_BIO_HG38_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_Broad_WGS),]

# RNA-Seq per seqcenter (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
kuAllFinalNonzero_BIO_HG38_HiSeq_UNC <- kuAllFinalNonzero_BIO_HG38_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_UNC),]
kuAllFinalNonzero_BIO_HG38_HiSeq_CMS <- kuAllFinalNonzero_BIO_HG38_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_CMS),]

# WGS all
metaKUAllFinalNonzero_BIO_HiSeq_WGS <- metaKUAllFinalNonzero_BIO_HiSeq %>% filter(experimental_strategy == "WGS") %>% droplevels()
kuAllFinalNonzero_BIO_HG38_HiSeq_WGS <- kuAllFinalNonzero_BIO_HG38_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_WGS),]

# RNA all
metaKUAllFinalNonzero_BIO_HiSeq_RNA <- metaKUAllFinalNonzero_BIO_HiSeq %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
kuAllFinalNonzero_BIO_HG38_HiSeq_RNA <- kuAllFinalNonzero_BIO_HG38_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_RNA),]

#--------------------Subset T2T count data by seq center--------------------#
# WGS per seqcenter (note that Broad has both WGS and RNA, so separate objects are made for both)
kuAllFinalNonzero_BIO_T2T_HiSeq_HMS <- kuAllFinalNonzero_BIO_T2T_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_HMS),]
kuAllFinalNonzero_BIO_T2T_HiSeq_BCM <- kuAllFinalNonzero_BIO_T2T_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_BCM),]
kuAllFinalNonzero_BIO_T2T_HiSeq_MDA <- kuAllFinalNonzero_BIO_T2T_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_MDA),]
kuAllFinalNonzero_BIO_T2T_HiSeq_WashU <- kuAllFinalNonzero_BIO_T2T_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_WashU),]
kuAllFinalNonzero_BIO_T2T_HiSeq_Broad_WGS <- kuAllFinalNonzero_BIO_T2T_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_Broad_WGS),]

# RNA-Seq per seqcenter (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
kuAllFinalNonzero_BIO_T2T_HiSeq_UNC <- kuAllFinalNonzero_BIO_T2T_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_UNC),]
kuAllFinalNonzero_BIO_T2T_HiSeq_CMS <- kuAllFinalNonzero_BIO_T2T_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_CMS),]

# WGS all
metaKUAllFinalNonzero_BIO_HiSeq_WGS <- metaKUAllFinalNonzero_BIO_HiSeq %>% filter(experimental_strategy == "WGS") %>% droplevels()
kuAllFinalNonzero_BIO_T2T_HiSeq_WGS <- kuAllFinalNonzero_BIO_T2T_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_WGS),]

# RNA all
metaKUAllFinalNonzero_BIO_HiSeq_RNA <- metaKUAllFinalNonzero_BIO_HiSeq %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
kuAllFinalNonzero_BIO_T2T_HiSeq_RNA <- kuAllFinalNonzero_BIO_T2T_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_RNA),]

#--------------------Subset Pan count data by seq center--------------------#
# WGS per seqcenter (note that Broad has both WGS and RNA, so separate objects are made for both)
kuAllFinalNonzero_BIO_Pan_HiSeq_HMS <- kuAllFinalNonzero_BIO_Pan_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_HMS),]
kuAllFinalNonzero_BIO_Pan_HiSeq_BCM <- kuAllFinalNonzero_BIO_Pan_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_BCM),]
kuAllFinalNonzero_BIO_Pan_HiSeq_MDA <- kuAllFinalNonzero_BIO_Pan_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_MDA),]
kuAllFinalNonzero_BIO_Pan_HiSeq_WashU <- kuAllFinalNonzero_BIO_Pan_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_WashU),]
kuAllFinalNonzero_BIO_Pan_HiSeq_Broad_WGS <- kuAllFinalNonzero_BIO_Pan_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_Broad_WGS),]

# RNA-Seq per seqcenter (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
kuAllFinalNonzero_BIO_Pan_HiSeq_UNC <- kuAllFinalNonzero_BIO_Pan_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_UNC),]
kuAllFinalNonzero_BIO_Pan_HiSeq_CMS <- kuAllFinalNonzero_BIO_Pan_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_CMS),]

# WGS all
metaKUAllFinalNonzero_BIO_HiSeq_WGS <- metaKUAllFinalNonzero_BIO_HiSeq %>% filter(experimental_strategy == "WGS") %>% droplevels()
kuAllFinalNonzero_BIO_Pan_HiSeq_WGS <- kuAllFinalNonzero_BIO_Pan_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_WGS),]

# RNA all
metaKUAllFinalNonzero_BIO_HiSeq_RNA <- metaKUAllFinalNonzero_BIO_HiSeq %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
kuAllFinalNonzero_BIO_Pan_HiSeq_RNA <- kuAllFinalNonzero_BIO_Pan_HiSeq[rownames(metaKUAllFinalNonzero_BIO_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset HG38 count data
  kuAllFinalNonzero_BIO_HG38_HiSeq_HMS,
  kuAllFinalNonzero_BIO_HG38_HiSeq_BCM,
  kuAllFinalNonzero_BIO_HG38_HiSeq_MDA,
  kuAllFinalNonzero_BIO_HG38_HiSeq_WashU,
  kuAllFinalNonzero_BIO_HG38_HiSeq_Broad_WGS,
  kuAllFinalNonzero_BIO_HG38_HiSeq_UNC,
  kuAllFinalNonzero_BIO_HG38_HiSeq_CMS,
  kuAllFinalNonzero_BIO_HG38_HiSeq_WGS,
  kuAllFinalNonzero_BIO_HG38_HiSeq_RNA,
  
  # Subset T2T count data
  kuAllFinalNonzero_BIO_T2T_HiSeq_HMS,
  kuAllFinalNonzero_BIO_T2T_HiSeq_BCM,
  kuAllFinalNonzero_BIO_T2T_HiSeq_MDA,
  kuAllFinalNonzero_BIO_T2T_HiSeq_WashU,
  kuAllFinalNonzero_BIO_T2T_HiSeq_Broad_WGS,
  kuAllFinalNonzero_BIO_T2T_HiSeq_UNC,
  kuAllFinalNonzero_BIO_T2T_HiSeq_CMS,
  kuAllFinalNonzero_BIO_T2T_HiSeq_WGS,
  kuAllFinalNonzero_BIO_T2T_HiSeq_RNA,
  
  # Subset Pan count data
  kuAllFinalNonzero_BIO_Pan_HiSeq_HMS,
  kuAllFinalNonzero_BIO_Pan_HiSeq_BCM,
  kuAllFinalNonzero_BIO_Pan_HiSeq_MDA,
  kuAllFinalNonzero_BIO_Pan_HiSeq_WashU,
  kuAllFinalNonzero_BIO_Pan_HiSeq_Broad_WGS,
  kuAllFinalNonzero_BIO_Pan_HiSeq_UNC,
  kuAllFinalNonzero_BIO_Pan_HiSeq_CMS,
  kuAllFinalNonzero_BIO_Pan_HiSeq_WGS,
  kuAllFinalNonzero_BIO_Pan_HiSeq_RNA,
  
  # Subset metadata
  metaKUAllFinalNonzero_BIO_HiSeq_HMS,
  metaKUAllFinalNonzero_BIO_HiSeq_BCM,
  metaKUAllFinalNonzero_BIO_HiSeq_MDA,
  metaKUAllFinalNonzero_BIO_HiSeq_WashU,
  metaKUAllFinalNonzero_BIO_HiSeq_Broad_WGS,
  metaKUAllFinalNonzero_BIO_HiSeq_UNC,
  metaKUAllFinalNonzero_BIO_HiSeq_CMS,
  metaKUAllFinalNonzero_BIO_HiSeq_WGS,
  metaKUAllFinalNonzero_BIO_HiSeq_RNA,
  
  # Full metadata
  metaKUAllFinalNonzero_BIO_HiSeq,
  file = "Interim_data/KU_All_BIO_data_for_ml_tcga_by_seq_center_15Oct23.RData")

