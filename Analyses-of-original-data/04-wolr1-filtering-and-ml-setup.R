#-----------------------------------------------------------------------------
# 04-wolr1-filtering-and-ml-setup.R
# Copyright (c) 2023--, Greg Poore
# Purposes:
# - Filter WoLr1 data based on Conterminator results
#-----------------------------------------------------------------------------

#----------------------------------------------------------#
# Load environments
#----------------------------------------------------------#

# Load dependencies
require(doMC)
require(plyr)
require(dplyr)
require(tidyr)
require(reshape2)
require(ggpubr)
require(ggsci)
require(ggrepel)
require(ggpmisc)

numCores <- detectCores()
registerDoMC(cores=numCores)

abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

#----------------------------------------------------------#
# Load data
#----------------------------------------------------------#
# Load WoLr1 data from original paper
load("Supporting_data/Shogun_data/Shogun_All_Quantile_Nov202019_Final.RData", verbose = TRUE)
# shogunMetadataQCFilt
# mergedShogunDataQCFilt
# snmDataShogunQCFilt

dim(mergedShogunDataQCFilt) # 13517  1594

# Load WoLr1 ranks
wolr1Ranks <- read.csv("Supporting_data/WoLr1-filtering/wolr1-ranks.tsv",
                       sep = "\t", stringsAsFactors = FALSE, row.names = 1)
dim(wolr1Ranks) # 10575     7

# Load WoLr1 OGUs with Conterminator hits
wolr1Conterminator <- read.csv("Supporting_data/WoLr1-filtering/wolr1-ogus-conterminator.csv",
                               stringsAsFactors = FALSE) %>% pull(OGUs)

# Subset WoLr1 ranks by Conterminator hits
sum(wolr1Conterminator %in% rownames(wolr1Ranks)) # 114
sum(rownames(wolr1Ranks) %in% wolr1Conterminator) # 114 --> 1.08% of total OGUs

wolr1ConterminatorRanks <- wolr1Ranks[wolr1Conterminator,]
dim(wolr1ConterminatorRanks) # 114 7

# Subset to affected genera
wolr1ConterminatorRanksUniqueGenera <- sort(unique(wolr1ConterminatorRanks$genus[nzchar(wolr1ConterminatorRanks$genus)]))
length(wolr1ConterminatorRanksUniqueGenera) # 84

# Subset to affected species
wolr1ConterminatorRanksUniqueSpecies <- sort(unique(wolr1ConterminatorRanks$species[nzchar(wolr1ConterminatorRanks$species)]))
length(wolr1ConterminatorRanksUniqueSpecies) # 106

## Save as supp table
require(tibble)
wolr1ConterminatorRanks %>%
  rownames_to_column("genomeID") %>%
  write.csv("Supp_Tables/WoLr1DB_conterminator_1Nov23.csv",
            row.names = FALSE)

#----------------------------------------------------------#
# Identify and remove features in Shogun **genus** data
#----------------------------------------------------------#
shogunFeats <- colnames(mergedShogunDataQCFilt)
head(shogunFeats,30)
sum(grepl("\\.g__",shogunFeats)) # 1326 --> 83.19%
sum(grepl("\\.g__[A-z]+",shogunFeats)) # 1240 --> 77.79%

# Subset to known genus level hits
mergedShogunDataQCFilt_Genus <- mergedShogunDataQCFilt
shogunFeats_Genus <- colnames(mergedShogunDataQCFilt_Genus)
shogunFeats_GenusForm <- gsub("^k__.+","",gsub("^k__.+\\.g__","",shogunFeats_Genus))

wolr1ConterminatorRanksUniqueGeneraForm <- gsub(" ",".",wolr1ConterminatorRanksUniqueGenera)
genusMatchesIdx <- unique(grep(paste(wolr1ConterminatorRanksUniqueGeneraForm,collapse="|"), 
                        shogunFeats_GenusForm, value=FALSE))

mergedShogunDataQCFilt_Genus_NoHu <- mergedShogunDataQCFilt_Genus[,-c(genusMatchesIdx)]
dim(mergedShogunDataQCFilt_Genus) # 13517  1240 | 1594
dim(mergedShogunDataQCFilt_Genus_NoHu) # 13517  1158 | 1512 --> 82 dropped (6.61% of 1240, 5.14% of 1594)

save(mergedShogunDataQCFilt_Genus,
     mergedShogunDataQCFilt_Genus_NoHu,
     shogunMetadataQCFilt,
     file = "Interim_data/merged_shogun_data_genus_level_and_Conterminator_filtered_18Sept23.RData")

#----------------------------------------------------------#
# Find overlap with WIS genera
#----------------------------------------------------------#

# Find overlap with full WIS hit list
load("Interim_data/wisGeneraBactKnownUniqueRev_21Sept23.RData", verbose = TRUE)
# wisGeneraBactKnownUniqueRev
length(wisGeneraBactKnownUniqueRev) # 218
int_WIS_All_Wolr1Conterminator <- intersect(wisGeneraBactKnownUniqueRev, 
                                            wolr1ConterminatorRanksUniqueGeneraForm)
length(int_WIS_All_Wolr1Conterminator) # 26 --> 11.93% (26 / 218)

# Find overlap within intersection of WIS and Kraken data
shogunGenusWISFeats <- intersect(shogunFeats_GenusForm,
                                 wisGeneraBactKnownUniqueRev)
length(shogunGenusWISFeats) # 189

int_WIS_TCGA_Wolr1Conterminator <- intersect(shogunGenusWISFeats, 
                                                wolr1ConterminatorRanksUniqueGeneraForm)
length(int_WIS_TCGA_Wolr1Conterminator) # 26 --> 13.76% (26 / 189)

#----------------------------------------------------------#
# Subset to Illumina HiSeq
#----------------------------------------------------------#

shogunMetadataQCFilt %>% count(data_submitting_center_label, experimental_strategy)
shogunMetadataQCFilt %>% count(platform) # HiSeq accounts for 94.83% of samples
shogunMetadataQCFilt %>% count(platform, data_submitting_center_label)

# Subset metadata to Illumina HiSeq samples and remove the only 2 Hopkins samples
shogunMetadataQCFilt_HiSeq <- shogunMetadataQCFilt %>%
  filter(platform == "Illumina HiSeq") %>% droplevels()

dim(shogunMetadataQCFilt_HiSeq) # 12818    41

shogunMetadataQCFilt_HiSeq %>% count(platform, data_submitting_center_label)

## Subset count data to Illumina HiSeq samples
# Genus
mergedShogunDataQCFilt_Genus_HiSeq <- mergedShogunDataQCFilt_Genus[rownames(shogunMetadataQCFilt_HiSeq),]
mergedShogunDataQCFilt_Genus_NoHu_HiSeq <- mergedShogunDataQCFilt_Genus_NoHu[rownames(shogunMetadataQCFilt_HiSeq),]

dim(mergedShogunDataQCFilt_Genus_HiSeq) # 12818  1512
dim(mergedShogunDataQCFilt_Genus_NoHu_HiSeq) # 12818  1512

#----------------------------------------------------------#
# Subset full data to individual seq centers
# shogunMetadataQCFilt_HiSeq
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
shogunMetadataQCFilt_HiSeq_HMS <- shogunMetadataQCFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
shogunMetadataQCFilt_HiSeq_BCM <- shogunMetadataQCFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
shogunMetadataQCFilt_HiSeq_MDA <- shogunMetadataQCFilt_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
shogunMetadataQCFilt_HiSeq_WashU <- shogunMetadataQCFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
shogunMetadataQCFilt_HiSeq_Broad_WGS <- shogunMetadataQCFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
shogunMetadataQCFilt_HiSeq_UNC <- shogunMetadataQCFilt_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
shogunMetadataQCFilt_HiSeq_CMS <- shogunMetadataQCFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()
shogunMetadataQCFilt_HiSeq_Broad_RNA <- shogunMetadataQCFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

#--------------------Subset **genus** raw count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
mergedShogunDataQCFilt_Genus_HiSeq_HMS <- mergedShogunDataQCFilt_Genus_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_HMS),]
mergedShogunDataQCFilt_Genus_HiSeq_BCM <- mergedShogunDataQCFilt_Genus_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_BCM),]
mergedShogunDataQCFilt_Genus_HiSeq_MDA <- mergedShogunDataQCFilt_Genus_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_MDA),]
mergedShogunDataQCFilt_Genus_HiSeq_WashU <- mergedShogunDataQCFilt_Genus_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_WashU),]
mergedShogunDataQCFilt_Genus_HiSeq_Broad_WGS <- mergedShogunDataQCFilt_Genus_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
mergedShogunDataQCFilt_Genus_HiSeq_UNC <- mergedShogunDataQCFilt_Genus_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_UNC),]
mergedShogunDataQCFilt_Genus_HiSeq_CMS <- mergedShogunDataQCFilt_Genus_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_CMS),]
mergedShogunDataQCFilt_Genus_HiSeq_Broad_RNA <- mergedShogunDataQCFilt_Genus_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_Broad_RNA),]

#--------------------Subset **genus** filtered count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
mergedShogunDataQCFilt_Genus_NoHu_HiSeq_HMS <- mergedShogunDataQCFilt_Genus_NoHu_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_HMS),]
mergedShogunDataQCFilt_Genus_NoHu_HiSeq_BCM <- mergedShogunDataQCFilt_Genus_NoHu_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_BCM),]
mergedShogunDataQCFilt_Genus_NoHu_HiSeq_MDA <- mergedShogunDataQCFilt_Genus_NoHu_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_MDA),]
mergedShogunDataQCFilt_Genus_NoHu_HiSeq_WashU <- mergedShogunDataQCFilt_Genus_NoHu_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_WashU),]
mergedShogunDataQCFilt_Genus_NoHu_HiSeq_Broad_WGS <- mergedShogunDataQCFilt_Genus_NoHu_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
mergedShogunDataQCFilt_Genus_NoHu_HiSeq_UNC <- mergedShogunDataQCFilt_Genus_NoHu_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_UNC),]
mergedShogunDataQCFilt_Genus_NoHu_HiSeq_CMS <- mergedShogunDataQCFilt_Genus_NoHu_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_CMS),]
mergedShogunDataQCFilt_Genus_NoHu_HiSeq_Broad_RNA <- mergedShogunDataQCFilt_Genus_NoHu_HiSeq[rownames(shogunMetadataQCFilt_HiSeq_Broad_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset genus raw count data
  mergedShogunDataQCFilt_Genus_HiSeq_HMS,
  mergedShogunDataQCFilt_Genus_HiSeq_BCM,
  mergedShogunDataQCFilt_Genus_HiSeq_MDA,
  mergedShogunDataQCFilt_Genus_HiSeq_WashU,
  mergedShogunDataQCFilt_Genus_HiSeq_Broad_WGS,
  mergedShogunDataQCFilt_Genus_HiSeq_UNC,
  mergedShogunDataQCFilt_Genus_HiSeq_CMS,
  mergedShogunDataQCFilt_Genus_HiSeq_Broad_RNA,
  
  # Subset genus filtered count data
  mergedShogunDataQCFilt_Genus_NoHu_HiSeq_HMS,
  mergedShogunDataQCFilt_Genus_NoHu_HiSeq_BCM,
  mergedShogunDataQCFilt_Genus_NoHu_HiSeq_MDA,
  mergedShogunDataQCFilt_Genus_NoHu_HiSeq_WashU,
  mergedShogunDataQCFilt_Genus_NoHu_HiSeq_Broad_WGS,
  mergedShogunDataQCFilt_Genus_NoHu_HiSeq_UNC,
  mergedShogunDataQCFilt_Genus_NoHu_HiSeq_CMS,
  mergedShogunDataQCFilt_Genus_NoHu_HiSeq_Broad_RNA,
  
  # Subset metadata
  shogunMetadataQCFilt_HiSeq_HMS,
  shogunMetadataQCFilt_HiSeq_BCM,
  shogunMetadataQCFilt_HiSeq_MDA,
  shogunMetadataQCFilt_HiSeq_WashU,
  shogunMetadataQCFilt_HiSeq_Broad_WGS,
  shogunMetadataQCFilt_HiSeq_UNC,
  shogunMetadataQCFilt_HiSeq_CMS,
  shogunMetadataQCFilt_HiSeq_Broad_RNA,

  file = "Interim_data/shogun_genus_raw_and_NoHu_data_for_ml_tcga_by_seq_center_18Sept23.RData")

# Scripts: SXXX

#----------------------------------------------------------#
# Plot results
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals
source("00-functions.R") # for barplotSummaryPerf() and barplotPerf() functions

## Load data
mlPerfAll10k_Allcancer_Shogun_NoHu <- read.csv("Supporting_scripts/S17-ML-10k-tcga-Shogun-raw-vs-NoHu-seqcenter/rep_perfML_10k_tcga_Shogun_raw_NoHu_ALL_18Sept23.csv", 
                                               stringsAsFactors = FALSE) %>%
  select(-X) %>%
  filter(!grepl("Species",datasetName)) %>%
  distinct() %>%
  rename(AUROC=auroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

source("00-functions.R")
#----------------All----------------#
# Note: The results are *not* split by seqcenters, so only use summary plot FXN
barplotSummaryPerfNoHu(inputData = mlPerfAll10k_Allcancer_Shogun_NoHu,
                       seqCenterAbbrev="All",
                       sampleTypeInput = "Blood Derived Normal",
                       prefix2Remove = "mergedShogunDataQCFilt_Genus_",
                       intFlag = FALSE,
                       plotWidthSingle = 8,
                       plotWidthCombined = 8,
                       fileNameString = "Shogun_huVSnohu")

barplotSummaryPerfNoHu(inputData = mlPerfAll10k_Allcancer_Shogun_NoHu,
                       seqCenterAbbrev="All",
                       sampleTypeInput = "Primary Tumor",
                       prefix2Remove = "mergedShogunDataQCFilt_Genus_",
                       intFlag = FALSE,
                       plotWidthSingle = 12,
                       plotWidthCombined = 12,
                       outputData = TRUE,
                       fileNameString = "Shogun_huVSnohu")

barplotSummaryPerfNoHu(inputData = mlPerfAll10k_Allcancer_Shogun_NoHu,
                       seqCenterAbbrev="All",
                       sampleTypeInput = "Primary Tumor vs Solid Tissue Normal",
                       prefix2Remove = "mergedShogunDataQCFilt_Genus_",
                       intFlag = FALSE,
                       plotWidthSingle = 5,
                       plotWidthCombined = 5,
                       fileNameString = "Shogun_huVSnohu")

