# 02.1-wis-feature-subsetting.R
# Author: Greg Poore
# Date: July 24, 2022
# Purpose: Intersect Kraken TCGA data with WIS taxa

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
# Import TCGA data
#----------------------------------------------------------#

load("Input_data/tcgaVbDataAndMetadataAndSNM_consolidated_Nov23.RData", verbose = TRUE)

#----------------------------------------------------------#
# Import WIS data summarized at genus level
#----------------------------------------------------------#
load("Input_data/wis-bacteria-fungi-genera-species-bio-24July22.RData", verbose = TRUE)

wisGenera <- data.frame(tax_table(psWzBacteriaAndFungi_genus_Bio))
wisGeneraBact <- wisGenera %>% filter(kingdom == "Bacteria") %>% droplevels()
wisGeneraBactKnown <- wisGeneraBact %>% filter(!grepl("Unknown",genus)) %>% droplevels()
length(unique(wisGeneraBactKnown$genus)) # 216

wisGeneraBactKnownUnique <- unique(wisGeneraBactKnown$genus)

#----------------------------------------------------------#
# Intersect Kraken-derived TCGA data and WIS genera
#----------------------------------------------------------#

krakenGenera <- gsub("k__.+\\.g__","",colnames(vbDataBarnDFReconciled))

# Testing special cases
grep("Shigella",krakenGenera, value = TRUE)
grep("Escherichia",krakenGenera, value = TRUE)
grep("Clostridium",krakenGenera, value = TRUE)

# Notes on merging taxa @ genus level:
# - WIS combined "Escherichia/Shigella" and both are found in krakenGenera --> allow both
# - WIS has several versions of Clostridium but Kraken only has one --> allow one

wisGeneraBactKnownUniqueRev <- c(wisGeneraBactKnownUnique,
                                 "Escherichia",
                                 "Shigella")

# save(wisGeneraBactKnownUniqueRev,
#      file = "Interim_data/wisGeneraBactKnownUniqueRev_21Sept23.RData")

# Check overlap number
sum(krakenGenera %in% wisGeneraBactKnownUniqueRev) # 184
sum(wisGeneraBactKnownUniqueRev %in% krakenGenera) # 184

# Intersect features
krakenGeneraWISInteresected <- intersect(krakenGenera, wisGeneraBactKnownUniqueRev)

# Create new data frames of intersected dataframes
tcgaGenusKrakenAll <- vbDataBarnDFReconciled
colnames(tcgaGenusKrakenAll) <- gsub("k__.+\\.g__","",colnames(tcgaGenusKrakenAll))
tcgaGenusKrakenQC <- vbDataBarnDFReconciledQC
colnames(tcgaGenusKrakenQC) <- gsub("k__.+\\.g__","",colnames(tcgaGenusKrakenQC))

tcgaGenusKrakenAllFiltWIS <- tcgaGenusKrakenAll[, krakenGeneraWISInteresected]
tcgaGenusKrakenQCFiltWIS <- tcgaGenusKrakenQC[, krakenGeneraWISInteresected]

dim(tcgaGenusKrakenAllFiltWIS) # 18116 184
dim(tcgaGenusKrakenQCFiltWIS) # 17625 184

save(tcgaGenusKrakenAllFiltWIS,
     tcgaGenusKrakenQCFiltWIS,
     metadataSamplesAll,
     metadataSamplesAllQC,
     file = "Interim_data/tcga-wis-overlapping-data-and-metadata-subset-24July22.RData")

# Other notes:
# - The QC data will be primarily used
# - No samples have zero counts after the feature subsetting

sum(rowSums(tcgaGenusKrakenAllFiltWIS)==0) # 0
sum(rowSums(tcgaGenusKrakenQCFiltWIS)==0) # 0

#--------------------------------------------------------------------------------------------------------------------#
# Separate raw data into seq center-experimental strategy groups (to preclude needing batch correction)
#--------------------------------------------------------------------------------------------------------------------#
metadataSamplesAllQC %>% count(data_submitting_center_label, experimental_strategy)
metadataSamplesAllQC %>% count(platform) # HiSeq accounts for 91.27% of samples
metadataSamplesAllQC %>% count(platform, data_submitting_center_label)

# Subset metadata to Illumina HiSeq samples
metadataSamplesAll_HiSeq <- metadataSamplesAll %>%
  filter(platform == "Illumina HiSeq") %>% droplevels()
metadataSamplesAllQC_HiSeq <- metadataSamplesAllQC %>%
  filter(platform == "Illumina HiSeq") %>% droplevels()

# Subset count data to Illumina HiSeq samples
tcgaGenusKrakenAllFiltWIS_HiSeq <- tcgaGenusKrakenAllFiltWIS[rownames(metadataSamplesAll_HiSeq),]
tcgaGenusKrakenQCFiltWIS_HiSeq <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq),]

dim(tcgaGenusKrakenAllFiltWIS_HiSeq) # 16243 184
dim(tcgaGenusKrakenQCFiltWIS_HiSeq) # 16087 184

save(tcgaGenusKrakenAllFiltWIS_HiSeq,
     tcgaGenusKrakenQCFiltWIS_HiSeq,
     metadataSamplesAll_HiSeq,
     metadataSamplesAllQC_HiSeq,
     file = "Interim_data/tcga-hiseq-wis-overlapping-data-and-metadata-subset-25July22.RData")

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_HMS <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
metadataSamplesAllQC_HiSeq_BCM <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Baylor College of Medicine") %>% droplevels()
metadataSamplesAllQC_HiSeq_MDA <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% droplevels()
metadataSamplesAllQC_HiSeq_WashU <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Washington University School of Medicine") %>% droplevels()
metadataSamplesAllQC_HiSeq_Broad_WGS <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_UNC <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "University of North Carolina") %>% droplevels()
metadataSamplesAllQC_HiSeq_CMS <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% droplevels()
metadataSamplesAllQC_HiSeq_Broad_RNA <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
tcgaGenusKrakenQCFiltWIS_HMS <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_HMS),]
tcgaGenusKrakenQCFiltWIS_BCM <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_BCM),]
tcgaGenusKrakenQCFiltWIS_MDA <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_MDA),]
tcgaGenusKrakenQCFiltWIS_WashU <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_WashU),]
tcgaGenusKrakenQCFiltWIS_Broad_WGS <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
tcgaGenusKrakenQCFiltWIS_UNC <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_UNC),]
tcgaGenusKrakenQCFiltWIS_CMS <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_CMS),]
tcgaGenusKrakenQCFiltWIS_Broad_RNA <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------Subset metadata and count data by WGS (for multi-class classification)--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_WGS <- metadataSamplesAllQC_HiSeq %>% filter(experimental_strategy == "WGS") %>% droplevels()
tcgaGenusKrakenQCFiltWIS_WGS <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_WGS),]

# Save data for ML
save(# Subset raw count data
  tcgaGenusKrakenQCFiltWIS_WGS,
  # Subset metadata
  metadataSamplesAllQC_HiSeq_WGS,
  file = "Interim_data/data_for_multiclass_ml_tcga_wgs24July22.RData")

# Scripts: SXXX

#----------------------------------------------------#
# Save data for ML
#----------------------------------------------------#

save(# Subset raw count data
  tcgaGenusKrakenQCFiltWIS_HMS,
  tcgaGenusKrakenQCFiltWIS_BCM,
  tcgaGenusKrakenQCFiltWIS_MDA,
  tcgaGenusKrakenQCFiltWIS_WashU,
  tcgaGenusKrakenQCFiltWIS_Broad_WGS,
  tcgaGenusKrakenQCFiltWIS_UNC,
  tcgaGenusKrakenQCFiltWIS_CMS,
  tcgaGenusKrakenQCFiltWIS_Broad_RNA,
  
  # Subset metadata
  metadataSamplesAllQC_HiSeq_HMS,
  metadataSamplesAllQC_HiSeq_BCM,
  metadataSamplesAllQC_HiSeq_MDA,
  metadataSamplesAllQC_HiSeq_WashU,
  metadataSamplesAllQC_HiSeq_Broad_WGS,
  metadataSamplesAllQC_HiSeq_UNC,
  metadataSamplesAllQC_HiSeq_CMS,
  metadataSamplesAllQC_HiSeq_Broad_RNA,
  file = "Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy25July22.RData")

# Scripts: SXXX

