# 05-KU-T2T-processing.R
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
# Import TCGA data
# Loaded in 04-read-depth-and-zebra-plots.R
#----------------------------------------------------------#

load("Interim_data/ku_hg38_T2T_Pan_NoHu_14Oct23.RData", verbose = TRUE)

# Check for any zero-sum samples and remove -- almost all are RNA-Seq
metaKUT2TZeroSum <- metaKUT2T[which(rowSums(kuT2TFinal)==0),]
table(metaKUT2TZeroSum$experimental_strategy) # WGS=2, RNA-Seq=7458

metaKUT2TFinalNonzero_Full <- droplevels(metaKUT2T[which(rowSums(kuT2TFinal)!=0),])
kuT2TFinalNonzero_Full <- kuT2TFinal[rownames(metaKUT2TFinalNonzero_Full),]
dim(metaKUT2TFinalNonzero_Full) # 8052 44

#----------------------------------------------------------#
# Import filtered TCGA data (filtered by read counts and unique k-mer minimums)
# Loaded in 04-read-depth-and-zebra-plots.R
#----------------------------------------------------------#

## Metadata (from mycobiome paper)
load("Supporting_files/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData", verbose=T)
dim(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts) # 15512 43

metaMyco <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts

# Format myco metadata and resave
metaKU <- metaMyco %>% rownames_to_column("mycoIDs")
metaKUFilenames <- gsub("\\.filtered\\.$","",metaKU$run_prefix)
rownames(metaKU) <- metaKUFilenames

## KrakenUniq data
kuT2TFilt <- read.csv("Input_data/t2t-only/krakenuniq/combined_kraken_report_tcga_sample_ids_genus_KRAKENUNIQ_FILTER_pese.tsv",
                          sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(kuT2TFilt)[1] <- "Classification"
dim(kuT2TFilt) # 2132 15513
kuT2TFilt[1:3,1:5]

# Format features
kuT2TFiltFeatures <- gsub("^k__.+g__","",kuT2TFilt$Classification)
length(kuT2TFiltFeatures) # 2126
length(unique(kuT2TFiltFeatures)) # 2126

kuT2TFilt$Classification <- kuT2TFiltFeatures
kuT2TFilt2 <- kuT2TFilt %>% column_to_rownames("Classification")
kuT2TFilt2[1:3,1:3]

kuT2TFilt2Tr <- data.frame(t(kuT2TFilt2)) %>%
  select(-Homo)
dim(kuT2TFilt2Tr) # 15512 2125
kuT2TFilt2Tr[1:3,1:3]

# Format KrakenUniq sample IDs
kuT2TFilt2TrFor <- kuT2TFilt2Tr
rownames(kuT2TFilt2TrFor) <- gsub("\\.filtered$","",rownames(kuT2TFilt2TrFor))

# Check overlap both ways
sum(metaKUFilenames %in% rownames(kuT2TFilt2TrFor)) # 15512
sum(rownames(kuT2TFilt2TrFor) %in% metaKUFilenames) # 15512

# Reorder metadata
metaKUOrdT2TFilt <- metaKU[rownames(kuT2TFilt2TrFor),]

# Use myco sample IDs for conciseness
metaKUOrdT2TFiltFor <- metaKUOrdT2TFilt %>%
  rownames_to_column("kuIDs") %>%
  column_to_rownames("mycoIDs")

kuT2TFilt2TrFor2 <- kuT2TFilt2TrFor
identical(rownames(kuT2TFilt2TrFor2), metaKUOrdT2TFiltFor$kuIDs) # TRUE
rownames(kuT2TFilt2TrFor2) <- rownames(metaKUOrdT2TFiltFor)

# Check final results
kuT2TFilt2TrFor2[1:3,1:3]
identical(rownames(metaKUOrdT2TFiltFor),rownames(kuT2TFilt2TrFor2)) # TRUE

# Save as clean tables
metaKUT2TFiltFinal <- metaKUOrdT2TFiltFor
kuT2TFiltFinal <- kuT2TFilt2TrFor2[,colSums(kuT2TFilt2TrFor2)!=0]

# Check for any zero-sum samples and remove -- almost all are RNA-Seq
metaKUT2TFiltFinalZeroSum <- metaKUT2TFiltFinal[which(rowSums(kuT2TFiltFinal)==0),]
table(metaKUT2TFiltFinalZeroSum$experimental_strategy) # WGS=1168, RNA-Seq=8381

# Subset nonzero samples and slightly clean up object names
metaKUT2TFinalNonzero_Filt <- droplevels(metaKUT2TFiltFinal[which(rowSums(kuT2TFiltFinal)!=0),])
kuT2TFinalNonzero_Filt <- kuT2TFiltFinal[rownames(metaKUT2TFinalNonzero_Filt),]
dim(metaKUT2TFinalNonzero_Filt) # 5963   44
dim(kuT2TFinalNonzero_Filt) # 5963  671

#----------------------------------------------------------#
# Import WIS data summarized at genus level
#----------------------------------------------------------#
load("Supporting_files/wis-bacteria-fungi-genera-species-bio-24July22.RData", verbose = TRUE)

wisGenera <- data.frame(tax_table(psWzBacteriaAndFungi_genus_Bio))
wisGeneraKnown <- wisGenera %>% filter(!grepl("Unknown",genus)) %>% droplevels()
length(unique(wisGeneraKnown$genus)) # 438
wisGeneraKnownUnique <- unique(wisGeneraKnown$genus)

#----------------------------------------------------------#
# Intersect Kraken-derived TCGA data and WIS genera
#----------------------------------------------------------#

kuT2TGenera <- colnames(kuT2TFinal)
kuT2TFiltGenera <- colnames(kuT2TFiltFinal)

# Testing special cases
grep("Shigella",kuT2TGenera, value = TRUE)
grep("Escherichia",kuT2TGenera, value = TRUE)
grep("Clostridium",kuT2TGenera, value = TRUE)

# Notes on merging taxa @ genus level:
# - WIS combined "Escherichia/Shigella" and both are found in kuT2TGenera --> allow both
# - WIS has several versions of Clostridium but Kraken only has one --> allow one
wisGeneraKnownUniqueRev <- c(wisGeneraKnownUnique, "Escherichia", "Shigella")

# Check overlap number
sum(kuT2TGenera %in% wisGeneraKnownUniqueRev) # 184
sum(kuT2TFiltGenera %in% wisGeneraKnownUniqueRev) # 171

# Intersect features
kuT2TGeneraWISInteresected <- intersect(kuT2TGenera, wisGeneraKnownUniqueRev)
kuT2TFiltGeneraWISInteresected <- intersect(kuT2TFiltGenera, wisGeneraKnownUniqueRev)

# Create WIS intersected dataframe
kuT2TFinalWIS <- kuT2TFinal[, kuT2TGeneraWISInteresected]
kuT2TFiltFinalWIS <- kuT2TFiltFinal[, kuT2TFiltGeneraWISInteresected]
dim(kuT2TFinalWIS) # 15512 184
dim(kuT2TFiltFinalWIS) # 15512 171

# Check for zero-sum samples
sum(rowSums(kuT2TFinal)==0) # 7460
sum(rowSums(kuT2TFinalWIS)==0) # 7965

sum(rowSums(kuT2TFiltFinal)==0) # 9549
sum(rowSums(kuT2TFiltFinalWIS)==0) # 10026

# Check how many are due to lower RNA-Seq depth
metaKUT2TWISZeroSum <- metaKUT2T[which(rowSums(kuT2TFinalWIS)==0),]
table(metaKUT2TWISZeroSum$experimental_strategy) # WGS=3, RNA-Seq=7962

all(rownames(metaKUT2T) == rownames(kuT2TFinalWIS)) # TRUE
metaKUT2TFinalNonzero_WIS <- droplevels(metaKUT2T[which(rowSums(kuT2TFinalWIS)!=0),])
kuT2TFinalNonzero_WIS <- kuT2TFinalWIS[rownames(metaKUT2TFinalNonzero_WIS),]

all(rownames(metaKUT2TFiltFinal) == rownames(kuT2TFiltFinalWIS)) # TRUE
metaKUT2TFiltFinalNonzero_WIS <- droplevels(metaKUT2TFiltFinal[which(rowSums(kuT2TFiltFinalWIS)!=0),])
kuT2TFiltFinalNonzero_WIS <- kuT2TFiltFinalWIS[rownames(metaKUT2TFiltFinalNonzero_WIS),]

sum(rowSums(kuT2TFinalNonzero_Full)==0) # 0
sum(rowSums(kuT2TFinalNonzero_WIS)==0) # 0
sum(rowSums(kuT2TFiltFinalNonzero_WIS)==0) # 0

#----------------------------------------------------------#
# Subset to taxa pipeline (hereon "BIO")
#----------------------------------------------------------#

## Identify which viruses:
# - Passed the pathogenic filter threshold in ≥1 sample
# - Had a high number of total counts
# - Are known to infect humans or bacteria that reside in/on humans (excluding C acnes)
filtViruses <- kuT2TFinalNonzero_Filt[,grep("virus",colnames(kuT2TFinalNonzero_Filt))] # 46 viruses
filtVirusesOrd <- colSums(filtViruses)[order(colSums(filtViruses))]
# Light literature search to confirm biological associations
filtVirusesOrdPassing <- c("Mastadenovirus",
                           "Alphapapillomavirus",
                           "Lymphocryptovirus",
                           "Tequatrovirus", # phage infecting E coli, Shigella, Yersinia, Enterobacteria
                           "Roseolovirus",
                           "Simplexvirus",
                           "Betapolyomavirus",
                           "Cytomegalovirus",
                           "Skunavirus", # phage of Lactococcus
                           "Nonagvirus", # phage of E coli and Salmonella
                           "Sashavirus", # phage of Salmonella 
                           "Alphatorquevirus", # PMID: 32401076
                           "Moineauvirus", # phage of Streptococcus
                           "Gammaretrovirus", # PMID: 20219088
                           "Ceduovirus", # phage of Lactococcus
                           "Sextaecvirus", # phage of Staphyloccocus
                           "Deltapolyomavirus", # eg human polyomavirus 6
                           "Eponavirus", # eg Eponavirus epona infects Faecalibacterium
                           "Salivirus",
                           "Taranisvirus", # phage of Faecalibacterium
                           "Phietavirus", # phage of Staph aureus
                           "Brigitvirus", # PMID: 37069161
                           "Betacoronavirus",
                           "Toutatisvirus", # phage of Faecalibacterium
                           "Lentivirus",
                           "Harbinvirus", # phage of Firmicutes
                           "Slopekvirus", # phage of Klebsiella
                           "Oengusvirus", # phage of Faecalibacterium
                           "Erythroparvovirus", # eg Parvovirus B19
                           "Septimatrevirus", # phage of Pseudonomas
                           "Sepunavirus", # phage of Staphylococcus
                           "Phikzvirus", # phage of Pseudomonas
                           "Tunavirus", # phage of E coli and Shigella
                           "Rhadinovirus", # eg Kaposi sarcoma
                           "Elvirus", # phage of Pseudomonas
                           "Chivirus", # phage of Salmonella
                           "Jiaodavirus", # phage of Klebsiella
                           "Webervirus", # phage of Klebsiella
                           "Vectrevirus", # phage of Escherichia
                           "Alphapolyomavirus", # eg Merkel cell polyomavirus
                           "Fromanvirus" # phage of Mycobacterium
                           )
filtVirusesOrdNotPassing <- filtVirusesOrd[which(!(names(filtVirusesOrd) %in% filtVirusesOrdPassing))]

## Load non-viral biological taxa and intersect with KU filtered data
load("Interim_data/taxa_filtering_pipeline_13Oct23.RData", verbose=TRUE)

sum(unique(taxRS210_ff_combTaxaSpecies$genus) %in% colnames(kuT2TFinalNonzero_Filt)) # 253
sum(colnames(kuT2TFinalNonzero_Filt) %in% unique(taxRS210_ff_combTaxaSpecies$genus)) # 253

# Combine viral and bacteria genera
bioTaxa <- c(intersect(colnames(kuT2TFinalNonzero_Filt),
                       unique(taxRS210_ff_combTaxaSpecies$genus)),
             filtVirusesOrdPassing)

kuT2TFinal_BIO <- kuT2TFinal[,bioTaxa]
dim(kuT2TFinal_BIO) # 15512   294

# Check how many are due to lower RNA-Seq depth
metaKUT2TFinalBIOZeroSum <- metaKUT2T[which(rowSums(kuT2TFinal_BIO)==0),]
table(metaKUT2TFinalBIOZeroSum$experimental_strategy) # WGS=2, RNA-Seq=7683

all(rownames(metaKUT2T) == rownames(kuT2TFinal_BIO)) # TRUE
metaKUT2TFinalNonzero_BIO <- droplevels(metaKUT2T[which(rowSums(kuT2TFinal_BIO)!=0),])
kuT2TFinalNonzero_BIO <- kuT2TFinal_BIO[rownames(metaKUT2TFinalNonzero_BIO),]
dim(kuT2TFinalNonzero_BIO) # 7827  294

##------------------Create plots of dropped samples------------------##
# Percent of dropped samples per experimental strategy
metaKUT2T %>%
  mutate(nonzeroBIO = ifelse(rowSums(kuT2TFinal_BIO)!=0, TRUE, FALSE)) %>%
  count(nonzeroBIO, experimental_strategy) %>%
  mutate(experimental_strategy = factor(experimental_strategy, levels = c("WGS","RNA-Seq"))) %>%
  mutate(nonzeroBIO = factor(nonzeroBIO, levels = c(TRUE,FALSE))) %>%
  mutate(totSamples = c(10776,4736,10776,4736),
         perc = round(100*n/totSamples,2)) %>%
  ggbarplot(x = "experimental_strategy",
            y = "perc",
            palette = c("#0072B5FF","#BC3C29FF"),
            xlab = "",
            ylab = "Percent of samples (%)",
            label = TRUE, 
            lab.pos = "out",
            lab.size = 2.5,
            ylim = c(0,110),
            position = position_dodge(0.9),
            fill = "nonzeroBIO") +
  rotate_x_text(30)
ggsave("Figures/kuT2T_BIO_barplot_nonzero_samples_WGS_RNA_18Oct.jpeg",
       dpi = "retina", units = "in", width = 1.5, height = 5)

# Relate back to original file sizes
metaKUT2T %>%
  mutate(nonzeroSum = (rowSums(kuT2TFinal_BIO)!=0)) %>%
  mutate(nonzeroSum = ifelse(nonzeroSum,"Nonzero","Zero")) %>%
  ggboxplot(x = "nonzeroSum",
            y = "bam_total_reads",
            fill = "nonzeroSum",
            palette = c("#0072B5FF","#BC3C29FF"),
            notch = TRUE,
            xlab = "",
            outlier.size = 0.5,
            ylab = "Original TCGA bam file read counts",
            legend = "none") +
  scale_y_log10() +
  rotate_x_text(30) +
  stat_compare_means(label = "p.signif", comparisons = list(c("Nonzero","Zero")))
ggsave(filename = "Figures/kuT2T_BIO_bamTotalReads_barplot_nonzero_vs_zero_samples_WGS_RNA_18Oct.jpeg",
       units = "in", width = 2, height = 5, dpi = "retina")

##------------------Read count plots------------------##

source("00-Functions.R") # for plotReadCounts()
plotReadCounts(metaData = metaKUT2TFinalNonzero_BIO,
               countData = kuT2TFinalNonzero_BIO,
               dataString = "kuT2T_BIO",
               plotWidths = c(4,4.5,2,5,4,8,4))

#----------------------------------------------------------#
# Standardized naming
#----------------------------------------------------------#

## Metadata:
# - metaKUT2TFinalNonzero_Full
# - metaKUT2TFinalNonzero_Filt
# - metaKUT2TFinalNonzero_WIS
# - metaKUT2TFinalNonzero_BIO

## Count data:
# - kuT2TFinalNonzero_Full
# - kuT2TFinalNonzero_Filt
# - kuT2TFinalNonzero_WIS
# - kuT2TFinalNonzero_BIO

## Create version of Filt that only uses BIO features
kuT2TFinalNonzero_BIOFiltTmp <- kuT2TFinalNonzero_Filt[,intersect(colnames(kuT2TFinalNonzero_Filt),
                                                               colnames(kuT2TFinalNonzero_BIO))]
kuT2TFinalNonzero_BIOFilt <- kuT2TFinalNonzero_BIOFiltTmp[rowSums(kuT2TFinalNonzero_BIOFiltTmp)!=0,]
metaKUT2TFinalNonzero_BIOFilt <- droplevels(metaKUT2TFinalNonzero_BIO[rownames(kuT2TFinalNonzero_BIOFilt),])
dim(kuT2TFinalNonzero_BIOFilt) #  5957  294 --> 6 samples were dropped from Filt
dim(metaKUT2TFinalNonzero_BIOFilt) #  5957  44

# save( # Metadata
#   metaKUT2TFinalNonzero_Full,
#   metaKUT2TFinalNonzero_Filt,
#   metaKUT2TFinalNonzero_WIS,
#   metaKUT2TFinalNonzero_BIO,
#   metaKUT2TFinalNonzero_BIOFilt,
#   # Count data
#   kuT2TFinalNonzero_Full,
#   kuT2TFinalNonzero_Filt,
#   kuT2TFinalNonzero_WIS,
#   kuT2TFinalNonzero_BIO,
#   kuT2TFinalNonzero_BIOFilt,
#   file = "Interim_data/kuMetaCounts_Full_Filt_WIS_BIO_14Oct23.RData"
# )

load("Interim_data/kuMetaCounts_Full_Filt_WIS_BIO_14Oct23.RData")

# Save for supp tables
write.csv(colnames(kuT2TFinalNonzero_BIO),
          file = "Supp_tables/kuT2T_BIO_features.csv")
kuT2TFinalNonzero_BIO %>%
  rownames_to_column("sampleid") %>%
  write.csv("Supp_tables/kuT2TFinalNonzero_BIO_Fullset.csv",
            row.names = FALSE)
metaKUT2TFinalNonzero_BIO %>%
  rownames_to_column("sampleid") %>%
  write.csv("Supp_tables/metaKUT2TFinalNonzero_BIO_Fullset.csv",
            row.names = FALSE)

#----------------------------------------------------------#
# Subset to Illumina HiSeq
#----------------------------------------------------------#

metaKUT2TFinalNonzero_Full %>% count(data_submitting_center_label, experimental_strategy)
metaKUT2TFinalNonzero_Full %>% count(cgc_platform) # HiSeq accounts for 94.76% of samples (7628 / 8052)
metaKUT2TFinalNonzero_Full %>% count(cgc_platform, data_submitting_center_label)
metaKUT2TFinalNonzero_Full %>% count(portion_is_ffpe, data_submitting_center_label)

# Subset metadata to Illumina HiSeq samples and remove the only 2 Hopkins samples; also remove MDA RPPA 11 samples
metaKUT2TFinalNonzero_Full_HiSeq <- metaKUT2TFinalNonzero_Full %>%
  filter(!is.na(data_submitting_center_label)) %>%
  filter(portion_is_ffpe == "NO") %>%
  filter(!(data_submitting_center_label %in% c("Johns Hopkins / University of Southern California",
                                               "MD Anderson - RPPA Core Facility (Proteomics)"))) %>%
  filter(cgc_platform == "Illumina HiSeq") %>% droplevels()
# WIS
metaKUT2TFinalNonzero_WIS_HiSeq <- metaKUT2TFinalNonzero_Full_HiSeq %>%
  rownames_to_column("sampleID") %>%
  filter(sampleID %in% rownames(metaKUT2TFinalNonzero_WIS)) %>%
  column_to_rownames("sampleID") %>% droplevels()
# BIO
metaKUT2TFinalNonzero_BIO_HiSeq <- metaKUT2TFinalNonzero_Full_HiSeq %>%
  rownames_to_column("sampleID") %>%
  filter(sampleID %in% rownames(metaKUT2TFinalNonzero_BIO)) %>%
  column_to_rownames("sampleID") %>% droplevels()
# Filt
metaKUT2TFinalNonzero_Filt_HiSeq <- metaKUT2TFinalNonzero_Full_HiSeq %>%
  rownames_to_column("sampleID") %>%
  filter(sampleID %in% rownames(metaKUT2TFinalNonzero_Filt)) %>%
  column_to_rownames("sampleID") %>% droplevels()
# BIOFilt
metaKUT2TFinalNonzero_BIOFilt_HiSeq <- metaKUT2TFinalNonzero_Full_HiSeq %>%
  rownames_to_column("sampleID") %>%
  filter(sampleID %in% rownames(metaKUT2TFinalNonzero_BIOFilt)) %>%
  column_to_rownames("sampleID") %>% droplevels()

# Subset count data to Illumina HiSeq samples
kuT2TFinalNonzero_Full_HiSeq <- kuT2TFinalNonzero_Full[rownames(metaKUT2TFinalNonzero_Full_HiSeq),]
kuT2TFinalNonzero_WIS_HiSeq <- kuT2TFinalNonzero_WIS[rownames(metaKUT2TFinalNonzero_WIS_HiSeq),]
kuT2TFinalNonzero_BIO_HiSeq <- kuT2TFinalNonzero_BIO[rownames(metaKUT2TFinalNonzero_BIO_HiSeq),]
kuT2TFinalNonzero_Filt_HiSeq <- kuT2TFinalNonzero_Filt[rownames(metaKUT2TFinalNonzero_Filt_HiSeq),]
kuT2TFinalNonzero_BIOFilt_HiSeq <- kuT2TFinalNonzero_BIOFilt[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq),]

dim(kuT2TFinalNonzero_Full_HiSeq) # 7586 2125
dim(kuT2TFinalNonzero_WIS_HiSeq) # 7082  184
dim(kuT2TFinalNonzero_BIO_HiSeq) # 7362  294
dim(kuT2TFinalNonzero_Filt_HiSeq) # 5549  671
dim(kuT2TFinalNonzero_BIOFilt_HiSeq) # 5544  294

#----------------------------------------------------------#
# Subset data to individual seq centers: Full data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_Full_HiSeq_HMS <- metaKUT2TFinalNonzero_Full_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaKUT2TFinalNonzero_Full_HiSeq_BCM <- metaKUT2TFinalNonzero_Full_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_Full_HiSeq_MDA <- metaKUT2TFinalNonzero_Full_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaKUT2TFinalNonzero_Full_HiSeq_WashU <- metaKUT2TFinalNonzero_Full_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_Full_HiSeq_Broad_WGS <- metaKUT2TFinalNonzero_Full_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metaKUT2TFinalNonzero_Full_HiSeq_UNC <- metaKUT2TFinalNonzero_Full_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaKUT2TFinalNonzero_Full_HiSeq_CMS <- metaKUT2TFinalNonzero_Full_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_Full_HiSeq_HMS <- kuT2TFinalNonzero_Full_HiSeq[rownames(metaKUT2TFinalNonzero_Full_HiSeq_HMS),]
kuT2TFinalNonzero_Full_HiSeq_BCM <- kuT2TFinalNonzero_Full_HiSeq[rownames(metaKUT2TFinalNonzero_Full_HiSeq_BCM),]
kuT2TFinalNonzero_Full_HiSeq_MDA <- kuT2TFinalNonzero_Full_HiSeq[rownames(metaKUT2TFinalNonzero_Full_HiSeq_MDA),]
kuT2TFinalNonzero_Full_HiSeq_WashU <- kuT2TFinalNonzero_Full_HiSeq[rownames(metaKUT2TFinalNonzero_Full_HiSeq_WashU),]
kuT2TFinalNonzero_Full_HiSeq_Broad_WGS <- kuT2TFinalNonzero_Full_HiSeq[rownames(metaKUT2TFinalNonzero_Full_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_Full_HiSeq_UNC <- kuT2TFinalNonzero_Full_HiSeq[rownames(metaKUT2TFinalNonzero_Full_HiSeq_UNC),]
kuT2TFinalNonzero_Full_HiSeq_CMS <- kuT2TFinalNonzero_Full_HiSeq[rownames(metaKUT2TFinalNonzero_Full_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_Full_HiSeq_WGS <- metaKUT2TFinalNonzero_Full_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
kuT2TFinalNonzero_Full_HiSeq_WGS <- kuT2TFinalNonzero_Full_HiSeq[rownames(metaKUT2TFinalNonzero_Full_HiSeq_WGS),]

# RNA (note that Broad has both RNA and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_Full_HiSeq_RNA <- metaKUT2TFinalNonzero_Full_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
kuT2TFinalNonzero_Full_HiSeq_RNA <- kuT2TFinalNonzero_Full_HiSeq[rownames(metaKUT2TFinalNonzero_Full_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_Full_HiSeq_HMS,
  kuT2TFinalNonzero_Full_HiSeq_BCM,
  kuT2TFinalNonzero_Full_HiSeq_MDA,
  kuT2TFinalNonzero_Full_HiSeq_WashU,
  kuT2TFinalNonzero_Full_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_Full_HiSeq_UNC,
  kuT2TFinalNonzero_Full_HiSeq_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_Full_HiSeq_WGS,
  kuT2TFinalNonzero_Full_HiSeq_RNA,
  
  # Subset metadata
  metaKUT2TFinalNonzero_Full_HiSeq_HMS,
  metaKUT2TFinalNonzero_Full_HiSeq_BCM,
  metaKUT2TFinalNonzero_Full_HiSeq_MDA,
  metaKUT2TFinalNonzero_Full_HiSeq_WashU,
  metaKUT2TFinalNonzero_Full_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_Full_HiSeq_UNC,
  metaKUT2TFinalNonzero_Full_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_Full_HiSeq_WGS,
  metaKUT2TFinalNonzero_Full_HiSeq_RNA,
  
  # Full metadata
  metaKUT2TFinalNonzero_Full_HiSeq,
  file = "Interim_data/KU_T2T_Full_data_for_ml_tcga_by_seq_center_14Oct23.RData")

#----------------------------------------------------------#
# Subset data to individual seq centers: WIS data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_WIS_HiSeq_HMS <- metaKUT2TFinalNonzero_WIS_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaKUT2TFinalNonzero_WIS_HiSeq_BCM <- metaKUT2TFinalNonzero_WIS_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_WIS_HiSeq_MDA <- metaKUT2TFinalNonzero_WIS_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaKUT2TFinalNonzero_WIS_HiSeq_WashU <- metaKUT2TFinalNonzero_WIS_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_WIS_HiSeq_Broad_WGS <- metaKUT2TFinalNonzero_WIS_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metaKUT2TFinalNonzero_WIS_HiSeq_UNC <- metaKUT2TFinalNonzero_WIS_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaKUT2TFinalNonzero_WIS_HiSeq_CMS <- metaKUT2TFinalNonzero_WIS_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_WIS_HiSeq_HMS <- kuT2TFinalNonzero_WIS_HiSeq[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_HMS),]
kuT2TFinalNonzero_WIS_HiSeq_BCM <- kuT2TFinalNonzero_WIS_HiSeq[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_BCM),]
kuT2TFinalNonzero_WIS_HiSeq_MDA <- kuT2TFinalNonzero_WIS_HiSeq[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_MDA),]
kuT2TFinalNonzero_WIS_HiSeq_WashU <- kuT2TFinalNonzero_WIS_HiSeq[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_WashU),]
kuT2TFinalNonzero_WIS_HiSeq_Broad_WGS <- kuT2TFinalNonzero_WIS_HiSeq[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_WIS_HiSeq_UNC <- kuT2TFinalNonzero_WIS_HiSeq[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_UNC),]
kuT2TFinalNonzero_WIS_HiSeq_CMS <- kuT2TFinalNonzero_WIS_HiSeq[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_WIS_HiSeq_WGS <- metaKUT2TFinalNonzero_WIS_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
kuT2TFinalNonzero_WIS_HiSeq_WGS <- kuT2TFinalNonzero_WIS_HiSeq[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_WGS),]

# RNA (note that Broad has both RNA and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_WIS_HiSeq_RNA <- metaKUT2TFinalNonzero_WIS_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
kuT2TFinalNonzero_WIS_HiSeq_RNA <- kuT2TFinalNonzero_WIS_HiSeq[rownames(metaKUT2TFinalNonzero_WIS_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_WIS_HiSeq_HMS,
  kuT2TFinalNonzero_WIS_HiSeq_BCM,
  kuT2TFinalNonzero_WIS_HiSeq_MDA,
  kuT2TFinalNonzero_WIS_HiSeq_WashU,
  kuT2TFinalNonzero_WIS_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_WIS_HiSeq_UNC,
  kuT2TFinalNonzero_WIS_HiSeq_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_WIS_HiSeq_WGS,
  kuT2TFinalNonzero_WIS_HiSeq_RNA,
  
  # Subset metadata
  metaKUT2TFinalNonzero_WIS_HiSeq_HMS,
  metaKUT2TFinalNonzero_WIS_HiSeq_BCM,
  metaKUT2TFinalNonzero_WIS_HiSeq_MDA,
  metaKUT2TFinalNonzero_WIS_HiSeq_WashU,
  metaKUT2TFinalNonzero_WIS_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_WIS_HiSeq_UNC,
  metaKUT2TFinalNonzero_WIS_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_WIS_HiSeq_WGS,
  metaKUT2TFinalNonzero_WIS_HiSeq_RNA,
  
  # WIS metadata
  metaKUT2TFinalNonzero_WIS_HiSeq,
  file = "Interim_data/KU_T2T_WIS_data_for_ml_tcga_by_seq_center_14Oct23.RData")

#----------------------------------------------------------#
# Subset data to individual seq centers: BIO data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_BIO_HiSeq_HMS <- metaKUT2TFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIO_HiSeq_BCM <- metaKUT2TFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIO_HiSeq_MDA <- metaKUT2TFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIO_HiSeq_WashU <- metaKUT2TFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIO_HiSeq_Broad_WGS <- metaKUT2TFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metaKUT2TFinalNonzero_BIO_HiSeq_UNC <- metaKUT2TFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIO_HiSeq_CMS <- metaKUT2TFinalNonzero_BIO_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_BIO_HiSeq_HMS <- kuT2TFinalNonzero_BIO_HiSeq[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_HMS),]
kuT2TFinalNonzero_BIO_HiSeq_BCM <- kuT2TFinalNonzero_BIO_HiSeq[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_BCM),]
kuT2TFinalNonzero_BIO_HiSeq_MDA <- kuT2TFinalNonzero_BIO_HiSeq[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_MDA),]
kuT2TFinalNonzero_BIO_HiSeq_WashU <- kuT2TFinalNonzero_BIO_HiSeq[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_WashU),]
kuT2TFinalNonzero_BIO_HiSeq_Broad_WGS <- kuT2TFinalNonzero_BIO_HiSeq[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_BIO_HiSeq_UNC <- kuT2TFinalNonzero_BIO_HiSeq[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_UNC),]
kuT2TFinalNonzero_BIO_HiSeq_CMS <- kuT2TFinalNonzero_BIO_HiSeq[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_BIO_HiSeq_WGS <- metaKUT2TFinalNonzero_BIO_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
kuT2TFinalNonzero_BIO_HiSeq_WGS <- kuT2TFinalNonzero_BIO_HiSeq[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_WGS),]

# RNA (note that Broad has both RNA and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_BIO_HiSeq_RNA <- metaKUT2TFinalNonzero_BIO_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
kuT2TFinalNonzero_BIO_HiSeq_RNA <- kuT2TFinalNonzero_BIO_HiSeq[rownames(metaKUT2TFinalNonzero_BIO_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_BIO_HiSeq_HMS,
  kuT2TFinalNonzero_BIO_HiSeq_BCM,
  kuT2TFinalNonzero_BIO_HiSeq_MDA,
  kuT2TFinalNonzero_BIO_HiSeq_WashU,
  kuT2TFinalNonzero_BIO_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_BIO_HiSeq_UNC,
  kuT2TFinalNonzero_BIO_HiSeq_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_BIO_HiSeq_WGS,
  kuT2TFinalNonzero_BIO_HiSeq_RNA,
  
  # Subset metadata
  metaKUT2TFinalNonzero_BIO_HiSeq_HMS,
  metaKUT2TFinalNonzero_BIO_HiSeq_BCM,
  metaKUT2TFinalNonzero_BIO_HiSeq_MDA,
  metaKUT2TFinalNonzero_BIO_HiSeq_WashU,
  metaKUT2TFinalNonzero_BIO_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_BIO_HiSeq_UNC,
  metaKUT2TFinalNonzero_BIO_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_BIO_HiSeq_WGS,
  metaKUT2TFinalNonzero_BIO_HiSeq_RNA,
  
  # BIO metadata
  metaKUT2TFinalNonzero_BIO_HiSeq,
  file = "Interim_data/KU_T2T_BIO_data_for_ml_tcga_by_seq_center_14Oct23.RData")

#----------------------------------------------------------#
# Subset data to individual seq centers: Filt data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_Filt_HiSeq_HMS <- metaKUT2TFinalNonzero_Filt_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaKUT2TFinalNonzero_Filt_HiSeq_BCM <- metaKUT2TFinalNonzero_Filt_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_Filt_HiSeq_MDA <- metaKUT2TFinalNonzero_Filt_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaKUT2TFinalNonzero_Filt_HiSeq_WashU <- metaKUT2TFinalNonzero_Filt_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_Filt_HiSeq_Broad_WGS <- metaKUT2TFinalNonzero_Filt_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metaKUT2TFinalNonzero_Filt_HiSeq_UNC <- metaKUT2TFinalNonzero_Filt_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaKUT2TFinalNonzero_Filt_HiSeq_CMS <- metaKUT2TFinalNonzero_Filt_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_Filt_HiSeq_HMS <- kuT2TFinalNonzero_Filt_HiSeq[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_HMS),]
kuT2TFinalNonzero_Filt_HiSeq_BCM <- kuT2TFinalNonzero_Filt_HiSeq[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_BCM),]
kuT2TFinalNonzero_Filt_HiSeq_MDA <- kuT2TFinalNonzero_Filt_HiSeq[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_MDA),]
kuT2TFinalNonzero_Filt_HiSeq_WashU <- kuT2TFinalNonzero_Filt_HiSeq[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_WashU),]
kuT2TFinalNonzero_Filt_HiSeq_Broad_WGS <- kuT2TFinalNonzero_Filt_HiSeq[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_Filt_HiSeq_UNC <- kuT2TFinalNonzero_Filt_HiSeq[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_UNC),]
kuT2TFinalNonzero_Filt_HiSeq_CMS <- kuT2TFinalNonzero_Filt_HiSeq[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_Filt_HiSeq_WGS <- metaKUT2TFinalNonzero_Filt_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
kuT2TFinalNonzero_Filt_HiSeq_WGS <- kuT2TFinalNonzero_Filt_HiSeq[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_WGS),]

# RNA (note that Broad has both RNA and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_Filt_HiSeq_RNA <- metaKUT2TFinalNonzero_Filt_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
kuT2TFinalNonzero_Filt_HiSeq_RNA <- kuT2TFinalNonzero_Filt_HiSeq[rownames(metaKUT2TFinalNonzero_Filt_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_Filt_HiSeq_HMS,
  kuT2TFinalNonzero_Filt_HiSeq_BCM,
  kuT2TFinalNonzero_Filt_HiSeq_MDA,
  kuT2TFinalNonzero_Filt_HiSeq_WashU,
  kuT2TFinalNonzero_Filt_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_Filt_HiSeq_UNC,
  kuT2TFinalNonzero_Filt_HiSeq_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_Filt_HiSeq_WGS,
  kuT2TFinalNonzero_Filt_HiSeq_RNA,
  
  # Subset metadata
  metaKUT2TFinalNonzero_Filt_HiSeq_HMS,
  metaKUT2TFinalNonzero_Filt_HiSeq_BCM,
  metaKUT2TFinalNonzero_Filt_HiSeq_MDA,
  metaKUT2TFinalNonzero_Filt_HiSeq_WashU,
  metaKUT2TFinalNonzero_Filt_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_Filt_HiSeq_UNC,
  metaKUT2TFinalNonzero_Filt_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_Filt_HiSeq_WGS,
  metaKUT2TFinalNonzero_Filt_HiSeq_RNA,
  
  # Filt metadata
  metaKUT2TFinalNonzero_Filt_HiSeq,
  file = "Interim_data/KU_T2T_Filt_data_for_ml_tcga_by_seq_center_14Oct23.RData")

#----------------------------------------------------------#
# Subset data to individual seq centers: BIOFilt data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_BIOFilt_HiSeq_HMS <- metaKUT2TFinalNonzero_BIOFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIOFilt_HiSeq_BCM <- metaKUT2TFinalNonzero_BIOFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIOFilt_HiSeq_MDA <- metaKUT2TFinalNonzero_BIOFilt_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIOFilt_HiSeq_WashU <- metaKUT2TFinalNonzero_BIOFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIOFilt_HiSeq_Broad_WGS <- metaKUT2TFinalNonzero_BIOFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metaKUT2TFinalNonzero_BIOFilt_HiSeq_UNC <- metaKUT2TFinalNonzero_BIOFilt_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaKUT2TFinalNonzero_BIOFilt_HiSeq_CMS <- metaKUT2TFinalNonzero_BIOFilt_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
kuT2TFinalNonzero_BIOFilt_HiSeq_HMS <- kuT2TFinalNonzero_BIOFilt_HiSeq[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_HMS),]
kuT2TFinalNonzero_BIOFilt_HiSeq_BCM <- kuT2TFinalNonzero_BIOFilt_HiSeq[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_BCM),]
kuT2TFinalNonzero_BIOFilt_HiSeq_MDA <- kuT2TFinalNonzero_BIOFilt_HiSeq[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_MDA),]
kuT2TFinalNonzero_BIOFilt_HiSeq_WashU <- kuT2TFinalNonzero_BIOFilt_HiSeq[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_WashU),]
kuT2TFinalNonzero_BIOFilt_HiSeq_Broad_WGS <- kuT2TFinalNonzero_BIOFilt_HiSeq[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
kuT2TFinalNonzero_BIOFilt_HiSeq_UNC <- kuT2TFinalNonzero_BIOFilt_HiSeq[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_UNC),]
kuT2TFinalNonzero_BIOFilt_HiSeq_CMS <- kuT2TFinalNonzero_BIOFilt_HiSeq[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_BIOFilt_HiSeq_WGS <- metaKUT2TFinalNonzero_BIOFilt_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
kuT2TFinalNonzero_BIOFilt_HiSeq_WGS <- kuT2TFinalNonzero_BIOFilt_HiSeq[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_WGS),]

# RNA (note that Broad has both RNA and RNA, so separate objects are made for both)
metaKUT2TFinalNonzero_BIOFilt_HiSeq_RNA <- metaKUT2TFinalNonzero_BIOFilt_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
kuT2TFinalNonzero_BIOFilt_HiSeq_RNA <- kuT2TFinalNonzero_BIOFilt_HiSeq[rownames(metaKUT2TFinalNonzero_BIOFilt_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  kuT2TFinalNonzero_BIOFilt_HiSeq_HMS,
  kuT2TFinalNonzero_BIOFilt_HiSeq_BCM,
  kuT2TFinalNonzero_BIOFilt_HiSeq_MDA,
  kuT2TFinalNonzero_BIOFilt_HiSeq_WashU,
  kuT2TFinalNonzero_BIOFilt_HiSeq_Broad_WGS,
  kuT2TFinalNonzero_BIOFilt_HiSeq_UNC,
  kuT2TFinalNonzero_BIOFilt_HiSeq_CMS,
  
  # WGS and RNA subset data
  kuT2TFinalNonzero_BIOFilt_HiSeq_WGS,
  kuT2TFinalNonzero_BIOFilt_HiSeq_RNA,
  
  # Subset metadata
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_HMS,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_BCM,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_MDA,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_WashU,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_Broad_WGS,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_UNC,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_WGS,
  metaKUT2TFinalNonzero_BIOFilt_HiSeq_RNA,
  
  # BIOFilt metadata
  metaKUT2TFinalNonzero_BIOFilt_HiSeq,
  file = "Interim_data/KU_T2T_BIOFilt_data_for_ml_tcga_by_seq_center_14Oct23.RData")

#----------------------------------------------------------#
# Batch correction with ConQuR
# NOTE: Done on cloud -- see supplementary scripts
#----------------------------------------------------------#
# Load data from above
load("Interim_data/KU_T2T_Full_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("Interim_data/KU_T2T_WIS_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("Interim_data/KU_T2T_BIO_data_for_ml_tcga_by_seq_center_14Oct23.RData", verbose = TRUE)
load("Interim_data/KU_T2T_Filt_data_for_ml_tcga_by_seq_center_14Oct23.RData")
load("Interim_data/KU_T2T_BIOFilt_data_for_ml_tcga_by_seq_center_14Oct23.RData")

# # Load ConQuR-corrected data
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_WIS_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIO_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Filt_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_BIOFilt_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/krakenuniq/kuT2TFinalNonzero_Full_HiSeq_RNA_CQ_CMS.RData")

#----------------------------------------------------------#
# Calculate alpha and beta diversity per batch using phyloseq
#----------------------------------------------------------#
require(ggsci)
require(ggpubr)

source("00-Functions.R") # for alphaBetaFXN() and runAlphaBetaSeqCenter()
# Done
runAlphaBetaSeqCenter(metaString = "metaKUT2TFinalNonzero_BIO_HiSeq",
                      countString = "kuT2TFinalNonzero_BIO_HiSeq",
                      dataStringInput = "kuT2T_BIO",
                      useTaxTableInput = FALSE)

#----------------------------------------------------#
# Export per-batch data to run Qiime2
#----------------------------------------------------#
require(biomformat)
require(rhdf5)

# ## Create and write taxa file
# taxaFile <- data.frame(`Feature ID` = colnames(kuT2TFinalNonzero_Full_HiSeq_WGS),
#                        Taxon = colnames(kuT2TFinalNonzero_Full_HiSeq_WGS),
#                        check.names = FALSE)
# write.table(taxaFile,
#             file = "Qiime_data_and_scripts/Qiime_input_data/ku-taxa.txt",
#             quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

source("00-Functions.R") # for export2Qiime()
export2Qiime(metaString = "metaKUT2TFinalNonzero_BIO_HiSeq",
             countString = "kuT2TFinalNonzero_BIO_HiSeq",
             dataString = "kuT2T_BIO")

#----------------------------------------------------#
# Export ConQuR-corrected data to run Qiime2
#----------------------------------------------------#

export2QiimeCQ <- function(metaString = "metaKUT2TFinalNonzero_BIO_HiSeq",
                           countString = "kuT2TFinalNonzero_BIO_HiSeq",
                           dataString = "kuT2T_BIO"){
  
  require(biomformat)
  require(rhdf5)
  
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

source("00-Functions.R") # for export2QiimeCQ()
export2QiimeCQ(metaString = "metaKUT2TFinalNonzero_BIO_HiSeq",
             countString = "kuT2TFinalNonzero_BIO_HiSeq",
             dataString = "kuT2T_BIO")

#----------------------------------------------------#
# Calculate per-center differential abundances
#----------------------------------------------------#

# # Load data from above
# load("Interim_data/KU_T2T_Full_data_for_ml_tcga_by_seq_center_10Oct23.RData")
# load("Interim_data/KU_T2T_WIS_data_for_ml_tcga_by_seq_center_10Oct23.RData")
# load("Interim_data/KU_T2T_BIO_data_for_ml_tcga_by_seq_center_10Oct23.RData")
# load("Interim_data/KU_T2T_Filt_data_for_ml_tcga_by_seq_center_10Oct23.RData")

source("00-Functions.R") # for runAncomBC_1VsAll_OGUs()
# BIO
runAncomBC_1VsAll_OGUs(metaString = "metaKUT2TFinalNonzero_BIO_HiSeq",
                       countString = "kuT2TFinalNonzero_BIO_HiSeq",
                       dataString = "kuT2T_BIO",
                       taxTable = NA,
                       makeTaxTable = TRUE,
                       qvalCutoff = 0.05,
                       sampleTypes = c("Primary Tumor","Blood Derived Normal"),
                       SeqCenters = c("HMS","BCM","MDA","WashU","Broad_WGS","UNC","CMS"),
                       taxaPlotLabel = "genus")

metaKUT2TFinalNonzero_BIO_HiSeq_WGS_CQ <- metaKUT2TFinalNonzero_BIO_HiSeq_WGS
runAncomBC_1VsAll_OGUs(metaString = "metaKUT2TFinalNonzero_BIO_HiSeq_WGS",
                       countString = "kuT2TFinalNonzero_BIO_HiSeq_WGS",
                       dataString = "kuT2T_BIO_WGS_CQ",
                       taxTable = NA,
                       makeTaxTable = TRUE,
                       qvalCutoff = 0.05,
                       showTopXFlag = FALSE,
                       sampleTypes = c("Primary Tumor","Blood Derived Normal"),
                       SeqCenters = c("CQ"),
                       taxaPlotLabel = "genus")

metaKUT2TFinalNonzero_BIO_HiSeq_RNA_CQ <- metaKUT2TFinalNonzero_BIO_HiSeq_RNA
runAncomBC_1VsAll_OGUs(metaString = "metaKUT2TFinalNonzero_BIO_HiSeq_RNA",
                       countString = "kuT2TFinalNonzero_BIO_HiSeq_RNA",
                       dataString = "kuT2T_BIO_RNA_CQ",
                       taxTable = NA,
                       makeTaxTable = TRUE,
                       qvalCutoff = 0.05,
                       showTopXFlag = FALSE,
                       sampleTypes = c("Primary Tumor"),
                       SeqCenters = c("CQ"),
                       taxaPlotLabel = "genus")
#----------------------------------------------------#
# Replot ANCOM-BC results
#----------------------------------------------------#

source("00-Functions.R") # for ancomReplot()
# BCM
ancomReplot(seqCenter = "BCM",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO")
ancomReplot(seqCenter = "BCM",
            sampleType = "Blood Derived Normal",
            dataString = "kuT2T_BIO")
# HMS
ancomReplot(seqCenter = "HMS",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            fontSize = 7,
            pointSize = 0.5)
ancomReplot(seqCenter = "HMS",
            sampleType = "Blood Derived Normal",
            dataString = "kuT2T_BIO",
            fontSize = 7,
            pointSize = 0.5)
# WashU
ancomReplot(seqCenter = "WashU",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
ancomReplot(seqCenter = "WashU",
            sampleType = "Blood Derived Normal",
            dataString = "kuT2T_BIO",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
# Broad_WGS
ancomReplot(seqCenter = "Broad_WGS",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            fontSize = 7,
            pointSize = 0.5)
ancomReplot(seqCenter = "Broad_WGS",
            sampleType = "Blood Derived Normal",
            dataString = "kuT2T_BIO",
            fontSize = 7,
            pointSize = 0.5)
# MDA
ancomReplot(seqCenter = "MDA",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
ancomReplot(seqCenter = "MDA",
            sampleType = "Blood Derived Normal",
            dataString = "kuT2T_BIO",
            numticks = 2,
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
# CMS
ancomReplot(seqCenter = "CMS",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)

#----------------------------------------------------#
# Replot alpha diversity results
#----------------------------------------------------#

source("00-Functions.R") # for alphaReplot()
# BCM
alphaReplot(seqCenter = "BCM",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            alphaDivType = "Observed",
            plotWidth = 6,
            plotHeight = 3,
            fontSize = 8)
# HMS
alphaReplot(seqCenter = "HMS",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# Broad_WGS
alphaReplot(seqCenter = "Broad_WGS",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# WashU
alphaReplot(seqCenter = "WashU",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# MDA
alphaReplot(seqCenter = "MDA",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# CMS
alphaReplot(seqCenter = "CMS",
            sampleType = "Primary Tumor",
            dataString = "kuT2T_BIO",
            alphaDivType = "Observed",
            plotWidth = 4,
            plotHeight = 3,
            fontSize = 8,
            ptOnlyFlag = TRUE)
