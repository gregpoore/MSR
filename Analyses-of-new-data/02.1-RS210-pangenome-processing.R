# 02-RS210-pangenome-processing.R
# Author: Greg Poore
# Date: Oct 7, 2023
# Purposes:
# - Load RS210 pangenome-filtered data and format

#-------------------------------#
# Load dependencies
require(doMC)
require(plyr)
require(dplyr)
require(tibble)
require(reshape2)
require(phyloseq)
require(biomformat)
require(rhdf5)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import TCGA metadata
#----------------------------------------------------------#
## Metadata (from mycobiome paper)
load("Supporting_files/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData", verbose=T)
dim(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts) # 15512 43

metaMyco <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts

#----------------------------------------------------------#
# Import TCGA RS210 pangenome data (without transcript removal)
#----------------------------------------------------------#
## RS210 pangenome data
rs210PanBIOM <- read_biom(biom_file = "Input_data/pangenome/RS210/tcga-rs210-pangenome-merged-feature-table-pese.biom")
rs210Pan <- t(as(biom_data(rs210PanBIOM), "matrix"))
dim(rs210Pan) # 15510 23377
rs210Pan[1:3,1:3]

rs210PanSampleIDs <- rownames(rs210Pan)
rs210PanSampleIDs2 <- gsub("\\.filtered\\.R1.+","",rs210PanSampleIDs)
rs210PanSampleIDs3 <- gsub("\\.R1\\.trimmed.+","",rs210PanSampleIDs2)
rownames(rs210Pan) <- rs210PanSampleIDs3
rs210Pan[1:3,1:3]

rs210PanOGUs <- colnames(rs210Pan)
# save(rs210PanOGUs,
#      file = "Interim_data/rs210PanOGUs_13Oct23.RData")

#----------------------------------------------------------#
# Align metadata and count data
# Adjustments made:
# - WGS filenames in RS210 data had ".filtered" suffix but RNA filenames did not
# - Removed ".filtered" suffix from all RS210 filenames
# - After ordering, used the myco metadata sample IDs for conciseness
#----------------------------------------------------------#
# Format myco metadata and resave
metaRS <- metaMyco %>% rownames_to_column("mycoIDs")
metaRSFilenames <- gsub("\\.filtered\\.$","",metaRS$run_prefix)
rownames(metaRS) <- metaRSFilenames

# # Format RS210 sample IDs
# Check overlap both ways
sum(metaRSFilenames %in% rownames(rs210Pan)) # 15510
sum(rownames(rs210Pan) %in% metaRSFilenames) # 15510

# Reorder metadata
metaRSOrd <- metaRS[rownames(rs210Pan),]

# Use myco sample IDs for conciseness
metaRSOrdFor <- metaRSOrd %>%
  rownames_to_column("rsIDs") %>%
  column_to_rownames("mycoIDs")

rs210Pan2 <- rs210Pan
identical(rownames(rs210Pan2), metaRSOrdFor$rsIDs) # TRUE
rownames(rs210Pan2) <- rownames(metaRSOrdFor)

# Check final results
head(metaRSOrdFor)
rs210Pan2[1:3,1:3]
identical(rownames(metaRSOrdFor),rownames(rs210Pan2)) # TRUE

# Save as clean tables
metaRSPrefinal <- metaRSOrdFor
rs210PanPrefinal <- rs210Pan2

# Check for any zero-sum samples and remove -- none are removed
metaRSPrefinalZeroSum <- metaRSPrefinal[which(rowSums(rs210PanPrefinal)==0),]
table(metaRSPrefinalZeroSum$experimental_strategy) # 0

metaRSPrefinalNonzero <- droplevels(metaRSPrefinal[which(rowSums(rs210PanPrefinal)!=0),])
rs210PanPrefinalNonzero <- rs210PanPrefinal[which(rowSums(rs210PanPrefinal)!=0),]
dim(metaRSPrefinalNonzero) # 15510    44

# save(metaRSFinal,
#      metaRSFinalNonzero,
#      rs210PanFinal,
#      rs210PanFinalNonzero,
#      file = "Interim_data/metadata_count_data_full_RS210_13Sept23.RData")

#----------------------------------------------------------#
# Replace TCGA RS210 pangenome RNA-Seq data with post-transcript removal data
#----------------------------------------------------------#
##---------Subset above data to WGS only---------##
metaRSPrefinalNonzero_WGS <- metaRSPrefinalNonzero %>%
  filter(experimental_strategy == "WGS") %>% droplevels()
rs210PanPrefinalNonzero_WGS <- rs210PanPrefinalNonzero[rownames(metaRSPrefinalNonzero_WGS),]
rs210PanPrefinalNonzero_WGSFilt <- rs210PanPrefinalNonzero_WGS[,colSums(rs210PanPrefinalNonzero_WGS)!=0]

##---------RS210 transcript removal---------##
rs210TranscriptBIOM <- read_biom(biom_file = "Input_data/pangenome/RS210/tcga-rs210-pangenome-merged-feature-table-transcript-removal.biom")
rs210Transcript <- t(as(biom_data(rs210TranscriptBIOM), "matrix"))
dim(rs210Transcript) # 9346 22888
rs210Transcript[1:3,1:3]

rs210TranscriptSampleIDs <- rownames(rs210Transcript)
rs210TranscriptSampleIDs2 <- gsub("\\.filtered\\.R1.+","",rs210TranscriptSampleIDs)
rs210TranscriptSampleIDs3 <- gsub("\\.R1\\.trimmed.+","",rs210TranscriptSampleIDs2)
rownames(rs210Transcript) <- rs210TranscriptSampleIDs3
rs210Transcript[1:3,1:3]

##---------Reformat metadata for RNA-Seq data---------##

# Check overlap both ways
sum(metaRSFilenames %in% rownames(rs210Transcript)) # 9346
sum(rownames(rs210Transcript) %in% metaRSFilenames) # 9346

# Reorder metadata
metaRSOrdTranscript <- metaRS[rownames(rs210Transcript),]

# Use myco sample IDs for conciseness
metaRSOrdTranscriptFor <- metaRSOrdTranscript %>%
  rownames_to_column("rsIDs") %>%
  column_to_rownames("mycoIDs")

rs210Transcript2 <- rs210Transcript
identical(rownames(rs210Transcript2), metaRSOrdTranscriptFor$rsIDs) # TRUE
rownames(rs210Transcript2) <- rownames(metaRSOrdTranscriptFor)

# Check final results
head(metaRSOrdTranscriptFor)
rs210Transcript2[1:3,1:3]
identical(rownames(metaRSOrdTranscriptFor),rownames(rs210Transcript2)) # TRUE

# Save as clean tables
metaRSTranscriptPrefinal <- metaRSOrdTranscriptFor
rs210TranscriptPrefinal <- rs210Transcript2

# Check for any zero-sum samples and remove -- none are removed
metaRSTranscriptPrefinalZeroSum <- metaRSTranscriptPrefinal[which(rowSums(rs210TranscriptPrefinal)==0),]
table(metaRSTranscriptPrefinalZeroSum$experimental_strategy) # 0

metaRSTranscriptPrefinalNonzero <- droplevels(metaRSTranscriptPrefinal[which(rowSums(rs210TranscriptPrefinal)!=0),])
rs210TranscriptPrefinalNonzero <- rs210TranscriptPrefinal[which(rowSums(rs210TranscriptPrefinal)!=0),]
dim(metaRSTranscriptPrefinalNonzero) # 9346   44

##---------Merge WGS and RNA-Seq data---------##
# Merge using smartbind
require(gtools)
rs210PanFinal <- smartbind(rs210PanPrefinalNonzero_WGSFilt, 
                          rs210TranscriptPrefinalNonzero)
rs210PanFinal[is.na(rs210PanFinal)] <- 0
rownames(rs210PanFinal) <- c(rownames(rs210PanPrefinalNonzero_WGSFilt),
                            rownames(rs210TranscriptPrefinalNonzero))
dim(rs210PanFinal) # 14080 23381

metaRSFinal <- rbind(metaRSPrefinalNonzero_WGS,
                            metaRSTranscriptPrefinalNonzero)
dim(metaRSFinal) # 14080    44

# Verify that no samples have zero counts
all(rowSums(rs210PanFinal)!=0) # TRUE
metaRSFinalNonzero <- metaRSFinal

# Verify identical row/sample order
all(rownames(rs210PanFinal) == rownames(metaRSFinal)) # TRUE

# Subset features
rs210PanFinalOGUs <- colnames(rs210PanFinal)

# save(rs210PanFinalOGUs,
#      file = "Interim_data/rs210PanFinalOGUs_WGS_RNA_Transcript_13Oct23.RData")
# save(rs210PanFinal,
#      rs210PanFinalOGUs,
#      metaRSFinal,
#      file = "Interim_data/rs210PanFinal_WGS_RNA_Transcript_13Oct23.RData")

#----------------------------------------------------------#
# Import data subsets:
# - WIS data summarized at species level (NOTE: This was already intersected using the RS210 lineage map)
# - Taxa filtering pipeline results
#----------------------------------------------------------#
wisRS210 <- read.csv("Supporting_data/taxa_RS210_WIS_All_Biological_Samples_Hits_16Aug23.csv",
                     stringsAsFactors = FALSE)

rs210Tax <- read.csv("Supporting_data/RS210_lineages.csv",
                     stringsAsFactors = FALSE)

load("Interim_data/taxa_filtering_pipeline_13Oct23.RData", verbose = TRUE)
# taxRS210_ff_combTaxaSpecies
# taxRS210_ff_combTaxaSpeciesZebra
# taxRS210_ff_combTaxaSpeciesZebraWithViruses
# taxRS210_ff_combTaxaSpeciesZebraWithViruses_9090
# taxRS210_ff_combTaxaSpeciesZebraWithViruses_7575
# taxRS210_ff_combTaxaSpeciesZebraWithViruses_5080
# taxRS210_ff_combTaxaSpeciesZebraWithViruses_5075
# taxRS210_ff_combTaxaSpeciesZebraWithViruses_5050
# taxRS210_ff_virusCov
# taxRS210_ff
# zebraTCGA_AllFilt

#----------------------------------------------------------#
# Intersect RS210-derived TCGA data with feature subsets
#----------------------------------------------------------#

#-----------------WIS-----------------#

# Create WIS intersected dataframe
rs210PanFinalWIS <- rs210PanFinal[, wisRS210$OGU]
dim(rs210PanFinalWIS) # 14080   813

# Check for zero-sum samples
sum(rowSums(rs210PanFinal)==0) # 0
sum(rowSums(rs210PanFinalWIS)==0) # 6899

# Check how many are due to lower RNA-Seq depth
metaRSFinalWISZeroSum <- metaRSFinal[which(rowSums(rs210PanFinalWIS)==0),]
table(metaRSFinalWISZeroSum$experimental_strategy) # WGS=68, RNA-Seq=6831

metaRSFinalWISNonzero <- droplevels(metaRSFinal[which(rowSums(rs210PanFinalWIS)!=0),])
rs210PanFinalWISNonzero <- rs210PanFinalWIS[which(rowSums(rs210PanFinalWIS)!=0),]

sum(rowSums(rs210PanFinal)==0) # 0
sum(rowSums(rs210PanFinalWISNonzero)==0) # 0

#-----------------Taxa filtering: 9090-----------------#
rs210PanFinal9090 <- rs210PanFinal[, taxRS210_ff_combTaxaSpeciesZebraWithViruses_9090$OGU]
metaRSFinal[which(rowSums(rs210PanFinal9090)==0),] %>%
  count(experimental_strategy) # WGS=220, RNA-Seq=6882

metaRSFinal9090_Nonzero <- droplevels(metaRSFinal[which(rowSums(rs210PanFinal9090)!=0),])
rs210PanFinal9090_Nonzero <- rs210PanFinal9090[which(rowSums(rs210PanFinal9090)!=0),]

#-----------------Taxa filtering: 7575-----------------#
rs210PanFinal7575 <- rs210PanFinal[, taxRS210_ff_combTaxaSpeciesZebraWithViruses_7575$OGU]
metaRSFinal[which(rowSums(rs210PanFinal7575)==0),] %>%
  count(experimental_strategy) # WGS=114, RNA-Seq=6881

metaRSFinal7575_Nonzero <- droplevels(metaRSFinal[which(rowSums(rs210PanFinal7575)!=0),])
rs210PanFinal7575_Nonzero <- rs210PanFinal7575[which(rowSums(rs210PanFinal7575)!=0),]

#-----------------Taxa filtering: 5050-----------------#
rs210PanFinal5050 <- rs210PanFinal[, taxRS210_ff_combTaxaSpeciesZebraWithViruses_5050$OGU]
metaRSFinal[which(rowSums(rs210PanFinal5050)==0),] %>%
  count(experimental_strategy) # WGS=80, RNA-Seq=6881

metaRSFinal5050_Nonzero <- droplevels(metaRSFinal[which(rowSums(rs210PanFinal5050)!=0),])
rs210PanFinal5050_Nonzero <- rs210PanFinal5050[which(rowSums(rs210PanFinal5050)!=0),]

#-----------------Taxa filtering: 5080-----------------#
rs210PanFinal5080 <- rs210PanFinal[, taxRS210_ff_combTaxaSpeciesZebraWithViruses_5080$OGU]
metaRSFinal[which(rowSums(rs210PanFinal5080)==0),] %>%
  count(experimental_strategy) # WGS=82, RNA-Seq=6881

metaRSFinal5080_Nonzero <- droplevels(metaRSFinal[which(rowSums(rs210PanFinal5080)!=0),])
rs210PanFinal5080_Nonzero <- rs210PanFinal5080[which(rowSums(rs210PanFinal5080)!=0),]

## Save for supp tables
rs210PanFinal5050_Nonzero %>%
  rownames_to_column("sampleid") %>%
  write.csv("Supp_tables/rs210PanFinal5050_Nonzero.csv",
            row.names = FALSE)
metaRSFinal5050_Nonzero %>%
  rownames_to_column("sampleid") %>%
  write.csv("Supp_tables/metaRSFinal5050_Nonzero.csv",
            row.names = FALSE)

#----------------------------------------------------------#
# Examining read distributions
#----------------------------------------------------------#

##------------------Create plots of dropped samples------------------##
# Percent of dropped samples per experimental strategy
all(rownames(rs210PanFinal5050) == rownames(metaRSFinal)) # TRUE

metaRSFinal %>%
  mutate(nonzero5050 = ifelse(rowSums(rs210PanFinal5050)!=0, TRUE, FALSE)) %>%
  count(nonzero5050, experimental_strategy) %>%
  mutate(experimental_strategy = factor(experimental_strategy, levels = c("WGS","RNA-Seq"))) %>%
  mutate(nonzero5050 = factor(nonzero5050, levels = c(TRUE,FALSE))) %>%
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
            fill = "nonzero5050") +
  rotate_x_text(30)
ggsave("Figures/rs210Pan_5050_barplot_nonzero_samples_WGS_RNA_18Oct.jpeg",
       dpi = "retina", units = "in", width = 1.5, height = 5)

# Relate back to original file sizes
metaRSFinal %>%
  mutate(nonzeroSum = (rowSums(rs210PanFinal5050)!=0)) %>%
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
ggsave(filename = "Figures/rs210Pan_5050_bamTotalReads_barplot_nonzero_vs_zero_samples_WGS_RNA_18Oct.jpeg",
       units = "in", width = 2, height = 5, dpi = "retina")

##------------------Read count plots------------------##

source("00-Functions.R") # for plotReadCounts()
plotReadCounts(metaData = metaRSFinal5050_Nonzero,
               countData = rs210PanFinal5050_Nonzero,
               dataString = "rs210Pan_Filt5050")

#----------------------------------------------------------#
# Subset to Illumina HiSeq
#----------------------------------------------------------#

metaRSFinalNonzero %>% count(data_submitting_center_label, experimental_strategy)
metaRSFinalNonzero %>% count(cgc_platform) # HiSeq accounts for 97.0% of samples (13657/14080)
metaRSFinalNonzero %>% count(cgc_platform, data_submitting_center_label)
metaRSFinalWISNonzero %>% count(cgc_platform, data_submitting_center_label)
metaRSFinalWISNonzero %>% count(portion_is_ffpe, data_submitting_center_label)

# Subset metadata to Illumina HiSeq samples
# There are only 2 Hopkins HiSeq samples --> remove
# There are only 10 MDA RPPA center HiSeq samples --> remove
# Only 2 samples are YES for portion_is_ffpe --> remove ones that are not NO
metaRSFinalNonzero_HiSeq <- metaRSFinalNonzero %>%
  filter(!is.na(data_submitting_center_label)) %>%
  filter(portion_is_ffpe == "NO") %>%
  filter(!(data_submitting_center_label %in% c("Johns Hopkins / University of Southern California",
           "MD Anderson - RPPA Core Facility (Proteomics)"))) %>%
  filter(cgc_platform == "Illumina HiSeq") %>% droplevels()
# WIS
metaRSFinalWISNonzero_HiSeq <- metaRSFinalNonzero_HiSeq %>%
  rownames_to_column("sampleID") %>%
  filter(sampleID %in% rownames(metaRSFinalWISNonzero)) %>%
  column_to_rownames("sampleID") %>% droplevels()
# Filt 9090
metaRSFinal9090_Nonzero_HiSeq <- metaRSFinalNonzero_HiSeq %>%
  rownames_to_column("sampleID") %>%
  filter(sampleID %in% rownames(metaRSFinal9090_Nonzero)) %>%
  column_to_rownames("sampleID") %>% droplevels()
# Filt 7575
metaRSFinal7575_Nonzero_HiSeq <- metaRSFinalNonzero_HiSeq %>%
  rownames_to_column("sampleID") %>%
  filter(sampleID %in% rownames(metaRSFinal7575_Nonzero)) %>%
  column_to_rownames("sampleID") %>% droplevels()
# Filt 5050
metaRSFinal5050_Nonzero_HiSeq <- metaRSFinalNonzero_HiSeq %>%
  rownames_to_column("sampleID") %>%
  filter(sampleID %in% rownames(metaRSFinal5050_Nonzero)) %>%
  column_to_rownames("sampleID") %>% droplevels()

# Subset count data to Illumina HiSeq samples
rs210PanFinalNonzero_HiSeq <- rs210PanFinal[rownames(metaRSFinalNonzero_HiSeq),]
rs210PanFinalWISNonzero_HiSeq <- rs210PanFinalWISNonzero[rownames(metaRSFinalWISNonzero_HiSeq),]
rs210PanFinal9090_Nonzero_HiSeq <- rs210PanFinal9090_Nonzero[rownames(metaRSFinal9090_Nonzero_HiSeq),]
rs210PanFinal7575_Nonzero_HiSeq <- rs210PanFinal7575_Nonzero[rownames(metaRSFinal7575_Nonzero_HiSeq),]
rs210PanFinal5050_Nonzero_HiSeq <- rs210PanFinal5050_Nonzero[rownames(metaRSFinal5050_Nonzero_HiSeq),]

dim(rs210PanFinalNonzero_HiSeq) # 13614 23381
dim(rs210PanFinalWISNonzero_HiSeq) # 6719  813
dim(rs210PanFinal9090_Nonzero_HiSeq) # 6517  508
dim(rs210PanFinal7575_Nonzero_HiSeq) # 6622  765
dim(rs210PanFinal5050_Nonzero_HiSeq) # 6655 1189

# save(rs210PanFinalNonzero_HiSeq,
#      rs210PanFinalWISNonzero_HiSeq,
#      rs210PanFinal9090_Nonzero_HiSeq,
#      rs210PanFinal7575_Nonzero_HiSeq,
#      rs210PanFinal5050_Nonzero_HiSeq,
#      metaRSFinalWISNonzero_HiSeq,
#      file = "Interim_data/rs210Pan_HiSeq_data_metadata_21Oct23.RData")

#----------------------------------------------------------#
# Subset full data to individual seq centers: Full data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaRSFinalNonzero_HiSeq_HMS <- metaRSFinalNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaRSFinalNonzero_HiSeq_BCM <- metaRSFinalNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaRSFinalNonzero_HiSeq_MDA <- metaRSFinalNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaRSFinalNonzero_HiSeq_WashU <- metaRSFinalNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaRSFinalNonzero_HiSeq_Broad_WGS <- metaRSFinalNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, but RNA-Seq is only GBM, only a WGS is made)
metaRSFinalNonzero_HiSeq_UNC <- metaRSFinalNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaRSFinalNonzero_HiSeq_CMS <- metaRSFinalNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rs210PanFinalNonzero_HiSeq_HMS <- rs210PanFinalNonzero_HiSeq[rownames(metaRSFinalNonzero_HiSeq_HMS),]
rs210PanFinalNonzero_HiSeq_BCM <- rs210PanFinalNonzero_HiSeq[rownames(metaRSFinalNonzero_HiSeq_BCM),]
rs210PanFinalNonzero_HiSeq_MDA <- rs210PanFinalNonzero_HiSeq[rownames(metaRSFinalNonzero_HiSeq_MDA),]
rs210PanFinalNonzero_HiSeq_WashU <- rs210PanFinalNonzero_HiSeq[rownames(metaRSFinalNonzero_HiSeq_WashU),]
rs210PanFinalNonzero_HiSeq_Broad_WGS <- rs210PanFinalNonzero_HiSeq[rownames(metaRSFinalNonzero_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
rs210PanFinalNonzero_HiSeq_UNC <- rs210PanFinalNonzero_HiSeq[rownames(metaRSFinalNonzero_HiSeq_UNC),]
rs210PanFinalNonzero_HiSeq_CMS <- rs210PanFinalNonzero_HiSeq[rownames(metaRSFinalNonzero_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS
metaRSFinalNonzero_HiSeq_WGS <- metaRSFinalNonzero_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
rs210PanFinalNonzero_HiSeq_WGS <- rs210PanFinalNonzero_HiSeq[rownames(metaRSFinalNonzero_HiSeq_WGS),]

# RNA
metaRSFinalNonzero_HiSeq_RNA <- metaRSFinalNonzero_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
rs210PanFinalNonzero_HiSeq_RNA <- rs210PanFinalNonzero_HiSeq[rownames(metaRSFinalNonzero_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  rs210PanFinalNonzero_HiSeq_HMS,
  rs210PanFinalNonzero_HiSeq_BCM,
  rs210PanFinalNonzero_HiSeq_MDA,
  rs210PanFinalNonzero_HiSeq_WashU,
  rs210PanFinalNonzero_HiSeq_Broad_WGS,
  rs210PanFinalNonzero_HiSeq_UNC,
  rs210PanFinalNonzero_HiSeq_CMS,
  
  # WGS and RNA subset data
  rs210PanFinalNonzero_HiSeq_WGS,
  rs210PanFinalNonzero_HiSeq_RNA,
  
  # Subset metadata
  metaRSFinalNonzero_HiSeq_HMS,
  metaRSFinalNonzero_HiSeq_BCM,
  metaRSFinalNonzero_HiSeq_MDA,
  metaRSFinalNonzero_HiSeq_WashU,
  metaRSFinalNonzero_HiSeq_Broad_WGS,
  metaRSFinalNonzero_HiSeq_UNC,
  metaRSFinalNonzero_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaRSFinalNonzero_HiSeq_WGS,
  metaRSFinalNonzero_HiSeq_RNA,
  
  # Full metadata
  metaRSFinalNonzero_HiSeq,
  file = "Interim_data/RS210_Full_data_for_ml_tcga_by_seq_center_13Oct23.RData")

#----------------------------------------------------------#
# Subset full data to individual seq centers: WIS data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaRSFinalWISNonzero_HiSeq_HMS <- metaRSFinalWISNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaRSFinalWISNonzero_HiSeq_BCM <- metaRSFinalWISNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaRSFinalWISNonzero_HiSeq_MDA <- metaRSFinalWISNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaRSFinalWISNonzero_HiSeq_WashU <- metaRSFinalWISNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaRSFinalWISNonzero_HiSeq_Broad_WGS <- metaRSFinalWISNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, but RNA-Seq is only GBM, only a WGS is made)
metaRSFinalWISNonzero_HiSeq_UNC <- metaRSFinalWISNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaRSFinalWISNonzero_HiSeq_CMS <- metaRSFinalWISNonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rs210PanFinalWISNonzero_HiSeq_HMS <- rs210PanFinalWISNonzero_HiSeq[rownames(metaRSFinalWISNonzero_HiSeq_HMS),]
rs210PanFinalWISNonzero_HiSeq_BCM <- rs210PanFinalWISNonzero_HiSeq[rownames(metaRSFinalWISNonzero_HiSeq_BCM),]
rs210PanFinalWISNonzero_HiSeq_MDA <- rs210PanFinalWISNonzero_HiSeq[rownames(metaRSFinalWISNonzero_HiSeq_MDA),]
rs210PanFinalWISNonzero_HiSeq_WashU <- rs210PanFinalWISNonzero_HiSeq[rownames(metaRSFinalWISNonzero_HiSeq_WashU),]
rs210PanFinalWISNonzero_HiSeq_Broad_WGS <- rs210PanFinalWISNonzero_HiSeq[rownames(metaRSFinalWISNonzero_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
rs210PanFinalWISNonzero_HiSeq_UNC <- rs210PanFinalWISNonzero_HiSeq[rownames(metaRSFinalWISNonzero_HiSeq_UNC),]
rs210PanFinalWISNonzero_HiSeq_CMS <- rs210PanFinalWISNonzero_HiSeq[rownames(metaRSFinalWISNonzero_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS
metaRSFinalWISNonzero_HiSeq_WGS <- metaRSFinalWISNonzero_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
rs210PanFinalWISNonzero_HiSeq_WGS <- rs210PanFinalWISNonzero_HiSeq[rownames(metaRSFinalWISNonzero_HiSeq_WGS),]

# RNA
metaRSFinalWISNonzero_HiSeq_RNA <- metaRSFinalWISNonzero_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
rs210PanFinalWISNonzero_HiSeq_RNA <- rs210PanFinalWISNonzero_HiSeq[rownames(metaRSFinalWISNonzero_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  rs210PanFinalWISNonzero_HiSeq_HMS,
  rs210PanFinalWISNonzero_HiSeq_BCM,
  rs210PanFinalWISNonzero_HiSeq_MDA,
  rs210PanFinalWISNonzero_HiSeq_WashU,
  rs210PanFinalWISNonzero_HiSeq_Broad_WGS,
  rs210PanFinalWISNonzero_HiSeq_UNC,
  rs210PanFinalWISNonzero_HiSeq_CMS,
  
  # WGS and RNA subset data
  rs210PanFinalWISNonzero_HiSeq_WGS,
  rs210PanFinalWISNonzero_HiSeq_RNA,
  
  # Subset metadata
  metaRSFinalWISNonzero_HiSeq_HMS,
  metaRSFinalWISNonzero_HiSeq_BCM,
  metaRSFinalWISNonzero_HiSeq_MDA,
  metaRSFinalWISNonzero_HiSeq_WashU,
  metaRSFinalWISNonzero_HiSeq_Broad_WGS,
  metaRSFinalWISNonzero_HiSeq_UNC,
  metaRSFinalWISNonzero_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaRSFinalWISNonzero_HiSeq_WGS,
  metaRSFinalWISNonzero_HiSeq_RNA,
  
  # Full metadata
  metaRSFinalWISNonzero_HiSeq,
  file = "Interim_data/RS210_WIS_data_for_ml_tcga_by_seq_center_13Oct23.RData")

#----------------------------------------------------------#
# Subset full data to individual seq centers: Filt 9090 data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaRSFinal9090_Nonzero_HiSeq_HMS <- metaRSFinal9090_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaRSFinal9090_Nonzero_HiSeq_BCM <- metaRSFinal9090_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaRSFinal9090_Nonzero_HiSeq_MDA <- metaRSFinal9090_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaRSFinal9090_Nonzero_HiSeq_WashU <- metaRSFinal9090_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaRSFinal9090_Nonzero_HiSeq_Broad_WGS <- metaRSFinal9090_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, but RNA-Seq is only GBM, only a WGS is made)
metaRSFinal9090_Nonzero_HiSeq_UNC <- metaRSFinal9090_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaRSFinal9090_Nonzero_HiSeq_CMS <- metaRSFinal9090_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rs210PanFinal9090_Nonzero_HiSeq_HMS <- rs210PanFinal9090_Nonzero_HiSeq[rownames(metaRSFinal9090_Nonzero_HiSeq_HMS),]
rs210PanFinal9090_Nonzero_HiSeq_BCM <- rs210PanFinal9090_Nonzero_HiSeq[rownames(metaRSFinal9090_Nonzero_HiSeq_BCM),]
rs210PanFinal9090_Nonzero_HiSeq_MDA <- rs210PanFinal9090_Nonzero_HiSeq[rownames(metaRSFinal9090_Nonzero_HiSeq_MDA),]
rs210PanFinal9090_Nonzero_HiSeq_WashU <- rs210PanFinal9090_Nonzero_HiSeq[rownames(metaRSFinal9090_Nonzero_HiSeq_WashU),]
rs210PanFinal9090_Nonzero_HiSeq_Broad_WGS <- rs210PanFinal9090_Nonzero_HiSeq[rownames(metaRSFinal9090_Nonzero_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
rs210PanFinal9090_Nonzero_HiSeq_UNC <- rs210PanFinal9090_Nonzero_HiSeq[rownames(metaRSFinal9090_Nonzero_HiSeq_UNC),]
rs210PanFinal9090_Nonzero_HiSeq_CMS <- rs210PanFinal9090_Nonzero_HiSeq[rownames(metaRSFinal9090_Nonzero_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS
metaRSFinal9090_Nonzero_HiSeq_WGS <- metaRSFinal9090_Nonzero_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
rs210PanFinal9090_Nonzero_HiSeq_WGS <- rs210PanFinal9090_Nonzero_HiSeq[rownames(metaRSFinal9090_Nonzero_HiSeq_WGS),]

# RNA
metaRSFinal9090_Nonzero_HiSeq_RNA <- metaRSFinal9090_Nonzero_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
rs210PanFinal9090_Nonzero_HiSeq_RNA <- rs210PanFinal9090_Nonzero_HiSeq[rownames(metaRSFinal9090_Nonzero_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  rs210PanFinal9090_Nonzero_HiSeq_HMS,
  rs210PanFinal9090_Nonzero_HiSeq_BCM,
  rs210PanFinal9090_Nonzero_HiSeq_MDA,
  rs210PanFinal9090_Nonzero_HiSeq_WashU,
  rs210PanFinal9090_Nonzero_HiSeq_Broad_WGS,
  rs210PanFinal9090_Nonzero_HiSeq_UNC,
  rs210PanFinal9090_Nonzero_HiSeq_CMS,
  
  # WGS and RNA subset data
  rs210PanFinal9090_Nonzero_HiSeq_WGS,
  rs210PanFinal9090_Nonzero_HiSeq_RNA,
  
  # Subset metadata
  metaRSFinal9090_Nonzero_HiSeq_HMS,
  metaRSFinal9090_Nonzero_HiSeq_BCM,
  metaRSFinal9090_Nonzero_HiSeq_MDA,
  metaRSFinal9090_Nonzero_HiSeq_WashU,
  metaRSFinal9090_Nonzero_HiSeq_Broad_WGS,
  metaRSFinal9090_Nonzero_HiSeq_UNC,
  metaRSFinal9090_Nonzero_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaRSFinal9090_Nonzero_HiSeq_WGS,
  metaRSFinal9090_Nonzero_HiSeq_RNA,
  
  # Full metadata
  metaRSFinal9090_Nonzero_HiSeq,
  file = "Interim_data/RS210_Filt9090_data_for_ml_tcga_by_seq_center_13Oct23.RData")

#----------------------------------------------------------#
# Subset full data to individual seq centers: Filt 7575 data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaRSFinal7575_Nonzero_HiSeq_HMS <- metaRSFinal7575_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaRSFinal7575_Nonzero_HiSeq_BCM <- metaRSFinal7575_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaRSFinal7575_Nonzero_HiSeq_MDA <- metaRSFinal7575_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaRSFinal7575_Nonzero_HiSeq_WashU <- metaRSFinal7575_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaRSFinal7575_Nonzero_HiSeq_Broad_WGS <- metaRSFinal7575_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, but RNA-Seq is only GBM, only a WGS is made)
metaRSFinal7575_Nonzero_HiSeq_UNC <- metaRSFinal7575_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaRSFinal7575_Nonzero_HiSeq_CMS <- metaRSFinal7575_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rs210PanFinal7575_Nonzero_HiSeq_HMS <- rs210PanFinal7575_Nonzero_HiSeq[rownames(metaRSFinal7575_Nonzero_HiSeq_HMS),]
rs210PanFinal7575_Nonzero_HiSeq_BCM <- rs210PanFinal7575_Nonzero_HiSeq[rownames(metaRSFinal7575_Nonzero_HiSeq_BCM),]
rs210PanFinal7575_Nonzero_HiSeq_MDA <- rs210PanFinal7575_Nonzero_HiSeq[rownames(metaRSFinal7575_Nonzero_HiSeq_MDA),]
rs210PanFinal7575_Nonzero_HiSeq_WashU <- rs210PanFinal7575_Nonzero_HiSeq[rownames(metaRSFinal7575_Nonzero_HiSeq_WashU),]
rs210PanFinal7575_Nonzero_HiSeq_Broad_WGS <- rs210PanFinal7575_Nonzero_HiSeq[rownames(metaRSFinal7575_Nonzero_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
rs210PanFinal7575_Nonzero_HiSeq_UNC <- rs210PanFinal7575_Nonzero_HiSeq[rownames(metaRSFinal7575_Nonzero_HiSeq_UNC),]
rs210PanFinal7575_Nonzero_HiSeq_CMS <- rs210PanFinal7575_Nonzero_HiSeq[rownames(metaRSFinal7575_Nonzero_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS
metaRSFinal7575_Nonzero_HiSeq_WGS <- metaRSFinal7575_Nonzero_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
rs210PanFinal7575_Nonzero_HiSeq_WGS <- rs210PanFinal7575_Nonzero_HiSeq[rownames(metaRSFinal7575_Nonzero_HiSeq_WGS),]

# RNA
metaRSFinal7575_Nonzero_HiSeq_RNA <- metaRSFinal7575_Nonzero_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
rs210PanFinal7575_Nonzero_HiSeq_RNA <- rs210PanFinal7575_Nonzero_HiSeq[rownames(metaRSFinal7575_Nonzero_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  rs210PanFinal7575_Nonzero_HiSeq_HMS,
  rs210PanFinal7575_Nonzero_HiSeq_BCM,
  rs210PanFinal7575_Nonzero_HiSeq_MDA,
  rs210PanFinal7575_Nonzero_HiSeq_WashU,
  rs210PanFinal7575_Nonzero_HiSeq_Broad_WGS,
  rs210PanFinal7575_Nonzero_HiSeq_UNC,
  rs210PanFinal7575_Nonzero_HiSeq_CMS,
  
  # WGS and RNA subset data
  rs210PanFinal7575_Nonzero_HiSeq_WGS,
  rs210PanFinal7575_Nonzero_HiSeq_RNA,
  
  # Subset metadata
  metaRSFinal7575_Nonzero_HiSeq_HMS,
  metaRSFinal7575_Nonzero_HiSeq_BCM,
  metaRSFinal7575_Nonzero_HiSeq_MDA,
  metaRSFinal7575_Nonzero_HiSeq_WashU,
  metaRSFinal7575_Nonzero_HiSeq_Broad_WGS,
  metaRSFinal7575_Nonzero_HiSeq_UNC,
  metaRSFinal7575_Nonzero_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaRSFinal7575_Nonzero_HiSeq_WGS,
  metaRSFinal7575_Nonzero_HiSeq_RNA,
  
  # Full metadata
  metaRSFinal7575_Nonzero_HiSeq,
  file = "Interim_data/RS210_Filt7575_data_for_ml_tcga_by_seq_center_13Oct23.RData")

#----------------------------------------------------------#
# Subset full data to individual seq centers: Filt 5050 data
#----------------------------------------------------------#

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaRSFinal5050_Nonzero_HiSeq_HMS <- metaRSFinal5050_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaRSFinal5050_Nonzero_HiSeq_BCM <- metaRSFinal5050_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaRSFinal5050_Nonzero_HiSeq_MDA <- metaRSFinal5050_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaRSFinal5050_Nonzero_HiSeq_WashU <- metaRSFinal5050_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaRSFinal5050_Nonzero_HiSeq_Broad_WGS <- metaRSFinal5050_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, but RNA-Seq is only GBM, only a WGS is made)
metaRSFinal5050_Nonzero_HiSeq_UNC <- metaRSFinal5050_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaRSFinal5050_Nonzero_HiSeq_CMS <- metaRSFinal5050_Nonzero_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rs210PanFinal5050_Nonzero_HiSeq_HMS <- rs210PanFinal5050_Nonzero_HiSeq[rownames(metaRSFinal5050_Nonzero_HiSeq_HMS),]
rs210PanFinal5050_Nonzero_HiSeq_BCM <- rs210PanFinal5050_Nonzero_HiSeq[rownames(metaRSFinal5050_Nonzero_HiSeq_BCM),]
rs210PanFinal5050_Nonzero_HiSeq_MDA <- rs210PanFinal5050_Nonzero_HiSeq[rownames(metaRSFinal5050_Nonzero_HiSeq_MDA),]
rs210PanFinal5050_Nonzero_HiSeq_WashU <- rs210PanFinal5050_Nonzero_HiSeq[rownames(metaRSFinal5050_Nonzero_HiSeq_WashU),]
rs210PanFinal5050_Nonzero_HiSeq_Broad_WGS <- rs210PanFinal5050_Nonzero_HiSeq[rownames(metaRSFinal5050_Nonzero_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
rs210PanFinal5050_Nonzero_HiSeq_UNC <- rs210PanFinal5050_Nonzero_HiSeq[rownames(metaRSFinal5050_Nonzero_HiSeq_UNC),]
rs210PanFinal5050_Nonzero_HiSeq_CMS <- rs210PanFinal5050_Nonzero_HiSeq[rownames(metaRSFinal5050_Nonzero_HiSeq_CMS),]

#--------------------Subset metadata and count data by WGS and RNA (for multi-class classification)--------------------#
# WGS
metaRSFinal5050_Nonzero_HiSeq_WGS <- metaRSFinal5050_Nonzero_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% 
  droplevels()
rs210PanFinal5050_Nonzero_HiSeq_WGS <- rs210PanFinal5050_Nonzero_HiSeq[rownames(metaRSFinal5050_Nonzero_HiSeq_WGS),]

# RNA
metaRSFinal5050_Nonzero_HiSeq_RNA <- metaRSFinal5050_Nonzero_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% 
  droplevels()
rs210PanFinal5050_Nonzero_HiSeq_RNA <- rs210PanFinal5050_Nonzero_HiSeq[rownames(metaRSFinal5050_Nonzero_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  rs210PanFinal5050_Nonzero_HiSeq_HMS,
  rs210PanFinal5050_Nonzero_HiSeq_BCM,
  rs210PanFinal5050_Nonzero_HiSeq_MDA,
  rs210PanFinal5050_Nonzero_HiSeq_WashU,
  rs210PanFinal5050_Nonzero_HiSeq_Broad_WGS,
  rs210PanFinal5050_Nonzero_HiSeq_UNC,
  rs210PanFinal5050_Nonzero_HiSeq_CMS,
  
  # WGS and RNA subset data
  rs210PanFinal5050_Nonzero_HiSeq_WGS,
  rs210PanFinal5050_Nonzero_HiSeq_RNA,
  
  # Subset metadata
  metaRSFinal5050_Nonzero_HiSeq_HMS,
  metaRSFinal5050_Nonzero_HiSeq_BCM,
  metaRSFinal5050_Nonzero_HiSeq_MDA,
  metaRSFinal5050_Nonzero_HiSeq_WashU,
  metaRSFinal5050_Nonzero_HiSeq_Broad_WGS,
  metaRSFinal5050_Nonzero_HiSeq_UNC,
  metaRSFinal5050_Nonzero_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaRSFinal5050_Nonzero_HiSeq_WGS,
  metaRSFinal5050_Nonzero_HiSeq_RNA,
  
  # Full metadata
  metaRSFinal5050_Nonzero_HiSeq,
  file = "Interim_data/RS210_Filt5050_data_for_ml_tcga_by_seq_center_13Oct23.RData")

#----------------------------------------------------------#
# Batch correction with ConQuR
# NOTE: See "Supporting_scripts/S06-RS210-ConQuR-correction/conqur-rs210-ALL-wgs-rna.R"
#----------------------------------------------------------#
# Load data from above
load("Interim_data/RS210_WIS_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/RS210_Filt9090_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/RS210_Filt7575_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/RS210_Filt5050_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/taxa_filtering_pipeline_13Oct23.RData")

# Load ConQuR-corrected data
load("Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ_CMS.RData")

#----------------------------------------------------------#
# Calculate alpha and beta diversity per batch using phyloseq
#----------------------------------------------------------#
require(ggsci)
require(ggpubr)

source("00-Functions.R") # for alphaBetaFXN() and runAlphaBetaSeqCenter()
# Done
runAlphaBetaSeqCenter(metaString = "metaRSFinal5050_Nonzero_HiSeq",
                      countString = "rs210PanFinal5050_Nonzero_HiSeq",
                      dataStringInput = "rs210Pan_Filt5050")

## On CQ data
alphaBetaFXN(metaData = metaRSFinal9090_Nonzero_HiSeq_WGS,
             countData = rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ,
             dataString = "rs210Pan_Filt9090_WGS_CQ")

#----------------------------------------------------------#
# Recalculate **species** level alpha and beta diversity per batch using phyloseq
#----------------------------------------------------------#
load("Interim_data/rs210Pan_HiSeq_data_metadata_21Oct23.RData")

source("00-Functions.R")
runAlphaBetaSeqCenter_TaxaLevel(metaData = metaRSFinal5050_Nonzero_HiSeq,
                                countData = rs210PanFinal5050_Nonzero_HiSeq,
                                dataStringInput = "rs210Pan_Filt5050",
                                taxaLevel = "Species")

#----------------------------------------------------#
# Export per-batch data to run Qiime2
#----------------------------------------------------#

## Create and write taxa file
taxaFile <- data.frame(`Feature ID` = taxRS210_ff_combTaxaSpeciesZebraWithViruses$OGU,
                       Taxon = taxRS210_ff_combTaxaSpeciesZebraWithViruses$OGU,
                       check.names = FALSE)
write.table(taxaFile,
            file = "Qiime_data_and_scripts/Qiime_input_data/rs210-taxa.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

source("00-Functions.R") # for export2Qiime()
export2Qiime(metaString = "metaRSFinal5050_Nonzero_HiSeq",
             countString = "rs210PanFinal5050_Nonzero_HiSeq",
             dataString = "rs210Pan_Filt5050")

#----------------------------------------------------#
# Export ConQuR-corrected data to run Qiime2
#----------------------------------------------------#

# Load ConQuR-corrected data
load("Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinalWISNonzero_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal9090_Nonzero_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal7575_Nonzero_HiSeq_RNA_CQ_CMS.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_WGS_CQ_BCM.RData")
load("Input_data/conqur-data/rs210/rs210PanFinal5050_Nonzero_HiSeq_RNA_CQ_CMS.RData")

source("00-Functions.R") # for export2QiimeCQ()
export2QiimeCQ(metaString = "metaRSFinal5050_Nonzero_HiSeq",
             countString = "rs210PanFinal5050_Nonzero_HiSeq",
             dataString = "rs210Pan_Filt5050")

#----------------------------------------------------#
# Calculate per-center differential abundances
#----------------------------------------------------#

# Load data from above
load("Interim_data/RS210_WIS_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/RS210_Filt9090_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/RS210_Filt7575_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/RS210_Filt5050_data_for_ml_tcga_by_seq_center_13Oct23.RData")
load("Interim_data/taxa_filtering_pipeline_13Oct23.RData")

source("00-Functions.R") # for runAncomBC_1VsAll_OGUs()
# Filt5050
runAncomBC_1VsAll_OGUs(metaString = "metaRSFinal5050_Nonzero_HiSeq",
                       countString = "rs210PanFinal5050_Nonzero_HiSeq",
                       dataString = "rs210Pan_Filt5050",
                       taxTable = taxRS210_ff_combTaxaSpeciesZebraWithViruses,
                       qvalCutoff = 0.05,
                       showTopXFlag = FALSE,
                       sampleTypes = c("Primary Tumor","Blood Derived Normal"),
                       SeqCenters = c("HMS","BCM","MDA","WashU","Broad_WGS","UNC","CMS"),
                       taxaPlotLabel = "genus")

metaRSFinal5050_Nonzero_HiSeq_WGS_CQ <- metaRSFinal5050_Nonzero_HiSeq_WGS
runAncomBC_1VsAll_OGUs(metaString = "metaRSFinal5050_Nonzero_HiSeq_WGS",
                       countString = "rs210PanFinal5050_Nonzero_HiSeq_WGS",
                       dataString = "rs210Pan_Filt5050_WGS_CQ",
                       taxTable = taxRS210_ff_combTaxaSpeciesZebraWithViruses,
                       qvalCutoff = 0.05,
                       showTopXFlag = FALSE,
                       sampleTypes = c("Primary Tumor","Blood Derived Normal"),
                       SeqCenters = c("CQ"),
                       taxaPlotLabel = "genus")

metaRSFinal5050_Nonzero_HiSeq_RNA_CQ <- metaRSFinal5050_Nonzero_HiSeq_RNA
runAncomBC_1VsAll_OGUs(metaString = "metaRSFinal5050_Nonzero_HiSeq_RNA",
                       countString = "rs210PanFinal5050_Nonzero_HiSeq_RNA",
                       dataString = "rs210Pan_Filt5050_RNA_CQ",
                       taxTable = taxRS210_ff_combTaxaSpeciesZebraWithViruses,
                       qvalCutoff = 0.05,
                       showTopXFlag = FALSE,
                       sampleTypes = c("Primary Tumor"),
                       SeqCenters = c("CQ"),
                       taxaPlotLabel = "genus")

#----------------------------------------------------#
# Replot ANCOM-BC results
#----------------------------------------------------#

source("00-Functions.R") # for ancomReplot()
ancomReplot(seqCenter = "HMS",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050",
            fontSize = 5,
            pointSize = 0.5)
ancomReplot(seqCenter = "HMS",
            sampleType = "Blood Derived Normal",
            dataString = "rs210Pan_Filt5050",
            fontSize = 5,
            pointSize = 0.5)

# BCM
ancomReplot(seqCenter = "BCM",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050")
ancomReplot(seqCenter = "BCM",
            sampleType = "Blood Derived Normal",
            dataString = "rs210Pan_Filt5050")
# WashU
ancomReplot(seqCenter = "WashU",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
ancomReplot(seqCenter = "WashU",
            sampleType = "Blood Derived Normal",
            dataString = "rs210Pan_Filt5050",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
# Broad_WGS
ancomReplot(seqCenter = "Broad_WGS",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050",
            fontSize = 7,
            pointSize = 0.5)
ancomReplot(seqCenter = "Broad_WGS",
            sampleType = "Blood Derived Normal",
            dataString = "rs210Pan_Filt5050",
            fontSize = 7,
            pointSize = 0.5)
# MDA
ancomReplot(seqCenter = "MDA",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
ancomReplot(seqCenter = "MDA",
            sampleType = "Blood Derived Normal",
            dataString = "rs210Pan_Filt5050",
            numticks = 2,
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
# CMS
ancomReplot(seqCenter = "CMS",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)

#----------------------------------------------------#
# Replot alpha diversity results
#----------------------------------------------------#

source("00-Functions.R") # for alphaReplot()
# HMS
alphaReplot(seqCenter = "HMS",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050_Species",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# BCM
alphaReplot(seqCenter = "BCM",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050_Species",
            alphaDivType = "Observed",
            plotWidth = 6,
            plotHeight = 3,
            fontSize = 8)
# Broad_WGS
alphaReplot(seqCenter = "Broad_WGS",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050_Species",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# WashU
alphaReplot(seqCenter = "WashU",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050_Species",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# MDA
alphaReplot(seqCenter = "MDA",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050_Species",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# CMS
alphaReplot(seqCenter = "CMS",
            sampleType = "Primary Tumor",
            dataString = "rs210Pan_Filt5050_Species",
            alphaDivType = "Observed",
            plotWidth = 4,
            plotHeight = 3,
            fontSize = 8,
            ptOnlyFlag = TRUE)




