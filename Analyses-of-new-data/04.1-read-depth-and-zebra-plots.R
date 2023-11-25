# 04-read-depth-and-zebra-plots.R
# Author: Greg Poore
# Date: Oct 7, 2023
# Purposes:
# - Analyze changes in nonhuman read depths between host depletion steps
# - Make radial barplots of Zebra coverages

#-------------------------------#
# Load dependencies
require(doMC)
require(plyr)
require(dplyr)
require(tidyr)
require(tibble)
require(reshape2)
require(phyloseq)
require(ggpubr)
require(ggsci)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import metadata and read depth data
# NOTES:
# - This metadata is taken from Narunsky-Haziza et al., which
# has fewer samples (n=15512) than Poore et al. due to sample dropout during 
# additional GRCh38 host depletion and read quality control.
# 
# - The metadata has columns bam_total_reads, bam_mapped_reads, 
# bam_unmapped_reads, and bam_ratio_unmapped that refer to the
# **original** TCGA bam files that were processed using idxstats
# on the Cancer Genomics Cloud (CGC) by SevenBridges (within the "tcga-kraken" project).
# The read counts are thus based on the raw data and not the GRCh38-depleted data.
#
# - The new read depth data only records non-human counts for the following:
# GRCh38, T2T-CHM13, and Pangenome depleted data. Note that these datasets represent
# **sequential** depletion, so T2T-CHM13 denotes host depletion with both GRCh38 and T2T-CHM13;
# similarly, Pangenome data denotes host depletion with GRCh38, T2T-CHM13, and Pangenome
#----------------------------------------------------------#

## Metadata (from mycobiome paper)
load("Supporting_files/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData", verbose=T)
dim(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts) # 15512 43

metaMyco <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts

# Format myco metadata and resave
metaRS <- metaMyco %>% rownames_to_column("mycoIDs")
metaRSFilenames <- gsub("\\.filtered\\.$","",metaRS$run_prefix)
rownames(metaRS) <- metaRSFilenames

## Import read depth data
readDepthsTCGA_WGS <- read.csv("Input_data/read-depths/pese_wgs_hd.tsv",
                               sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(fastq_fn = gsub("\\.filtered$","",fastq_fn))
dim(readDepthsTCGA_WGS) # 4736 5

readDepthsTCGA_RNA <- read.csv("Input_data/read-depths/pese_rna_transcript_hd.tsv",
                               sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(fastq_fn = gsub("\\.$","",fastq_fn)) %>%
  mutate(Pangenome_pese = Transcript_pese) %>%
  select(-Transcript_pese)
dim(readDepthsTCGA_RNA) # 10776 5

## Check overlap of rownames
sum(readDepthsTCGA_WGS$fastq_fn %in% rownames(metaRS)) # 4736
sum(readDepthsTCGA_RNA$fastq_fn %in% rownames(metaRS)) # 10776

## Combine WGS and RNA read depths
readDepthsTCGA_WGS_RNA <- rbind(readDepthsTCGA_WGS,
                                readDepthsTCGA_RNA) %>%
  column_to_rownames("fastq_fn")

metaRS_readDepths <- metaRS %>%
  rownames_to_column("sampleid") %>%
  mutate(hg19 = bam_unmapped_reads) %>%
  mutate(hg38 = readDepthsTCGA_WGS_RNA[sampleid,"HG38"]) %>%
  mutate(T2T = readDepthsTCGA_WGS_RNA[sampleid,"T2T_pese"]) %>%
  mutate(Pangenome = readDepthsTCGA_WGS_RNA[sampleid,"Pangenome_pese"]) %>%
  mutate(percChange_hg19_hg38 = 100*(hg38 - hg19)/hg19) %>%
  mutate(percChange_hg38_T2T = 100*(T2T - hg38)/hg38) %>%
  mutate(percChange_hg38_Pan = 100*(Pangenome - hg38)/hg38) %>%
  mutate(percChange_T2T_Pan = 100*(Pangenome - T2T)/T2T) %>%
  mutate(percTotal_hg19 = 100*hg19/bam_total_reads) %>%
  mutate(percTotal_hg38 = 100*hg38/bam_total_reads) %>%
  mutate(percTotal_T2T = 100*T2T/bam_total_reads) %>%
  mutate(percTotal_Pan = 100*Pangenome/bam_total_reads)

## Non-human read depths
summary(metaRS_readDepths$hg19)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 9.228e+05 5.622e+06 1.003e+07 2.671e+07 2.860e+07 1.935e+09 
summary(metaRS_readDepths$hg38)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 420      3096      4533   1762574    783830 335798226
summary(metaRS_readDepths$T2T)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 154      2762      3687   1031792    149875 240047438
summary(metaRS_readDepths$Pangenome)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 24       182       314    892493     77830 239659174 

## Percent non-human
summary(metaRS_readDepths$percTotal_hg19)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.349   3.510   4.705   7.134   7.675  95.816
summary(metaRS_readDepths$percTotal_hg38)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00025  0.00209  0.00349  0.52030  0.13091 93.42004
summary(metaRS_readDepths$percTotal_T2T)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00009  0.00181  0.00276  0.27644  0.03812 36.18289
summary(metaRS_readDepths$percTotal_Pan)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00002  0.00012  0.00019  0.23793  0.02027 34.15302 

#----------------------------------------------------------#
# Plots read depths by sample type
#----------------------------------------------------------#

#-----------------------Read depth changes-----------------------#
# hg19 vs hg38 vs T2T vs Pan
metaRS_readDepths %>%
  select(experimental_strategy, sample_type, 
         hg19,
         hg38, 
         T2T, 
         Pangenome,
         investigation) %>%
  mutate(investigation = gsub("^TCGA-","",investigation)) %>%
  mutate(experimental_strategy = factor(experimental_strategy, levels = c("WGS","RNA-Seq"))) %>%
  filter(sample_type %in% c("Primary Tumor","Blood Derived Normal", "Solid Tissue Normal")) %>%
  mutate(sample_type = case_when(
    sample_type == "Primary Tumor" ~ "PT",
    sample_type == "Blood Derived Normal" ~ "BDN",
    sample_type == "Solid Tissue Normal" ~ "STN"
  )) %>%
  reshape2::melt(id.vars = c("experimental_strategy","sample_type","investigation")) %>%
  ggboxplot(x = "sample_type",
            y = "value",
            fill = "variable",
            outlier.size = 0.5,
            palette = c("#0072B5FF","#20854EFF","#E18727FF","#BC3C29FF"),
            ylab = "Non-human read count",
            xlab = "",
            legend = "top",
            # add = "mean",
            facet.by = "experimental_strategy") +
  labs(fill = "Sequential host depletion:") +
  scale_y_log10(limits = c(1e1,10^9.5), breaks=c(1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9))
ggsave(filename = "Figures/readDepthChange_hg38_vs_T2T_vs_Pan_17Oct23.jpeg",
       dpi = "retina", units = "in", width = 4, height = 6) 

#-----------------------Percent non-human changes-----------------------#
# hg19 vs hg38 vs T2T vs Pan
metaRS_readDepths %>%
  select(experimental_strategy, sample_type, 
         percTotal_hg19,
         percTotal_hg38, 
         percTotal_T2T, 
         percTotal_Pan,
         investigation) %>%
  rename(hg19=percTotal_hg19, hg38=percTotal_hg38,
         T2T=percTotal_T2T, Pangenome=percTotal_Pan) %>%
  mutate(investigation = gsub("^TCGA-","",investigation)) %>%
  mutate(experimental_strategy = factor(experimental_strategy, levels = c("WGS","RNA-Seq"))) %>%
  filter(sample_type %in% c("Primary Tumor","Blood Derived Normal", "Solid Tissue Normal")) %>%
  mutate(sample_type = case_when(
    sample_type == "Primary Tumor" ~ "PT",
    sample_type == "Blood Derived Normal" ~ "BDN",
    sample_type == "Solid Tissue Normal" ~ "STN"
  )) %>%
  reshape2::melt(id.vars = c("experimental_strategy","sample_type","investigation")) %>%
  ggboxplot(x = "sample_type",
            y = "value",
            fill = "variable",
            outlier.size = 0.5,
            palette = c("#0072B5FF","#20854EFF","#E18727FF","#BC3C29FF"),
            ylab = "Percent of reads that are non-human (%)",
            xlab = "",
            legend = "top",
            # add = "mean",
            facet.by = "experimental_strategy") +
  labs(fill = "Sequential host depletion:") +
  scale_y_log10(limits = c(1e-5, 1e2), breaks = c(1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2))
ggsave(filename = "Figures/percNonhuman_hg19_hg38_vs_T2T_vs_Pan_17Oct23.jpeg",
       dpi = "retina", units = "in", width = 4, height = 6) 

#----------------------------------------------------------#
# Plots read depths by cancer type
#----------------------------------------------------------#

#-----------------------Read depth changes: WGS + RNA-----------------------#
# hg19 vs hg38 vs T2T vs Pan
metaRS_readDepths %>%
  select(experimental_strategy, sample_type, 
         hg19,
         hg38, 
         T2T, 
         Pangenome,
         investigation) %>%
  mutate(investigation = gsub("^TCGA-","",investigation)) %>%
  mutate(experimental_strategy = factor(experimental_strategy, levels = c("WGS","RNA-Seq"))) %>%
  filter(sample_type %in% c("Primary Tumor","Blood Derived Normal", "Solid Tissue Normal")) %>%
  mutate(sample_type = case_when(
    sample_type == "Primary Tumor" ~ "PT",
    sample_type == "Blood Derived Normal" ~ "BDN",
    sample_type == "Solid Tissue Normal" ~ "STN"
  )) %>%
  reshape2::melt(id.vars = c("experimental_strategy","sample_type","investigation")) %>%
  ggboxplot(x = "variable",
            y = "value",
            fill = "variable",
            outlier.size = 0.5,
            # add = "mean",
            palette = c("#0072B5FF","#20854EFF","#E18727FF","#BC3C29FF"),
            ylab = "Non-human read count",
            xlab = "",
            legend = "top") +
  labs(fill = "Sequential host depletion:") +
  scale_y_log10(limits = c(1e1,10^9.5), breaks=c(1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9)) +
  facet_wrap(~investigation, nrow = 6) +
  rotate_x_text(30)
ggsave(filename = "Figures/readDepthChangePerCT_hg38_vs_T2T_vs_Pan_17Oct23.jpeg",
       dpi = "retina", units = "in", width = 12, height = 14) 

#----------------------------------------------------------#
# Zebra plots
#----------------------------------------------------------#

#----------------------Radial barplots by sample type----------------------#
## Load RS210 taxa filt for taxRS210_ff_combTaxaSpecies object
load("Interim_data/taxa_filtering_pipeline_13Oct23.RData", verbose = TRUE)
load("Interim_data/RS210_Filt5050_data_for_ml_tcga_by_seq_center_13Oct23.RData")
dim(taxRS210_ff_combTaxaSpecies) # 3321 12

zebraTCGA_BDN <- read.csv("Input_data/zebra/kuFilepaths_BDN-transcript-removal.zebra_coverages.tsv", 
                          sep = "\t", stringsAsFactors = FALSE)
zebraTCGA_PT <- read.csv("Input_data/zebra/kuFilepaths_PT-transcript-removal.zebra_coverages.tsv", 
                          sep = "\t", stringsAsFactors = FALSE)
zebraTCGA_STN <- read.csv("Input_data/zebra/kuFilepaths_STN-transcript-removal.zebra_coverages.tsv", 
                         sep = "\t", stringsAsFactors = FALSE)

dim(zebraTCGA_BDN)
dim(zebraTCGA_PT)
dim(zebraTCGA_STN)
dim(zebraTCGA_AllFilt)

# Find intersecting genomes
zebraIntOGUs <- Reduce(intersect, list(zebraTCGA_BDN$gotu,
                                       zebraTCGA_PT$gotu,
                                       zebraTCGA_STN$gotu,
                                       rownames(zebraTCGA_AllFilt)))
zebraTCGA_BDNFilt <- zebraTCGA_BDN %>% filter(gotu %in% zebraIntOGUs) %>% column_to_rownames("gotu")
zebraTCGA_PTFilt <- zebraTCGA_PT %>% filter(gotu %in% zebraIntOGUs) %>% column_to_rownames("gotu")
zebraTCGA_STNFilt <- zebraTCGA_STN %>% filter(gotu %in% zebraIntOGUs) %>% column_to_rownames("gotu")
zebraTCGA_AllFilt2 <- zebraTCGA_AllFilt[rownames(zebraTCGA_AllFilt) %in% zebraIntOGUs,]

dim(zebraTCGA_BDNFilt)
dim(zebraTCGA_PTFilt)
dim(zebraTCGA_STNFilt)
dim(zebraTCGA_AllFilt2)

zebraTCGA_CombFilt <- zebraTCGA_AllFilt2 %>%
  rownames_to_column("gotu") %>%
  rename(strain = species) %>%
  mutate(strain = gsub("_"," ",strain)) %>%
  select(-covered_length, -total_length) %>%
  filter(strain != "Cutibacterium acnes") %>%
  rename(All = coverage_ratio) %>%
  mutate(All = round(100*All,1)) %>%
  mutate(BDN = round(100*zebraTCGA_BDNFilt[gotu,"coverage_ratio"],1)) %>%
  mutate(PT = round(100*zebraTCGA_PTFilt[gotu,"coverage_ratio"],1)) %>%
  mutate(STN = round(100*zebraTCGA_STNFilt[gotu,"coverage_ratio"],1)) %>%
  filter(gotu %in% taxRS210_ff_combTaxaSpecies$OGU) %>%
  arrange(desc(BDN), desc(PT))
dim(zebraTCGA_CombFilt) # 3301    6
head(zebraTCGA_CombFilt,20)

## Write function to plot Zebra coverages
plotZebra <- function(selectedSpp = "Streptococcus parasanguinis",
                      zebraInput = zebraTCGA_CombFilt,
                      rX = 2, # constant that controls spacing around center
                      dataString = "rs210Pan",
                      returnPlot = FALSE,
                      sppSize = 3
){
  zebraTCGA_CombFilt_Selected <- zebraInput %>%
    filter(strain == selectedSpp) %>%
    filter(!duplicated(strain)) %>%
    pivot_longer(cols = c("All","BDN","PT","STN"), names_to = "Data_type", values_to = "cov") %>%
    mutate(Text = paste0(cov,"%")) %>%
    bind_rows(data.frame(gotu=rep("",rX),strain=rep("",rX),
                         Data_type=letters[1:rX],cov=rep(0,rX),Text=rep("",rX))) %>%
    mutate(Data_type = factor(Data_type, levels = rev(Data_type)))
  
  ggplot(zebraTCGA_CombFilt_Selected, 
         aes(x = Data_type, y = cov, fill = Data_type, label = Text)) + 
    geom_bar(width = 0.9, stat="identity") + 
    coord_polar(theta = "y") +
    xlab("") + ylab("") +
    ylim(c(0,100)) +
    geom_text(data = zebraTCGA_CombFilt_Selected, hjust = 0, size = 3,
              aes(x = Data_type, y = 1, label = Text),angle=-13,color="black") +
    geom_text(label=gsub(" ","\n",zebraTCGA_CombFilt_Selected$strain[1]), 
              x=.5, y=.5, size=sppSize) +
    theme_minimal() +
    scale_fill_manual(values = c(
                                 "#0072B5FF", # purple > blue
                                 "#7876B1FF", # blue > purple
                                 "#E18727FF", # orange
                                 "#0072B5FF", # purple > blue
                                 "#BC3C29FF", # red
                                 "#7876B1FF"  # blue > purple
                                 )) +
    # Old fill manual (NB: This is very sensitive to changes due to the radial barplot construction)
    # scale_fill_manual(values = c("#7876B1FF",
    #                              "#0072B5FF",
    #                              "#20854EFF","#7876B1FF","#BC3C29FF","#0072B5FF")) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(-1,-1,-1.5,-1.5,"cm")) -> p
  ggsave(plot = p,
         filename = paste0("Figures/zebraPlots/",
                           gsub(" ","_",zebraTCGA_CombFilt_Selected$strain[1]),
                           "_",dataString,
                           ".jpeg"),
         dpi = 900, units = "in", width = 3.5, height = 3.5)
  if(returnPlot){
    return(p)
  }
}

plotZebra()
plotZebra("Streptococcus vestibularis")

zebraTCGA_CombFilt_PTSort <- zebraTCGA_CombFilt %>% arrange(desc(PT))
plotZebra("Fusobacterium nucleatum", zebraTCGA_CombFilt_PTSort)

## Grid top 36 BDN ranked
require(gridExtra)
topXMicrobes <- head(unique(zebraTCGA_CombFilt$strain),50)
pList <- list()
for(ii in 1:36){
  microbeX <- topXMicrobes[ii]
  pList[[ii]] <- plotZebra(microbeX,
                           sppSize = 4,
                           returnPlot = TRUE)
}
n <- length(pList)
nCol <- floor(sqrt(n))
pGrid <- do.call("grid.arrange", c(pList, ncol=nCol,clip=TRUE))
ggsave(plot = pGrid,
       filename = "Figures/zebraPlotsTop36.jpeg",
       dpi = 900, units = "in", width = 20, height = 20)

## Grid top 36 PT ranked
require(gridExtra)
topXMicrobes <- head(unique(zebraTCGA_CombFilt_PTSort$strain),50)
pList <- list()
for(ii in 1:36){
  microbeX <- topXMicrobes[ii]
  pList[[ii]] <- plotZebra(microbeX,
                           sppSize = 4,
                           returnPlot = TRUE)
}
n <- length(pList)
nCol <- floor(sqrt(n))
pGrid <- do.call("grid.arrange", c(pList, ncol=nCol,clip=TRUE))
ggsave(plot = pGrid,
       filename = "Figures/zebraPlotsTop36_PTranked.jpeg",
       dpi = 900, units = "in", width = 20, height = 20)

#----------------------Binned barplot colored by sample type----------------------#

zebraTCGA_CombFilt_Labeled <- zebraTCGA_AllFilt2 %>%
  rownames_to_column("gotu") %>%
  rename(strain = species) %>%
  mutate(strain = gsub("_"," ",strain)) %>%
  select(-covered_length, -total_length) %>%
  filter(strain != "Cutibacterium acnes") %>%
  rename(All = coverage_ratio) %>%
  mutate(BDN = zebraTCGA_BDNFilt[gotu,"coverage_ratio"]) %>%
  mutate(PT = zebraTCGA_PTFilt[gotu,"coverage_ratio"]) %>%
  mutate(STN = zebraTCGA_STNFilt[gotu,"coverage_ratio"]) %>%
  filter(gotu %in% taxRS210_ff_combTaxaSpecies$OGU) %>%
  arrange(desc(All)) %>%
  mutate(Label = case_when(
    (BDN >= PT) & (BDN >= STN) ~ "BDN",
    (STN > PT) & (STN > BDN) ~ "STN",
    (PT > STN) & (PT > BDN) ~ "PT",
  )) %>%
  filter(!duplicated(strain))

zebraTCGA_CombFilt_Labeled %>%
  mutate(All = 100*All) %>%
  mutate(rank = 1:length(gotu)) %>%
  ggplot(aes(x = rank, y=All, fill = Label)) +
  geom_bar(stat='identity') +
  scale_fill_nejm() +
  theme_pubr() +
  labs(x = "Unique species among filtered taxa",
       y = "Total aggregate genome coverage (%)",
       fill = "Most contributing sample type:") +
  theme(legend.position = "top")
ggsave(filename = "Figures/zebraRankedTotalAggregateCov_15Oct23.jpeg",
       dpi = "retina", units = "in", width = 4.5, height = 4)


## Write function to plot Zebra coverages
plotZebraWithoutAll <- function(selectedSpp = "Streptococcus parasanguinis",
                      zebraInput = zebraTCGA_CombFilt,
                      rX = 2, # constant that controls spacing around center
                      dataString = "rs210Pan",
                      returnPlot = FALSE,
                      sppSize = 3
){
  zebraTCGA_CombFilt_Selected <- zebraInput %>%
    filter(strain == selectedSpp) %>%
    filter(!duplicated(strain)) %>%
    select(-All) %>%
    pivot_longer(cols = c("BDN","PT","STN"), names_to = "Data_type", values_to = "cov") %>%
    mutate(Text = paste0(cov,"%")) %>%
    bind_rows(data.frame(gotu=rep("",rX),strain=rep("",rX),
                         Data_type=letters[1:rX],cov=rep(0,rX),Text=rep("",rX))) %>%
    mutate(Data_type = factor(Data_type, levels = rev(Data_type)))
  
  ggplot(zebraTCGA_CombFilt_Selected, 
         aes(x = Data_type, y = cov, fill = Data_type, label = Text)) + 
    geom_bar(width = 0.9, stat="identity") + 
    coord_polar(theta = "y") +
    xlab("") + ylab("") +
    ylim(c(0,100)) +
    geom_text(data = zebraTCGA_CombFilt_Selected, hjust = 0, size = 3,
              aes(x = Data_type, y = 1, label = Text),angle=-13,color="black") +
    geom_text(label=gsub(" ","\n",zebraTCGA_CombFilt_Selected$strain[1]), 
              x=.5, y=.5, size=sppSize) +
    theme_minimal() +
    scale_fill_manual(values = c(
      "#0072B5FF", # purple > blue
      "#7876B1FF", # blue > purple
      "#E18727FF", # orange
      "#0072B5FF", # purple > blue
      "#BC3C29FF", # red
      "#7876B1FF"  # blue > purple
    )) +
    # Old fill manual (NB: This is very sensitive to changes due to the radial barplot construction)
    # scale_fill_manual(values = c("#7876B1FF",
    #                              "#0072B5FF",
    #                              "#20854EFF","#7876B1FF","#BC3C29FF","#0072B5FF")) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(-1,-1,-1.5,-1.5,"cm")) -> p
  ggsave(plot = p,
         filename = paste0("Figures/zebraPlotsWithoutAll/",
                           gsub(" ","_",zebraTCGA_CombFilt_Selected$strain[1]),
                           "_",dataString,
                           ".jpeg"),
         dpi = 900, units = "in", width = 3.5, height = 3.5)
  if(returnPlot){
    return(p)
  }
}
plotZebraWithoutAll()
zebraTCGA_CombFilt_PTSort <- zebraTCGA_CombFilt %>% arrange(desc(PT))
plotZebraWithoutAll("Fusobacterium nucleatum", zebraTCGA_CombFilt_PTSort, sppSize = 4)

#----------------------Experimental strategy----------------------#

zebraTCGA_WGS <- read.csv("Input_data/zebra/specific_pangenome-adapter-filter-pese_13722_zebrasubset.zebra_coverages.tsv", 
                          sep = "\t", stringsAsFactors = FALSE) 
zebraTCGA_RNA <- read.csv("Input_data/zebra/specific_transcript-removal_13767_zebrasubset.zebra_coverages.tsv", 
                         sep = "\t", stringsAsFactors = FALSE) %>% column_to_rownames("gotu")

zebraTCGA_WGS_RNA <- zebraTCGA_WGS %>%
  filter(gotu %in% rownames(zebraTCGA_RNA)) %>%
  rename(cov_WGS = coverage_ratio) %>%
  mutate(cov_RNA = zebraTCGA_RNA[gotu,"coverage_ratio"])

cor.test(log10(zebraTCGA_WGS_RNA$cov_WGS),
         log10(zebraTCGA_WGS_RNA$cov_RNA), method = "spearman")

require(ggpmisc)
zebraTCGA_WGS_RNA %>%
  ggscatter(x = "cov_WGS",
            y = "cov_RNA",
            xlab = "log10(Per microbe WGS coverage ratio)",
            ylab = "log10(Per microbe RNA-Seq coverage)",
            size = 0.1) +
  scale_x_log10() +
  scale_y_log10() +
  stat_cor(method = "spearman",cor.coef.name = "rho") +
  coord_fixed()
ggsave(filename = "Figures/zebra_WGS_vs_RNA_spearman_TCGA_14Oct23.jpeg",
       dpi = "retina", units = "in", width = 4, height = 4)

## Repeat with filtered taxa
zebraTCGA_WGS_RNAFilt <- zebraTCGA_WGS_RNA %>%
  filter(gotu %in% taxRS210_ff_combTaxaSpecies$OGU)

cor.test(log10(zebraTCGA_WGS_RNAFilt$cov_WGS),
         log10(zebraTCGA_WGS_RNAFilt$cov_RNA), method = "pearson")

require(ggpmisc)
zebraTCGA_WGS_RNAFilt %>%
  ggscatter(x = "cov_WGS",
            y = "cov_RNA",
            xlab = "log10(Per microbe WGS coverage ratio)",
            ylab = "log10(Per microbe RNA-Seq coverage)",
            size = 0.1) +
  scale_x_log10() +
  scale_y_log10() +
  stat_cor(method = "spearman",cor.coef.name = "rho") +
  coord_fixed()
ggsave(filename = "Figures/zebraFilt_WGS_vs_RNA_spearman_TCGA_16Oct23.jpeg",
       dpi = "retina", units = "in", width = 4, height = 4)

#----------------------------------------------------------#
# Compare read depths to KrakenUniq counts
#----------------------------------------------------------#

# Format myco metadata and resave
metaKU <- metaMyco %>% rownames_to_column("mycoIDs")
metaKUFilenames <- gsub("\\.filtered\\.$","",metaKU$run_prefix)
rownames(metaKU) <- metaKUFilenames

#----------------Load hg38 data----------------#
## KrakenUniq data
kuHG38 <- read.csv("Input_data/hg38/krakenuniq/combined_kraken_report_tcga_sample_ids_genus.tsv",
                  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(kuHG38)[1] <- "Classification"
dim(kuHG38) # 2133 15513
kuHG38[1:3,1:5]

# Format features
kuHG38Features <- gsub("^k__.+g__","",kuHG38$Classification)
length(kuHG38Features) # 2133
length(unique(kuHG38Features)) # 2133

kuHG38$Classification <- kuHG38Features
kuHG382 <- kuHG38 %>% column_to_rownames("Classification")
kuHG382[1:3,1:3]

kuHG382Tr <- data.frame(t(kuHG382)) %>%
  select(-Homo)
dim(kuHG382Tr) # 15512 2131
kuHG382Tr[1:3,1:3]

# Format KrakenUniq sample IDs
kuHG382TrFor <- kuHG382Tr
rownames(kuHG382TrFor) <- gsub("\\.filtered$","",rownames(kuHG382TrFor))

# Check overlap both ways
sum(metaKUFilenames %in% rownames(kuHG382TrFor)) # 15512
sum(rownames(kuHG382TrFor) %in% metaKUFilenames) # 15512

# Reorder metadata
metaKUOrdHG38 <- metaKU[rownames(kuHG382TrFor),]

# Use myco sample IDs for conciseness
metaKUOrdHG38For <- metaKUOrdHG38 %>%
  rownames_to_column("kuIDs") %>%
  column_to_rownames("mycoIDs")

kuHG382TrFor2 <- kuHG382TrFor
identical(rownames(kuHG382TrFor2), metaKUOrdHG38For$kuIDs) # TRUE
rownames(kuHG382TrFor2) <- rownames(metaKUOrdHG38For)

kuHG382TrFor2[1:3,1:3]
identical(rownames(metaKUOrdHG38For),rownames(kuHG382TrFor2)) # TRUE

# Save as clean tables
metaKUHG38 <- metaKUOrdHG38For
kuHG38Final <- kuHG382TrFor2

#----------------Load T2T data----------------#

## KrakenUniq data
kuT2T <- read.csv("Input_data/t2t-only/krakenuniq/combined_kraken_report_tcga_sample_ids_genus_pese.tsv",
                  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(kuT2T)[1] <- "Classification"
dim(kuT2T) # 2126 15513
kuT2T[1:3,1:5]

# Format features
kuT2TFeatures <- gsub("^k__.+g__","",kuT2T$Classification)
length(kuT2TFeatures) # 2126
length(unique(kuT2TFeatures)) # 2126

kuT2T$Classification <- kuT2TFeatures
kuT2T2 <- kuT2T %>% column_to_rownames("Classification")
kuT2T2[1:3,1:3]

kuT2T2Tr <- data.frame(t(kuT2T2)) %>%
  select(-Homo)
dim(kuT2T2Tr) # 15512 2125
kuT2T2Tr[1:3,1:3]

# Format KrakenUniq sample IDs
kuT2T2TrFor <- kuT2T2Tr
rownames(kuT2T2TrFor) <- gsub("\\.filtered$","",rownames(kuT2T2TrFor))

# Check overlap both ways
sum(metaKUFilenames %in% rownames(kuT2T2TrFor)) # 15512
sum(rownames(kuT2T2TrFor) %in% metaKUFilenames) # 15512

# Reorder metadata
metaKUOrdT2T <- metaKU[rownames(kuT2T2TrFor),]

# Use myco sample IDs for conciseness
metaKUOrdT2TFor <- metaKUOrdT2T %>%
  rownames_to_column("kuIDs") %>%
  column_to_rownames("mycoIDs")

kuT2T2TrFor2 <- kuT2T2TrFor
identical(rownames(kuT2T2TrFor2), metaKUOrdT2TFor$kuIDs) # TRUE
rownames(kuT2T2TrFor2) <- rownames(metaKUOrdT2TFor)

kuT2T2TrFor2[1:3,1:3]
identical(rownames(metaKUOrdT2TFor),rownames(kuT2T2TrFor2)) # TRUE

# Save as clean tables
metaKUT2T <- metaKUOrdT2TFor
kuT2TFinal <- kuT2T2TrFor2

#----------------Load Pangenome non-transcript data----------------#

## KrakenUniq data
kuPan_Nontranscript <- read.csv("Input_data/pangenome/krakenuniq/combined_kraken_report_tcga_sample_ids_genus_pese.tsv",
                  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(kuPan_Nontranscript)[1] <- "Classification"
dim(kuPan_Nontranscript) # 2113 15513
kuPan_Nontranscript[1:3,1:5]

# Format features
kuPanNontranscriptFeatures <- gsub("^k__.+g__","",kuPan_Nontranscript$Classification)
length(kuPanNontranscriptFeatures) # 2113
length(unique(kuPanNontranscriptFeatures)) # 2113

kuPan_Nontranscript$Classification <- kuPanNontranscriptFeatures
kuPan_Nontranscript2 <- kuPan_Nontranscript %>% column_to_rownames("Classification")
kuPan_Nontranscript2[1:3,1:3]

kuPan_Nontranscript2Tr <- data.frame(t(kuPan_Nontranscript2)) %>%
  select(-Homo)
dim(kuPan_Nontranscript2Tr) # 15512 2112
kuPan_Nontranscript2Tr[1:3,1:3]

#----------------Load Pangenome transcript data----------------#

## KrakenUniq data
kuPan_Transcript <- read.csv("Input_data/pangenome/krakenuniq/combined_kraken_report_tcga_sample_ids_genus_transcript_removal.tsv",
                                sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(kuPan_Transcript)[1] <- "Classification"
dim(kuPan_Transcript) # 1936 10777
kuPan_Transcript[1:3,1:5]

# Format features
kuPanTranscriptFeatures <- gsub("^k__.+g__","",kuPan_Transcript$Classification)
length(kuPanTranscriptFeatures) # 1936
length(unique(kuPanTranscriptFeatures)) # 1936

kuPan_Transcript$Classification <- kuPanTranscriptFeatures
kuPan_Transcript2 <- kuPan_Transcript %>% column_to_rownames("Classification")
kuPan_Transcript2[1:3,1:3]

kuPan_Transcript2Tr <- data.frame(t(kuPan_Transcript2)) %>%
  select(-Homo)
dim(kuPan_Transcript2Tr) # 10776  1935
kuPan_Transcript2Tr[1:3,1:3]

#----------------Merge Pangenome non-transcript + transcript data----------------#

sum(rownames(kuPan_Nontranscript2Tr) %in% rownames(kuPan_Transcript2Tr)) # 10776
sum(rownames(kuPan_Transcript2Tr) %in% rownames(kuPan_Nontranscript2Tr)) # 10776

kuPan_Nontranscript2Tr_noRNA <- kuPan_Nontranscript2Tr[!(rownames(kuPan_Nontranscript2Tr) %in% 
                                                           rownames(kuPan_Transcript2Tr)),]
dim(kuPan_Nontranscript2Tr_noRNA) # 4736 2112

# Merge with smartbind
require(gtools)
kuPanPrefinal <- smartbind(kuPan_Nontranscript2Tr_noRNA, 
                          kuPan_Transcript2Tr)
kuPanPrefinal[is.na(kuPanPrefinal)] <- 0
rownames(kuPanPrefinal) <- c(rownames(kuPan_Nontranscript2Tr_noRNA),
                            rownames(kuPan_Transcript2Tr))
dim(kuPanPrefinal) # 15512  2112

# Format KrakenUniq sample IDs
kuPanPrefinalFor <- kuPanPrefinal
rownames(kuPanPrefinalFor) <- gsub("\\.filtered$","",rownames(kuPanPrefinalFor))

# Check overlap both ways
sum(metaKUFilenames %in% rownames(kuPanPrefinalFor)) # 15512
sum(rownames(kuPanPrefinalFor) %in% metaKUFilenames) # 15512

# Reorder metadata
metaKUOrdPan <- metaKU[rownames(kuPanPrefinalFor),]

# Use myco sample IDs for conciseness
metaKUOrdPanFor <- metaKUOrdPan %>%
  rownames_to_column("kuIDs") %>%
  column_to_rownames("mycoIDs")

kuPanPrefinalFor2 <- kuPanPrefinalFor
identical(rownames(kuPanPrefinalFor2), metaKUOrdPanFor$kuIDs) # TRUE
rownames(kuPanPrefinalFor2) <- rownames(metaKUOrdPanFor)

kuPanPrefinalFor2[1:3,1:3]
identical(rownames(metaKUOrdPanFor),rownames(kuPanPrefinalFor2)) # TRUE

# Save as clean tables
metaKUPan <- metaKUOrdPanFor
kuPanFinal <- kuPanPrefinalFor2

#----------------Save the above data----------------#

# save( # Metadata
#   metaKUHG38,
#   metaKUT2T,
#   metaKUPan,
#   # Count data
#   kuHG38Final,
#   kuT2TFinal,
#   kuPanFinal,
#   file = "Interim_data/ku_hg38_T2T_Pan_NoHu_14Oct23.RData")

#----------------Match data types and merge with read depth counts----------------#

load("Interim_data/ku_hg38_T2T_Pan_NoHu_14Oct23.RData")

all(rownames(metaKUHG38) == rownames(metaKUT2T)) # FALSE
sum(rownames(metaKUHG38) %in% rownames(metaKUT2T)) # 15512
sampleInt <- Reduce(intersect, list(rownames(metaKUHG38),
                                    rownames(metaKUT2T),
                                    rownames(metaKUPan)))
length(sampleInt) # 15512
sum(sampleInt %in% metaRS_readDepths$mycoIDs) # 15512
# --> need to reorder the count tables
metaKUPanOrd <- metaKUPan[sampleInt,]
kuHG38FinalOrd <- kuHG38Final[sampleInt,]
kuT2TFinalOrd <- kuT2TFinal[sampleInt,]
kuPanFinalOrd <- kuPanFinal[sampleInt,]

# Reorder and expand the read depth data frame
metaRS_readDepths_Ord <- metaRS_readDepths %>% column_to_rownames("mycoIDs")
metaRS_readDepths_Ord <- metaRS_readDepths_Ord[sampleInt,]
metaRS_readDepths_Ord_Microbe <- metaRS_readDepths_Ord %>%
  mutate(KU_HG38 = rowSums(kuHG38FinalOrd),
         KU_T2T = rowSums(kuT2TFinalOrd),
         KU_Pan = rowSums(kuPanFinalOrd)) %>%
  mutate(ratio_HG38_KU = KU_HG38/hg38,
         ratio_T2T_KU = KU_T2T/T2T,
         ratio_Pan_KU = KU_Pan/Pangenome)

#----------------Plot across levels of host depletion----------------#
metaRS_readDepths_Ord_Microbe %>%
  # filter(KU_Pan != 0) %>%
  mutate(HG38=100*ratio_HG38_KU, 
         T2T=100*ratio_T2T_KU,
         Pangenome=100*ratio_Pan_KU) %>%
  mutate(investigation = gsub("^TCGA-","",investigation)) %>%
  select(sample_type, investigation, experimental_strategy,
         HG38, T2T, Pangenome) %>%
  reshape2::melt(id.vars = c("sample_type", "investigation", "experimental_strategy")) %>%
  mutate(experimental_strategy = factor(experimental_strategy, levels = c("WGS","RNA-Seq"))) %>%
  ggboxplot(x = "variable",
            y = "value",
            fill = "variable",
            # add = "mean",
            palette = c("#0072B5FF","#20854EFF","#E18727FF","#BC3C29FF"),
            ylab = "Percent of non-human reads\nthat are microbial (%)",
            xlab = "",
            legend = "none") +
  labs(fill = "Sequential host depletion:") +
  scale_y_log10() +
  rotate_x_text(0)
ggsave(filename = "Figures/percMicrobial_noFacet_HG38_vs_T2T_vs_Pan_29Oct23.jpeg",
       dpi = "retina", units = "in", width = 2.5, height = 3) 

#----------------Plot correlation with input counts----------------#
## HG38 -- All
corCountsHG38 <- metaRS_readDepths_Ord_Microbe %>%
  # filter(KU_HG38!=0) %>%
  mutate(sample_type = gsub("Primary Blood Derived Cancer - Peripheral Blood",
                            "Primary Blood Cancer",sample_type)) %>%
  ggscatter(x = "hg38",
            y = "KU_HG38",
            xlab = "Non-human reads\nafter hg38 depletion",
            ylab = "Microbial reads\nafter hg38 depletion",
            size = 0.2) +
  # coord_fixed() +
  scale_x_log10() +
  scale_y_log10() +
  # geom_smooth(method='lm') +
  stat_cor(method = "spearman",cor.coef.name = "rho")
ggsave(plot = corCountsHG38,
       filename = "Figures/corInputVsMicrobialReads_HG38_All_14Oct23.jpeg",
       dpi = "retina", units = "in", width = 3, height = 3)

## T2T -- All
corCountsT2T <- metaRS_readDepths_Ord_Microbe %>%
  # filter(KU_T2T!=0) %>%
  mutate(sample_type = gsub("Primary Blood Derived Cancer - Peripheral Blood",
                            "Primary Blood Cancer",sample_type)) %>%
  ggscatter(x = "T2T",
            y = "KU_T2T",
            xlab = "Non-human reads\nafter T2T depletion",
            ylab = "Microbial reads\nafter T2T depletion",
            size = 0.2) +
  # coord_fixed() +
  scale_x_log10() +
  scale_y_log10() +
  # geom_smooth(method='lm') +
  stat_cor(method = "spearman",cor.coef.name = "rho")
ggsave(plot = corCountsT2T,
       filename = "Figures/corInputVsMicrobialReads_T2T_All_14Oct23.jpeg",
       dpi = "retina", units = "in", width = 3, height = 3)

## Pan -- All
corCountsPan <- metaRS_readDepths_Ord_Microbe %>%
  # filter(KU_Pan!=0) %>%
  mutate(sample_type = gsub("Primary Blood Derived Cancer - Peripheral Blood",
                            "Primary Blood Cancer",sample_type)) %>%
  ggscatter(x = "Pangenome",
            y = "KU_Pan",
            xlab = "Non-human reads\nafter Pangenome depletion",
            ylab = "Microbial reads\nafter Pangenome depletion",
            size = 0.2) +
  # coord_fixed() +
  scale_x_log10() +
  scale_y_log10() +
  # geom_smooth(method='lm') +
  stat_cor(method = "spearman",cor.coef.name = "rho")
ggsave(plot = corCountsPan,
       filename = "Figures/corInputVsMicrobialReads_Pangenome_All_14Oct23.jpeg",
       dpi = "retina", units = "in", width = 3, height = 3)

## Aggregate the "All" plots using gridExtra
require(gridExtra)
pCorList <- list(corCountsHG38, corCountsT2T, corCountsPan)
pCorGrid <- do.call("grid.arrange", c(pCorList, ncol=3,clip=TRUE))
ggsave(plot = pCorGrid,
       filename = "Figures/corInputVsMicrobialReads_Combined_HG38_T2T_Pangenome.jpeg",
       dpi = 900, units = "in", width = 10, height = 4)

#----------------------------------------------------------#
# Compare KrakenUniq counts across HG38 vs T2T vs Pangenome
#----------------------------------------------------------#
# Use data frame with microbe counts from above

#--------------Mean fold change per CT: HG38 vs T2T--------------#
kuFoldChange_HG38_T2T <- metaRS_readDepths_Ord_Microbe %>%
  filter(KU_T2T != 0) %>%
  droplevels() %>%
  # filter(sample_type == "Primary Tumor") %>%
  mutate(investigationShort = gsub("^TCGA-","",investigation)) %>%
  mutate(foldChange_KU_HG38_T2T = KU_HG38 / KU_T2T) %>%
  group_by(investigationShort) %>%
  summarise(meanFC = mean(foldChange_KU_HG38_T2T),
            sdFC = sd(foldChange_KU_HG38_T2T),
            seFC = sd(foldChange_KU_HG38_T2T)/sqrt(length(foldChange_KU_HG38_T2T)),
            medianFC = median(foldChange_KU_HG38_T2T),
            madFC = mad(foldChange_KU_HG38_T2T),
            numSamp = length(foldChange_KU_HG38_T2T)) %>% 
  filter(numSamp!=1) %>%
  arrange(meanFC) %>% data.frame()

library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# Plot fold changes
kuFoldChange_HG38_T2T %>%
  mutate(investigationShortOrd = reorder(investigationShort, -meanFC)) %>%
  ggplot(aes(x = investigationShortOrd,y=meanFC, fill=investigationShortOrd)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_text(aes(y = meanFC+seFC,label = round(meanFC,1)), 
            hjust = 1.2, angle = 90, size = 3.5) +
  geom_text(aes(y = 1,label = paste0("",numSamp)), 
            hjust = -0.2, vjust=0.5, angle = 90, size = 2, color="blue") +
  geom_errorbar(aes(ymin=meanFC-seFC, ymax=meanFC+seFC), width=.2,
                position=position_dodge(.9)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(size=8)) +
  # New code to reverse axes (next 2 lines)
  scale_y_continuous(trans = reverselog_trans(10), 
                     limits = c(500,1)) +
  scale_x_discrete(position = "top") +
  # scale_y_log10(limits = c(1,1e3)) +
  annotation_logticks(sides="l") +
  # geom_hline(yintercept = mean(kuFoldChange_HG38_T2T$meanFC), # =24.61567
  #            linetype=2, color = "darkred") +
  rotate_x_text(90) +
  theme(legend.position = "none", text = element_text(size = 16)) +
  labs(x = "",y="Mean fold change\nKrakenUniq reads (HG38/T2T)") +
  scale_fill_manual(values=rev(as.vector(pals::coolwarm(33))))
ggsave(filename = "Figures/kuMeanFoldChange_HG38_vs_T2T_16Oct23.jpeg",
       dpi = "retina", units = "in", width = 6, height = 6)

#--------------Mean fold change per CT: HG38 vs Pangenome--------------#
kuFoldChange_HG38_Pan <- metaRS_readDepths_Ord_Microbe %>%
  filter(KU_Pan != 0) %>% 
  droplevels() %>%
  # filter(sample_type == "Primary Tumor") %>%
  mutate(investigationShort = gsub("^TCGA-","",investigation)) %>%
  mutate(foldChange_KU_HG38_Pan = KU_HG38 / KU_Pan) %>%
  group_by(investigationShort) %>%
  summarise(meanFC = mean(foldChange_KU_HG38_Pan),
            sdFC = sd(foldChange_KU_HG38_Pan),
            seFC = sd(foldChange_KU_HG38_Pan)/sqrt(length(foldChange_KU_HG38_Pan)),
            medianFC = median(foldChange_KU_HG38_Pan),
            madFC = mad(foldChange_KU_HG38_Pan),
            numSamp = length(foldChange_KU_HG38_Pan)) %>% 
  filter(numSamp!=1) %>%
  arrange(meanFC) %>% data.frame()

# Plot fold changes
kuFoldChange_HG38_Pan %>%
  mutate(investigationShortOrd = reorder(investigationShort, -meanFC)) %>%
  ggplot(aes(x = investigationShortOrd,y=meanFC, fill=investigationShortOrd)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_text(aes(y = meanFC+seFC,label = round(meanFC,1)), 
            hjust = 1.2, angle = 90, size = 3.5) +
  geom_text(aes(y = 1,label = paste0("",numSamp)), 
            hjust = -0.2, vjust=0.5, angle = 90, size = 2, color="blue") +
  geom_errorbar(aes(ymin=meanFC-seFC, ymax=meanFC+seFC), width=.2,
                position=position_dodge(.9)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(size=8)) +
  scale_y_continuous(trans = reverselog_trans(10), 
                     limits = c(500,1)) +
  scale_x_discrete(position = "top") +
  # scale_y_log10(limits = c(1,1e3)) +
  annotation_logticks(sides="l") +
  # geom_hline(yintercept = mean(kuFoldChange_HG38_Pan$meanFC), # =32.25861
  #            linetype=2, color = "darkred") +
  rotate_x_text(90) +
  theme(legend.position = "none", text = element_text(size = 16)) +
  labs(x = "",y="Mean fold change\nKrakenUniq reads (HG38/Pangenome)") +
  scale_fill_manual(values=rev(as.vector(pals::coolwarm(28))))
ggsave(filename = "Figures/kuMeanFoldChange_HG38_vs_Pan_16Oct23.jpeg",
       dpi = "retina", units = "in", width = 6, height = 6)

#--------------Mean fold change per ES: HG38 vs T2T--------------#
kuFoldChangeES_HG38_T2T <- metaRS_readDepths_Ord_Microbe %>%
  filter(KU_T2T != 0) %>%
  droplevels() %>%
  # filter(sample_type == "Primary Tumor") %>%
  mutate(investigationShort = gsub("^TCGA-","",investigation)) %>%
  mutate(foldChange_KU_HG38_T2T = KU_HG38 / KU_T2T) %>%
  group_by(experimental_strategy) %>%
  summarise(meanFC = mean(foldChange_KU_HG38_T2T),
            sdFC = sd(foldChange_KU_HG38_T2T),
            seFC = sd(foldChange_KU_HG38_T2T)/sqrt(length(foldChange_KU_HG38_T2T)),
            medianFC = median(foldChange_KU_HG38_T2T),
            madFC = mad(foldChange_KU_HG38_T2T),
            numSamp = length(foldChange_KU_HG38_T2T)) %>% 
  filter(numSamp!=1) %>%
  arrange(meanFC) %>% data.frame()

kuFoldChangeES_HG38_T2T %>%
  # mutate(investigationShortOrd = reorder(investigationShort, -meanFC)) %>%
  ggplot(aes(x = experimental_strategy,y=meanFC, fill=experimental_strategy)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_text(aes(y = meanFC+seFC,label = round(meanFC,2)), 
            hjust = -0.2, angle = 0, size = 3.5) +
  geom_text(aes(y = 1,label = paste0("",numSamp)), 
            hjust = 0.5, vjust=-3, angle = 90, size = 2, color="blue") +
  geom_errorbar(aes(ymin=meanFC-seFC, ymax=meanFC+seFC), width=.2,
                position=position_dodge(.9)) +
  theme_pubr() + 
  coord_flip() +
  theme(legend.position = "none", text = element_text(size = 16)) +
  labs(x = "",y="Mean fold change\nKrakenUniq reads (hg38/T2T)") +
  scale_fill_nejm()
ggsave(filename = "Figures/kuMeanFoldChangeES_HG38_vs_T2T_17Oct23.jpeg",
       dpi = "retina", units = "in", width = 6, height = 2)

#--------------Mean fold change per ES: HG38 vs Pan--------------#
kuFoldChangeES_HG38_Pan <- metaRS_readDepths_Ord_Microbe %>%
  filter(KU_Pan != 0) %>%
  droplevels() %>%
  # filter(sample_type == "Primary Tumor") %>%
  mutate(investigationShort = gsub("^TCGA-","",investigation)) %>%
  mutate(foldChange_KU_HG38_Pan = KU_HG38 / KU_Pan) %>%
  group_by(experimental_strategy) %>%
  summarise(meanFC = mean(foldChange_KU_HG38_Pan),
            sdFC = sd(foldChange_KU_HG38_Pan),
            seFC = sd(foldChange_KU_HG38_Pan)/sqrt(length(foldChange_KU_HG38_Pan)),
            medianFC = median(foldChange_KU_HG38_Pan),
            madFC = mad(foldChange_KU_HG38_Pan),
            numSamp = length(foldChange_KU_HG38_Pan)) %>% 
  filter(numSamp!=1) %>%
  arrange(meanFC) %>% data.frame()

kuFoldChangeES_HG38_Pan %>%
  # mutate(investigationShortOrd = reorder(investigationShort, -meanFC)) %>%
  ggplot(aes(x = experimental_strategy,y=meanFC, fill=experimental_strategy)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_text(aes(y = meanFC+seFC,label = round(meanFC,2)), 
            hjust = -0.2, angle = 0, size = 3.5) +
  geom_text(aes(y = 1,label = paste0("",numSamp)), 
            hjust = 0.5, vjust=-3, angle = 90, size = 2, color="blue") +
  geom_errorbar(aes(ymin=meanFC-seFC, ymax=meanFC+seFC), width=.2,
                position=position_dodge(.9)) +
  theme_pubr() + 
  ylim(c(0,55)) +
  coord_flip() +
  theme(legend.position = "none", text = element_text(size = 16)) +
  labs(x = "",y="Mean fold change\nKrakenUniq reads (hg38/Pangenome)") +
  scale_fill_nejm()
ggsave(filename = "Figures/kuMeanFoldChangeES_HG38_vs_Pan_17Oct23.jpeg",
       dpi = "retina", units = "in", width = 6, height = 2)

