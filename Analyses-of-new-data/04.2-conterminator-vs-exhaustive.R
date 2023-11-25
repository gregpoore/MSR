# 04.2-conterminator-vs-exhaustive.R
# Author: Greg Poore
# Date: Oct 17, 2023
# Purposes:
# - Import fasta files from Conterminator and Exhaustive outputs
# - Compare number of genomes and regions affected by each

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
require(seqinr)
require(ggvenn)
require(rstatix)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import WoLr1 fasta files and format into data frames
#----------------------------------------------------------#
##---------------------Conterminator---------------------##
dbWoLr1_Conterminator <- read.fasta("Supporting_files/Database_scrubbing/WoLr1/hg38.t2t.pangenome.wol1-exk32.conterminator.fna")
length(dbWoLr1_Conterminator) # 369 --> NOTE: this is not a unique genome count

dbWoLr1_ConterminatorDf <- data.frame(gID = gsub("\\:.+","",getName(dbWoLr1_Conterminator)),
                                      startStop = gsub("^G.+\\:","",getName(dbWoLr1_Conterminator))) %>%
  separate("startStop", c("start","stop"),
           sep = "-", remove = TRUE, convert = TRUE) %>%
  mutate(length = stop-start+1)
length(unique(dbWoLr1_ConterminatorDf$gID)) # 114 --> this is the number of unique genomes affected

##---------------------Exhaustive---------------------##
dbWoLr1_Exhaustive <- read.fasta("Supporting_files/Database_scrubbing/WoLr1/hg38.t2t.pangenome.wol1-exk32.exhaustive.fna")
length(dbWoLr1_Exhaustive) # 676

dbWoLr1_ExhaustiveDf <- data.frame(gID = gsub("\\:.+","",getName(dbWoLr1_Exhaustive)),
                                      startStop = gsub("^G.+\\:","",getName(dbWoLr1_Exhaustive))) %>%
  separate("startStop", c("start","stop"),
           sep = "-", remove = TRUE, convert = TRUE) %>%
  mutate(length = stop-start+1)
length(unique(dbWoLr1_ExhaustiveDf$gID)) # 229 --> this is the number of unique genomes affected

##---------------------Compare genomes---------------------##
sum(dbWoLr1_ConterminatorDf$gID %in% dbWoLr1_ExhaustiveDf$gID) # 362
sum(dbWoLr1_ExhaustiveDf$gID %in% dbWoLr1_ConterminatorDf$gID) # 492

dbWoLr1_totalUnique <- unique(c(dbWoLr1_ConterminatorDf$gID,
                                dbWoLr1_ExhaustiveDf$gID))
length(dbWoLr1_totalUnique) # 236 --> 2.2% of 10575

dbWoLr1_List <- list(
  Conterminator = unique(dbWoLr1_ConterminatorDf$gID),
  Exhaustive = unique(dbWoLr1_ExhaustiveDf$gID)
)

require(ggvenn)
ggvenn(dbWoLr1_List, text_size = 3.5, show_percentage = FALSE) +
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF"))
ggsave("Figures/dbCompare_Venn_gIDs_WoLr1_Conterminator_Exhaustive.jpeg",
       dpi = "retina",units = "in", width = 3, height = 2)

##---------------------Compare distributions of lengths---------------------##
require(ggridges)
dbWoLr1_BothDf <- rbind(cbind(dbWoLr1_ConterminatorDf,method="Conterminator"),
                        cbind(dbWoLr1_ExhaustiveDf,method="Exhaustive")) %>%
  mutate(length = as.numeric(length))

## Unpaired comparisons
# Summarise lengths across genomes
dbWoLr1_BothDf_summarised <- dbWoLr1_BothDf %>%
  group_by(gID, method) %>%
  summarise(lengthComb = sum(length)) %>%
  arrange(gID)
# Calc unpaired Wilcoxon
dbWoLr1_BothDf_summarised %>% as.data.frame() %>% wilcox_test(lengthComb ~ method)
# Plot using ggridges
dbWoLr1_BothDf_summarised %>%
  ggplot(aes(x=log10(lengthComb),y=method,fill=stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01,
                               jittered_points = TRUE,
                               position = position_points_jitter(width = 0.05, height = 0),
                               point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  scale_fill_viridis_c(name = "Genome length", option = "C", alpha=0.7) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(x = "log10(length of contamination\nin affected genomes)",
       y = "")
ggsave(filename = "Figures/dbCompare_lengths_WoLr1_Conterminator_Exhaustive.jpeg",
       dpi = "retina", units = "in", width = 4, height = 3)

## Paired comparisons
dbWoLr1_Overlapping_gIDs <- dbWoLr1_BothDf_summarised$gID[duplicated(dbWoLr1_BothDf_summarised$gID)]
dbWoLr1_BothDf_summarised_paired <- dbWoLr1_BothDf_summarised %>%
  filter(gID %in% dbWoLr1_Overlapping_gIDs)
# Calc paired Wilcoxon
dbWoLr1_BothDf_summarised_paired %>%
  as.data.frame() %>%
  mutate(method = factor(method)) %>%
  wilcox_test(lengthComb ~ method, paired = TRUE) -> dbPairedWilcoxonWoLr1
# Plot paired differences
dbWoLr1_BothDf_summarised_paired %>%
  mutate(loglengthComb = log10(lengthComb)) %>%
  ggplot(aes(x=method,y=loglengthComb)) +
  geom_boxplot(aes(fill=method),notch = TRUE) +
  geom_point(size = 0.1, alpha=0.2) +
  geom_line(aes(group = gID),alpha=0.2,color="gray") +
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF")) +
  theme_pubr() + theme(legend.position = "none") +
  rotate_x_text(30) +
  labs(x = "", y = "log10(contamination length)") +
  # stat_compare_means(paired = TRUE)
  stat_pvalue_manual(dbPairedWilcoxonWoLr1, label = "p = {p}", y.position = 4.5)
ggsave(filename = "Figures/dbCompare_pairedLengths_WoLr1_Conterminator_Exhaustive.jpeg",
       dpi = "retina", units = "in", width = 2.5, height = 4.5)

#----------------------------------------------------------#
# Import RS210 fasta files and format into data frames
#----------------------------------------------------------#
##---------------------Conterminator---------------------##
dbRS210_Conterminator <- read.fasta("Supporting_files/Database_scrubbing/RS210/RS210_conterminator_T2T_GRCh38_pangenome.fna")
length(dbRS210_Conterminator) # 1973 --> NOTE: this is not a unique genome count

dbRS210_ConterminatorDf <- data.frame(gID = gsub("\\:.+","",getName(dbRS210_Conterminator)),
                                      startStop = gsub("^G.+\\:","",getName(dbRS210_Conterminator))) %>%
  separate("startStop", c("start","stop"),
           sep = "-", remove = TRUE, convert = TRUE) %>%
  mutate(length = stop-start+1)
length(unique(dbRS210_ConterminatorDf$gID)) # 463 --> this is the number of unique genomes affected

##---------------------Exhaustive---------------------##
dbRS210_Exhaustive <- read.fasta("Supporting_files/Database_scrubbing/RS210/rs210_hg38.t2t.pangenome_contamination.fna")
length(dbRS210_Exhaustive) # 8671

dbRS210_ExhaustiveDf <- data.frame(gID = gsub("\\:.+","",getName(dbRS210_Exhaustive)),
                                   startStop = gsub("^G.+\\:","",getName(dbRS210_Exhaustive))) %>%
  separate("startStop", c("start","stop"),
           sep = "-", remove = TRUE, convert = TRUE) %>%
  mutate(length = stop-start+1)
length(unique(dbRS210_ExhaustiveDf$gID)) # 1106 --> this is the number of unique genomes affected

##---------------------Compare genomes---------------------##
sum(dbRS210_ConterminatorDf$gID %in% dbRS210_ExhaustiveDf$gID) # 1963
sum(dbRS210_ExhaustiveDf$gID %in% dbRS210_ConterminatorDf$gID) # 2844

dbRS210_totalUnique <- unique(c(dbRS210_ConterminatorDf$gID,
                                dbRS210_ExhaustiveDf$gID))
length(dbRS210_totalUnique) # 1116 --> 3.76%

dbRS210_List <- list(
  Conterminator = unique(dbRS210_ConterminatorDf$gID),
  Exhaustive = unique(dbRS210_ExhaustiveDf$gID)
)

require(ggvenn)
ggvenn(dbRS210_List, text_size = 3.5, show_percentage = FALSE) +
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF"))
ggsave("Figures/dbCompare_Venn_gIDs_RS210_Conterminator_Exhaustive.jpeg",
       dpi = "retina",units = "in", width = 3, height = 2)

##---------------------Compare distributions of lengths---------------------##
require(ggridges)
dbRS210_BothDf <- rbind(cbind(dbRS210_ConterminatorDf,method="Conterminator"),
                        cbind(dbRS210_ExhaustiveDf,method="Exhaustive")) %>%
  mutate(length = as.numeric(length))

## Unpaired comparisons
# Summarise lengths across genomes
dbRS210_BothDf_summarised <- dbRS210_BothDf %>%
  group_by(gID, method) %>%
  summarise(lengthComb = sum(length)) %>%
  arrange(gID)
# Calc unpaired Wilcoxon
dbRS210_BothDf_summarised %>% as.data.frame() %>% wilcox_test(lengthComb ~ method)
# Plot using ggridges
dbRS210_BothDf_summarised %>%
  ggplot(aes(x=log10(lengthComb),y=method,fill=stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01,
                               jittered_points = TRUE,
                               position = position_points_jitter(width = 0.05, height = 0),
                               point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  scale_fill_viridis_c(name = "Genome length", option = "C", alpha=0.7) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(x = "log10(length of contamination\nin affected genomes)",
       y = "")
ggsave(filename = "Figures/dbCompare_lengths_RS210_Conterminator_Exhaustive.jpeg",
       dpi = "retina", units = "in", width = 4, height = 3)

## Paired comparisons
dbRS210_Overlapping_gIDs <- dbRS210_BothDf_summarised$gID[duplicated(dbRS210_BothDf_summarised$gID)]
dbRS210_BothDf_summarised_paired <- dbRS210_BothDf_summarised %>%
  filter(gID %in% dbRS210_Overlapping_gIDs)
# Calc paired Wilcoxon
dbRS210_BothDf_summarised_paired %>%
  as.data.frame() %>%
  mutate(method = factor(method)) %>%
  wilcox_test(lengthComb ~ method, paired = TRUE) -> dbPairedWilcoxonRS210
# Plot paired differences
dbRS210_BothDf_summarised_paired %>%
  mutate(loglengthComb = log10(lengthComb)) %>%
  ggplot(aes(x=method,y=loglengthComb)) +
  geom_boxplot(aes(fill=method),notch = TRUE) +
  geom_point(size = 0.1, alpha=0.2) +
  geom_line(aes(group = gID),alpha=0.2,color="gray") +
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF")) +
  theme_pubr() + theme(legend.position = "none") +
  rotate_x_text(30) +
  theme(text = element_text(size = 20)) +
  labs(x = "", y = "log10(contamination length)") +
  # stat_compare_means(paired = TRUE)
  stat_pvalue_manual(dbPairedWilcoxonRS210, label = "p = {p}", y.position = 5.3, size = 5)
ggsave(filename = "Figures/dbCompare_pairedLengths_RS210_Conterminator_Exhaustive.jpeg",
       dpi = "retina", units = "in", width = 3, height = 6)


