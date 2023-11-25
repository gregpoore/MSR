# 01-taxa-filtering-pipeline.R
# Author: Greg Poore
# Date: Oct 7, 2023
# Purposes:
# - Define taxa filtering pipeline

#-------------------------------#
# Load dependencies
require(doMC)
require(plyr)
require(dplyr)
require(tidyr)
require(tibble)
require(reshape2)
require(phyloseq)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import RS210 lineages and format
#----------------------------------------------------------#
# Load RS210 lineage
taxRS210 <- read.csv("Supporting_data/RS210_lineages.csv",
                     stringsAsFactors = FALSE, header = FALSE, col.names = c("OGU","Taxa")) %>%
  mutate(Taxa = gsub("[d|p|c|o|f|g|s]__","",Taxa))

# Format
taxRS210_ff <- taxRS210 %>%
  separate("Taxa", c("domain","phylum","class","order","family","genus","species"),
           sep = "; ", remove = TRUE) %>%
  mutate(species = gsub(" ATCC.+$","",species)) %>%
  mutate(species = gsub(" DSM.+$","",species)) %>%
  mutate(species = gsub(" sp\\..+$","",species)) %>%
  mutate(species = gsub(" NBRC.+$","",species)) %>%
  mutate(species = gsub("  "," ",species)) %>%
  mutate(species = gsub(" ","_",species))

#----------------------------------------------------------#
# Import WIS data summarized at species level
# See README in Supporting_data/Weizmann_data
#----------------------------------------------------------#

wzRDS <- readRDS("Supporting_data/Weizmann_data/bacterial_fungal_all_rank_count_hit_list.rds")
psWz_allRank <- wzRDS$all_rank
psWz_phylum <- wzRDS$phylum_phy
psWz_class <- wzRDS$class_phy
psWz_order <- wzRDS$order_phy
psWz_family <- wzRDS$family_phy
psWz_genus <- wzRDS$genus_phy
psWz_species <- wzRDS$species_phy

# Subset to non-control samples to compare overlap
psWz_allRank_Bio <- subset_samples(psWz_allRank, type.detail %in% c("normal", "nat", "tumor"))
psWz_phylum_Bio <- subset_samples(psWz_phylum, type.detail %in% c("normal", "nat", "tumor"))
psWz_class_Bio <- subset_samples(psWz_class, type.detail %in% c("normal", "nat", "tumor"))
psWz_order_Bio <- subset_samples(psWz_order, type.detail %in% c("normal", "nat", "tumor"))
psWz_family_Bio <- subset_samples(psWz_family, type.detail %in% c("normal", "nat", "tumor"))
psWz_genus_Bio <- subset_samples(psWz_genus, type.detail %in% c("normal", "nat", "tumor"))
psWz_species_Bio <- subset_samples(psWz_species, type.detail %in% c("normal", "nat", "tumor"))

wisSpeciesBioDf <- data.frame(tax_table(psWz_species_Bio)) %>%
  filter(!grepl("Unknown|other",species, ignore.case = TRUE)) %>%
  mutate(species = gsub("Escherichia/Shigella","Escherichia",species)) %>%
  mutate(species = gsub(" ","_",species))
dim(wisSpeciesBioDf) # 594   7

wisSpeciesBioKnown <- wisSpeciesBioDf$species

## Subset to tumor-only samples
psWz_species_Tumor <- subset_samples(psWz_species, type.detail %in% c("tumor"))
psWz_species_Tumor_OTU <- data.frame(otu_table(psWz_species_Tumor))
psWz_species_Tumor_OTU_nonzero <- psWz_species_Tumor_OTU[rowSums(psWz_species_Tumor_OTU)!=0,] # drops 172
psWz_species_Tumor_TAX <- data.frame(tax_table(psWz_species_Tumor))
psWz_species_Tumor_TAX_nonzero <- psWz_species_Tumor_TAX[rownames(psWz_species_Tumor_OTU_nonzero),]

wisSpeciesTumorDf <- psWz_species_Tumor_TAX_nonzero %>%
  filter(!grepl("Unknown|other",species, ignore.case = TRUE)) %>%
  mutate(species = gsub("Escherichia/Shigella","Escherichia",species)) %>%
  mutate(species = gsub(" ","_",species))
dim(wisSpeciesTumorDf) # 525   7

wisSpeciesTumorKnown <- wisSpeciesTumorDf$species
#----------------------------------------------------------#
# Intersect RS210 and WIS species
#----------------------------------------------------------#

taxRS210_ff_WIS_Bio <- taxRS210_ff[which(taxRS210_ff$species %in% wisSpeciesBioKnown),] %>%
  arrange(desc(domain)) # fungi first
dim(taxRS210_ff_WIS_Bio) # 813   8

taxRS210_ff_WIS_Tumor <- taxRS210_ff[which(taxRS210_ff$species %in% wisSpeciesTumorKnown),]  %>%
  arrange(desc(domain)) # fungi first
dim(taxRS210_ff_WIS_Tumor) # 807   8

taxRS210_ff_WIS_Tumor %>% count(domain)
taxRS210_ff_WIS_Tumor %>% pull(species) %>% unique()

#----------------------------------------------------------#
# Import external taxa tables
#----------------------------------------------------------#
# Metagenomes from multiple body sites (UNITN)
pasolli2019Cell <- read.csv("Supporting_data/pasolli2019-Cell-TableS4.csv",
                        stringsAsFactors = FALSE)
# UHGG gut metagenomes
almeida2020NBT <- read.csv("Supporting_data/almeida2020-NBT-Supp-Table2.csv",
                           stringsAsFactors = FALSE)
# List of 1513 bacterial pathogens
bartlett2022Micro <- read.csv("Supporting_data/bartlett2022-microbiology-bacterial-pathogens.csv",
                              stringsAsFactors = FALSE)
# List of ~1k cultivated gut organisms
rajilic2014FEMS <- read.csv("Supporting_data/Rajilic-Stojanovic2014-FEMS-MicroRev-TablesS1-S3.csv",
                            stringsAsFactors = FALSE)

## Format -- Pasolli / UNITN
pasolli2019Cell_EstimatedTaxa <- pasolli2019Cell %>% 
  filter(level_estimated_taxonomy == "Species") %>%
  select(estimated_taxonomy) %>%
  mutate(estimated_taxonomy = gsub("k__|p__|c__|o__|f__|g__|s__","",estimated_taxonomy)) %>%
  separate("estimated_taxonomy", c("domain","phylum","class","order","family","genus","species"),
           sep = "\\|", remove = TRUE) %>%
  mutate(species = gsub("_ATCC.+$","",species)) %>%
  mutate(species = gsub("_DSM.+$","",species)) %>%
  mutate(species = gsub("_sp_.+$","",species)) %>%
  mutate(species = gsub("_NBRC.+$","",species)) %>%
  mutate(species = gsub("__","_",species)) %>%
  mutate(species = gsub(" ","_",species)) %>%
  distinct()
dim(pasolli2019Cell_EstimatedTaxa) # 748 7

## Format -- Almeida / UHGG
almeida2020NBT_GTDBTaxa <- almeida2020NBT %>% 
  select(tax_GTDB) %>%
  mutate(tax_GTDB = gsub("d__|p__|c__|o__|f__|g__|s__","",tax_GTDB)) %>%
  separate("tax_GTDB", c("domain","phylum","class","order","family","genus","species"),
           sep = ";", remove = TRUE) %>%
  mutate(species = gsub(" ATCC.+$","",species)) %>%
  mutate(species = gsub(" DSM.+$","",species)) %>%
  mutate(species = gsub(" sp[0-9]+$","",species)) %>%
  mutate(species = gsub(" NBRC.+$","",species)) %>%
  mutate(species = gsub("  "," ",species)) %>%
  mutate(species = gsub(" ","_",species)) %>%
  filter(species != "") %>%
  distinct()
dim(almeida2020NBT_GTDBTaxa) # 1297    7

## Format -- Barlett / Pathogens
bartlett2022Micro_Taxa <- bartlett2022Micro %>% 
  filter(status == "established") %>%
  rename(domain=superkingdom) %>%
  mutate(species = paste0(genus,"_",species)) %>%
  select(-year, -status, -reference) %>%
  distinct()
dim(bartlett2022Micro_Taxa) # 1110 7

## Format -- Rajilic / Cultivated gut microbes
rajilic2014FEMS_Taxa <- rajilic2014FEMS %>%
  mutate(species = gsub(" ","_",species)) %>%
  select(-reference) %>%
  distinct()
dim(rajilic2014FEMS_Taxa) # 1008 5

#----------------------------------------------------------#
# Combine external taxa tables:
# - UNITN
# - UHGG
# - Pathogens
# - WIS
#----------------------------------------------------------#

combTaxaSpecies <- unique(c(pasolli2019Cell_EstimatedTaxa$species,
                          almeida2020NBT_GTDBTaxa$species,
                          bartlett2022Micro_Taxa$species,
                          wisSpeciesBioDf$species))
length(combTaxaSpecies) # 2862

combTaxaSpecies <- data.frame(species = combTaxaSpecies) %>%
  mutate(UNITN = species %in% pasolli2019Cell_EstimatedTaxa$species) %>%
  mutate(UHGG = species %in% almeida2020NBT_GTDBTaxa$species) %>%
  mutate(PATH = species %in% bartlett2022Micro_Taxa$species) %>%
  mutate(WIS = species %in% wisSpeciesBioDf$species)

## Intersect with RS210
taxRS210_ff_combTaxaSpecies <- taxRS210_ff[which(taxRS210_ff$species %in% combTaxaSpecies$species),] %>%
  arrange(desc(domain)) %>% # fungi first
  left_join(combTaxaSpecies, by  = "species")
dim(taxRS210_ff_combTaxaSpecies) # 3321    12

#----------------------------------------------------------#
# Add Zebra coverages
# NOTE: The databases above do *not* have viruses in them,
# so we will add back in viruses with >=50% coverage.
# NOTE: PhiX can infect E coli but is removed here
# NOTE: Cutibacterium acnes passes and is a biological spp but is removed here
#----------------------------------------------------------#

zebraTCGA_All <- read.csv("Input_data/zebra/13722_pangenome_13767_transcript_removal.zebra_coverages.tsv",
                          sep = "\t", stringsAsFactors = FALSE, row.names = 1) %>% 
  rename(species=strain) %>%
  mutate(species = gsub(" ","_",species))
dim(zebraTCGA_All) # 24120     4

# NOTE: Woltka table merging has a default threshold of ~1e-07
# for each feature to be detected (per sam file), but this threshold is not 
# used for Zebra. Thus, these lowly prevalent features were removed from the Zebra file
# by intersecting it with the Woltka-merged table (from script '02-RS210-pangenome-processing.R')
load("Interim_data/rs210PanFinalOGUs_WGS_RNA_Transcript_13Oct23.RData")
zebraTCGA_AllFilt <- zebraTCGA_All %>%
  rownames_to_column("OGU") %>%
  filter(OGU %in% rs210PanFinalOGUs) %>%
  column_to_rownames("OGU")
dim(zebraTCGA_AllFilt) # 23381     4

## Add coverage info to filtered RS210 table and filter
taxRS210_ff_combTaxaSpeciesZebra <- taxRS210_ff_combTaxaSpecies %>%
  mutate(cov_All = zebraTCGA_AllFilt[OGU,"coverage_ratio"]) %>%
  filter(!is.na(cov_All))

quantile(taxRS210_ff_combTaxaSpeciesZebra$cov_All)
# 0%                    25%          50%          75%         100% 
# 0.0001118008 0.0168275013 0.0900160541 0.5040450789 1.0000000000 

# Save supp table
taxRS210_ff_combTaxaSpeciesZebra %>%
  write.csv("Supp_tables/taxRS210_ff_combTaxaSpeciesZebra.csv",
            row.names = FALSE)

## Combine filtered non-viral and viral data

# Create viral table with known coverages
taxRS210_ff_virusCov <- taxRS210_ff %>%
  filter(domain == "Viruses") %>%
  filter(species != "Escherichia_virus_phiX174") %>% # Remove PhiX
  mutate(cov_All = zebraTCGA_AllFilt[OGU,"coverage_ratio"]) %>%
  filter(!is.na(cov_All)) %>%
  mutate(UNITN=FALSE, UHGG=FALSE, PATH=FALSE, WIS=FALSE)

# Check how many viruses or non-viruses survive the coverage cutoffs
taxRS210_ff_virusHighCov <- taxRS210_ff_virusCov %>%
  filter(cov_All >= 0.50) %>% # Can change this value
  mutate(UNITN=FALSE, UHGG=FALSE, PATH=FALSE, WIS=FALSE)
dim(taxRS210_ff_virusHighCov)

# Combine viral and non-viral data with known coverages
taxRS210_ff_combTaxaSpeciesZebraWithViruses <- rbind(taxRS210_ff_combTaxaSpeciesZebra,
                                                     taxRS210_ff_virusCov)

# Require 90% non-viral coverage and 90% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_9090 <- taxRS210_ff_combTaxaSpeciesZebraWithViruses %>%
  filter(species != "Cutibacterium_acnes") %>%
  filter( ( (domain!="Viruses") & (cov_All >= 0.90)) | ((domain=="Viruses") & (cov_All >= 0.90)) )
dim(taxRS210_ff_combTaxaSpeciesZebraWithViruses_9090) # 508  13

# Require 75% non-viral coverage and 75% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_7575 <- taxRS210_ff_combTaxaSpeciesZebraWithViruses %>%
  filter(species != "Cutibacterium_acnes") %>%
  filter( ( (domain!="Viruses") & (cov_All >= 0.75)) | ((domain=="Viruses") & (cov_All >= 0.75)) )
dim(taxRS210_ff_combTaxaSpeciesZebraWithViruses_7575) # 765  13

# Require 50% non-viral coverage and 80% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_5080 <- taxRS210_ff_combTaxaSpeciesZebraWithViruses %>%
  filter(species != "Cutibacterium_acnes") %>%
  filter( ( (domain!="Viruses") & (cov_All >= 0.50)) | ((domain=="Viruses") & (cov_All >= 0.80)) )
dim(taxRS210_ff_combTaxaSpeciesZebraWithViruses_5080) # 935  13

# Require 50% non-viral coverage and 75% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_5075 <- taxRS210_ff_combTaxaSpeciesZebraWithViruses %>%
  filter(species != "Cutibacterium_acnes") %>%
  filter( ( (domain!="Viruses") & (cov_All >= 0.50)) | ((domain=="Viruses") & (cov_All >= 0.75)) )
dim(taxRS210_ff_combTaxaSpeciesZebraWithViruses_5075) # 956  13

# Require 50% non-viral coverage and 50% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_5050 <- taxRS210_ff_combTaxaSpeciesZebraWithViruses %>%
  filter(species != "Cutibacterium_acnes") %>%
  filter( ( (domain!="Viruses") & (cov_All >= 0.50)) | ((domain=="Viruses") & (cov_All >= 0.50)) )
dim(taxRS210_ff_combTaxaSpeciesZebraWithViruses_5050) # 1189 13


## Identify number of unique species
# Require 90% non-viral coverage and 90% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_9090US <- taxRS210_ff_combTaxaSpeciesZebraWithViruses_9090 %>%
  pull(species) %>% unique()
length(taxRS210_ff_combTaxaSpeciesZebraWithViruses_9090US) # 249

# Require 75% non-viral coverage and 75% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_7575US <- taxRS210_ff_combTaxaSpeciesZebraWithViruses_7575 %>%
  pull(species) %>% unique()
length(taxRS210_ff_combTaxaSpeciesZebraWithViruses_7575US) # 381

# Require 50% non-viral coverage and 80% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_5080US <- taxRS210_ff_combTaxaSpeciesZebraWithViruses_5080 %>%
  pull(species) %>% unique()
length(taxRS210_ff_combTaxaSpeciesZebraWithViruses_5080US) # 435

# Require 50% non-viral coverage and 75% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_5075US <- taxRS210_ff_combTaxaSpeciesZebraWithViruses_5075 %>%
  pull(species) %>% unique()
length(taxRS210_ff_combTaxaSpeciesZebraWithViruses_5075US) # 456

# Require 50% non-viral coverage and 50% viral coverage
taxRS210_ff_combTaxaSpeciesZebraWithViruses_5050US <- taxRS210_ff_combTaxaSpeciesZebraWithViruses_5050 %>%
  pull(species) %>% unique()
length(taxRS210_ff_combTaxaSpeciesZebraWithViruses_5050US) # 689

# save(taxRS210_ff_combTaxaSpecies,
#      taxRS210_ff_combTaxaSpeciesZebra,
#      taxRS210_ff_combTaxaSpeciesZebraWithViruses,
#      taxRS210_ff_combTaxaSpeciesZebraWithViruses_9090,
#      taxRS210_ff_combTaxaSpeciesZebraWithViruses_7575,
#      taxRS210_ff_combTaxaSpeciesZebraWithViruses_5080,
#      taxRS210_ff_combTaxaSpeciesZebraWithViruses_5075,
#      taxRS210_ff_combTaxaSpeciesZebraWithViruses_5050,
#      taxRS210_ff_virusCov,
#      taxRS210_ff,
#      zebraTCGA_AllFilt,
#      file = "Interim_data/taxa_filtering_pipeline_13Oct23.RData")

# Save for supp tables
taxRS210_ff_combTaxaSpeciesZebraWithViruses_5050US %>%
  write.csv("Supp_tables/rs210clean_5050_unique_species.csv",
            row.names = FALSE)
taxRS210_ff_combTaxaSpeciesZebraWithViruses_7575US %>%
  write.csv("Supp_tables/rs210clean_7575_unique_species.csv",
            row.names = FALSE)
taxRS210_ff_combTaxaSpeciesZebraWithViruses_9090US %>%
  write.csv("Supp_tables/rs210clean_9090_unique_species.csv",
            row.names = FALSE)




