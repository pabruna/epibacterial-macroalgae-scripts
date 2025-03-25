
# Script: Microbiome Studies Processing and Merging
# Updated March 2025
# Description: This script integrates and normalizes multiple microbiome studies, processing feature tables,
# taxonomies, and metadata to build a phyloseq object. Includes TSS normalization, chloroplast/eukaryote filtering,
# taxonomy unification, and data aggregation. Output is a phyloseq object ready for advanced microbial analyses

# Workflow:
# 1. TSS normalization of DADA2-generated feature tables (conversion to ppm)
# 2. Metadata loading from consolidated NCBI BioSamples
# 3. Taxonomy processing (7-level standardization)
# 4. Study-specific feature-taxonomy merging (ASV reannotation)
# 5. Contaminant filtering (chloroplasts, eukaryotes)
# 6. ASV aggregation by full taxonomic lineage
# 7. Phyloseq object construction (OTU + TAX + METADATA)

# Output: "algas_ps-gg2nb.RData" - Phyloseq object for downstream analyses


###############################################################
library(dplyr) #v1.1.4
library(tidyverse) #v2.0.0
library(phyloseq) #v1.50.0


#normalizer function

tss_normalizer <- function(table) {
  df <- read.table(table, sep = "\t", header = TRUE)
  df[, -1] <- sapply(df[, -1], as.numeric)
  
  total_per_sample <- colSums(df[, -1])
  
  df_normalizado <- df %>%
    mutate(across(-1, ~ (. / total_per_sample[cur_column()]) * 1e6))
  
  return(df_normalizado)
}

### Import metadata
#The metadata was downloaded directly from NCBI, specifically from the BioSamples of each sample. All these metadata files were consolidated into a single file named `metadata.tsv`.

metadata_full <- read.table("./data/metadata.tsv", header = TRUE, sep = "\t")
rm(list=setdiff(ls(), c("metadata_full", "tss_normalizer")))

### Import feature-tables
# Each feature table was generated using DADA2 for each study. 
# Sequence length was estimated based on the quality of the reads and is summarized in Table S1. 
# Each feature table was validated with a rarefaction curve.


features_tables <- list.files(path = "./data/TSV_feature-tables/",
                              pattern = "*.tsv",
                              full.names = FALSE)

features_tables<- sub("\\.tsv$", "", features_tables)

tables <- list()

for (i in features_tables) {
  tables[[i]] <- tss_normalizer(paste0("./data/TSV_feature-tables/",i,".tsv"))
}

rm(list=setdiff(ls(), c("metadata_full", "tss_normalizer", "tables")))


### Import Taxonomy. 
#Each taxonomy table was generated using the Naïve Bayes classifier adjusted with specific primers in QIIME 2.

taxonomys <- list()
files_tax <- list.files(path = "./data/TSV_taxonomys/", pattern = "*.tsv", full.names = FALSE)
files_tax <- sub("\\.tsv$", "", files_tax)

for (i in files_tax) {
  temp <- read.table(paste0("./data/TSV_taxonomys/", i, ".tsv"), sep = "\t", header = TRUE)
  temp <- temp[, !(names(temp) %in% c("Confidence"))]
  
  split_data <- strsplit(as.character(temp$Taxon), ";")
  
  split_data <- lapply(split_data, function(x) {
    if (length(x) < 7) {
      x <- c(x, rep(NA, 7 - length(x)))
    } else if (length(x) > 7) {
      x <- x[1:7]
    }
    return(x)
  })
  
  new_data <- do.call(rbind, split_data)
  colnames(new_data) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  new_data <- cbind(Feature.ID = temp$Feature.ID, new_data)
  
  taxonomys[[i]] <- as.data.frame(new_data)
}

studies<-c("2", "4", "6", "7", "8", "9", "10", "15",
           "16", "17", "18", "19", "21", "22", "25", "27", "12")


mergeados <- list()
for (i in studies) {
  df1 <- taxonomys[[(paste0("study",i,"-taxonomy"))]]
  df2 <- tables[[(paste0(i,"-feature-table"))]]
  df1 <- rename(df1, "OTU_ID" = `Feature.ID`)
  merged <- df1 %>%
    right_join(df2, by = "OTU_ID") %>%
    mutate(Feature_ID = paste0("ASV",i,"_", row_number())) %>%
    select(Feature_ID, everything(), -OTU_ID)
  mergeados[[paste0("study",i)]] <- merged
}


rm(list=setdiff(ls(), c("metadata_full", "tss_normalizer", "tables", "mergeados", "taxonomys")))


#### Merge All!!

full_full<-bind_rows(mergeados) 



#unique(full_full$Phylum)[grepl("Chloroplast", unique(full_full$Phylum))]
#unique(full_full$Order)[grepl(c("Chloroplast"), unique(full_full$Order))]
#unique(full_full$Class)[grepl(c("Chloroplast"), unique(full_full$Class))]
#unique(full_full$Family)[grepl(c("Chloroplast"), unique(full_full$Family))]
#unique(full_full$Genus)[grepl(c("Chloroplast"), unique(full_full$Genus))]

full_full %>%
  filter(Domain %in% c("d__Bacteria", "d__Archaea" )) %>% 
  mutate_all(~str_replace_all(., " ", "")) -> full_full 


# Reemplazar NA con 0

full_full <- full_full %>%
  mutate(across(9:889, ~ ifelse(is.na(.), 0, .))) %>%
  mutate(across(9:889, as.numeric))
full_full <- full_full[, !(sapply(full_full, is.numeric) & colSums(full_full == 0, na.rm = TRUE) == nrow(full_full)), drop = FALSE]



full_full_sumarizado <- full_full %>%
  mutate(
    Domain = ifelse(is.na(Domain), "d__", Domain),
    Phylum = ifelse(is.na(Phylum), "p__", Phylum),
    Class = ifelse(is.na(Class), "c__", Class),
    Order = ifelse(is.na(Order), "o__", Order),
    Family = ifelse(is.na(Family), "f__", Family),
    Genus = ifelse(is.na(Genus), "g__", Genus),
    Species = ifelse(is.na(Species), "s__", Species)
  ) %>%
  mutate(full_tax = paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep = "|")) %>%
  group_by(full_tax) %>%
  summarise(across('SRR16971116':'SRR16090778', sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Domain = sapply(strsplit(full_tax, "\\|"), `[`, 1),
    Phylum = sapply(strsplit(full_tax, "\\|"), `[`, 2),
    Class = sapply(strsplit(full_tax, "\\|"), `[`, 3),
    Order = sapply(strsplit(full_tax, "\\|"), `[`, 4),
    Family = sapply(strsplit(full_tax, "\\|"), `[`, 5),
    Genus = sapply(strsplit(full_tax, "\\|"), `[`, 6),
    Species = sapply(strsplit(full_tax, "\\|"), `[`, 7)
  ) %>%
  select(-full_tax)

full_full_sumarizado %>%
  mutate(Feature_ID = paste0("ASV", row_number())) %>%
  relocate(Feature_ID,Domain, Phylum, Class, Order, Family, Genus,Species, .before = everything()) -> gg2_nb



# PHYLOSEQ OBJECT


taxonomia <- as.data.frame(gg2_nb[,1:8])
rownames(taxonomia)<-taxonomia$Feature_ID
tax <- as.matrix(taxonomia[,-1])


asv_table <- as.data.frame(gg2_nb[,c(1,9:889)])
rownames(asv_table) <- asv_table$Feature_ID
otu <- as.matrix(asv_table[-1])

samples_ok <-colnames(otu)

sam <- data.frame(metadata_full) %>%
  filter(sample_id %in% samples_ok)

rownames(sam) <- sam$sample_id
sam <- sam[,-1]


algas_ps<-phyloseq(otu_table(otu, taxa_are_rows=TRUE),
                   tax_table(tax),
                   sample_data(data.frame(sam)))


save(algas_ps, file = "algas_ps-gg2nb.RData")

