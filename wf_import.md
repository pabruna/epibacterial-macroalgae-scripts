---
title: "Workflow-central"
author: "Benjamin Leyton-Carcaman"
date: "2023-09-30"
output: html_document
---

### Librerias y funciones

```{r setup, include=TRUE}
knitr::opts_chunk$set(
	echo = TRUE,
	message = FALSE,
	warning = FALSE
)
library(dplyr)
library(tidyverse)
library(phyloseq)


tss_normalizer <- function(table) {
  df <- read.table(table, sep="\t", header=TRUE)
  df[, -1] <- sapply(df[, -1], as.numeric)
  
  # Normalizar por columna (muestra)
  total_per_sample <- colSums(df[, -1])
  df_normalizado <- df %>%
    mutate(across(-1, ~ (. / total_per_sample) * 1000000))
  
  return(df_normalizado)
}




```

### Import metadata

The metadata was downloaded directly from NCBI, specifically from the BioSamples of each sample. All these metadata files were consolidated into a single file named `metadata.tsv`.

```{r cars}

metadata_full <- read.table("metadata.tsv", header = TRUE, sep = "\t")

rm(list=setdiff(ls(), c("metadata_full", "tss_normalizer")))

```


### Import feature-tables

Each feature table was generated using DADA2 for each study. Sequence length was estimated based on the quality of the reads and is summarized in Table S1. Each feature table was validated with a rarefaction curve.



```{r pressure, echo=TRUE}

features_tables <- list.files(path = "TSV_feature-tables/",
                              pattern = "*.tsv",
                              full.names = FALSE)


features_tables<- sub("\\.tsv$", "", features_tables)

tables <- list()

for (i in features_tables) {
  tables[[i]] <- tss_normalizer(paste0("./TSV_feature-tables/",i,".tsv"))
}


rm(list=setdiff(ls(), c("metadata_full", "tss_normalizer", "tables")))
```


### Import Taxonomy
```{r}
taxonomys<-list()

files_tax<-list.files(path = "./TSV_taxonomys/", pattern = "*.tsv",
                      full.names = FALSE)
files_tax<-sub("\\.tsv$", "", files_tax)

for (i in files_tax){
  
  temp <- read.table(paste0("./TSV_taxonomys//", i, ".tsv"), sep = "\t", header = TRUE)
  temp <- temp[ , !(names(temp) %in% c("Confidence"))]
  
  split_data <- strsplit(as.character(temp$Taxon), ";")
  new_data <- do.call(rbind, split_data)
  
  colnames(new_data) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # Agregamos la columna "Feature ID" al dataframe resultante
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
dim(tables[["10-feature-table"]])
dim(taxonomys[["study10-taxonomy"]])
dim(mergeados[["study10"]])

rm(list=setdiff(ls(), c("metadata_full", "tss_normalizer", "tables", "mergeados", "taxonomys")))
```


```{r}

full_full<-bind_rows(mergeados)



#Comprobar weas no deseadas
unique(full_full$Phylum)[grepl("Mitochondria", unique(full_full$Phylum))]
unique(full_full$Order)[grepl(c("Chloroplast", "Mitochondria"), unique(full_full$Order))]
unique(full_full$Class)[grepl(c("Chloroplast", "Mitochondria"), unique(full_full$Class))]
unique(full_full$Family)[grepl(c("Chloroplast", "Mitochondria"), unique(full_full$Family))]
unique(full_full$Genus)[grepl(c("Chloroplast", "Mitochondria"), unique(full_full$Genus))]

full_full %>%
  filter(Domain %in% c("d__Bacteria", "d__Archaea" )) %>% #elima eucariotas
  mutate_all(~str_replace_all(., " ", "")) -> full_full

full_full <- full_full %>%
  replace(is.na(full_full), 0) %>%
  mutate_at(vars(9:889), as.numeric)



## eliminar muestras que sumen 0, es deir que no tengan ASVs, debido a los filtros

full_full <- full_full[, !(sapply(full_full, is.numeric) & colSums(full_full == 0, na.rm = TRUE) == nrow(full_full)), drop = FALSE]

full_full_sumarizado <- full_full %>%
  mutate_all(~ ifelse(. == "", "not_determined", .)) %>%
  group_by(Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(across('SRR16971116':'SRR16090778', sum, na.rm = TRUE), .groups = "drop")

full_full_sumarizado %>%
  mutate(Feature_ID = paste0("ASV", row_number())) %>%
  relocate(Feature_ID, .before = everything()) -> gg2_nb

#rm(list=setdiff(ls(), c("metadata", "gg2_nb")))

save(full_full_sumarizado, full_full, gg2_nb, file = "tmp.RData")
```


# Crear objeto phyloseq

```{r}

taxonomia <- as.data.frame(gg2_nb[,1:8])
#write.table(taxonomia, "./Phyloseq_input/tax-gg2-nb.txt", row.names = FALSE, quote = FALSE, sep = "\t")
rownames(taxonomia)<-taxonomia$Feature_ID
tax <- as.matrix(taxonomia[,-1])


asv_table <- as.data.frame(gg2_nb[,c(1,9:889)])
#write.table(asv_table, "./Phyloseq_input/table-gg2-nb.txt", row.names = FALSE, quote = FALSE, sep = "\t")
rownames(asv_table) <- asv_table$Feature_ID
otu <- as.matrix(asv_table[-1])

samples_ok <-colnames(otu)

sam <- data.frame(metadata_full) %>%
  filter(sample_id %in% samples_ok)
#write.table(sam, "./Phyloseq_input/metadata-gg2-nb.txt", row.names = FALSE, quote = FALSE, sep = "\t")

rownames(sam) <- sam$sample_id
sam <- sam[,-1]


algas_ps<-phyloseq(otu_table(otu, taxa_are_rows=TRUE),
                  tax_table(tax),
                  sample_data(data.frame(sam)))


save(algas_ps, file = "algas_ps-gg2nb.RData")
```

# Export to MA


```{r}
####### Export to microbialAnalyst

tax<-gg2_nb[,1:8]
colnames(tax)[1:2]<-c("#TAXONOMY","Kingdom")

otu<-gg2_nb[, c(1,9:889)]
colnames(otu)[1]<-"#NAME"

sam <- data.frame(metadata_full) %>%
  filter(sample_id %in% samples_ok)
colnames(sam)[1]<-"#NAME"


write.table(tax, "./Microbial-Analyst_input/tax-gg2-nb.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(otu, "./Microbial-Analyst_input/otu-gg2-nb.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(sam, "./Microbial-Analyst_input/metadata-gg2-nb.txt", row.names = FALSE, quote = FALSE, sep = "\t")
```
