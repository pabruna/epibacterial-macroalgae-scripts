# Script: Microbiome Shared Features Analysis Using UpSet Plots
# Updated: March 2025
# Description: Visualizes shared and unique microbial features across macroalgal 
# groups and geographic regions using UpSet plots

# Key Features:
# 1. Data preprocessing:
#    - Filters low-abundance taxa (sum < 10)
#    - Removes low-variance features (bottom 10% quantile)
#    - Aggregates data at three taxonomic levels:
#       * Phylum-level
#       * Family-level
#       * ASV-level (most granular)

# 2. Comparative analyses:
#    A) By Macroalgal Group:
#       - Chlorophyta (green algae)
#       - Phaeophyta (brown algae)
#       - Rhodophyta (red algae)
#       - Seawater controls
#    B) By Geographic Region:
#       - USA, UK, Canada, China, Chile, Australia

# 3. Visualization:
#    - UpSet plots showing intersections between groups
#    - Binary presence/absence threshold (≥4 counts)
#    - Consistent sorting for comparability
#    - Three taxonomic resolution levels for comprehensive analysis

# Output:
# - UpSet plots showing feature sharing patterns:
#   1. Among macroalgal groups at each taxonomic level
#   2. Among geographic regions at each taxonomic level
# - Reveals core vs. unique microbial communities

library(UpSetR) #v1.4.0
library(ComplexUpset) #v1.3.3
library(microbiome) #v1.28.0
library(phyloseq) #v1.50.0
library(ggplot2) #v3.5.1
library(dplyr) #v1.1.4
library(car) #v3.1-3
library(ggpubr) #v0.6.0
library(ggsci) #v3.2.0 

load("./algas_ps-gg2nb.RData")

algas_ps<- prune_taxa(taxa_sums(algas_ps) > 10, algas_ps)

#Filtro por baja varianza
otu_mat <- otu_table(algas_ps)
sd_taxa <- apply(otu_mat, 1, sd)
threshold_varianza <- quantile(sd_taxa, 0.1)
algas_ps <- prune_taxa(sd_taxa > threshold_varianza, algas_ps)

##

algas_asv <- algas_ps
algas_phylum <- tax_glom(algas_ps, "Phylum")
algas_family <- tax_glom(algas_ps, "Family")

## Upset a nivel de Phylum en Grupos de Macroalgas

phylum_ms <- merge_samples(algas_phylum, "Macroalgal_group", fun = sum)
phylum_ms_t <-as.data.frame(t(phylum_ms@otu_table))
phylum_ms_t[phylum_ms_t >= 4] <- 1

upset(data = phylum_ms_t, intersect = c("Seawater","Chlorophyta",
                                        "Rhodophyta", "Phaeophyta"),
      name = "Number of Phyla",  sort_sets=FALSE
)

## Upset a nivel de Familia en Grupos de Macroalgas

family_ms <- merge_samples(algas_family, "Macroalgal_group", fun = sum)
family_ms_t <-as.data.frame(t(family_ms@otu_table))
family_ms_t[family_ms_t >= 4] <- 1

upset(data = family_ms_t, intersect = c("Seawater","Chlorophyta",
                                                          "Rhodophyta", "Phaeophyta"),
      name = "Number of Families", sort_sets=FALSE
)


## Upset a nivel de ASV en grupos de Macroalgas

asv_ms <- merge_samples(algas_asv, "Macroalgal_group", fun = sum)
asv_ms_t <-as.data.frame(t(asv_ms@otu_table))
asv_ms_t[asv_ms_t >= 4] <- 1

upset(data = asv_ms_t, intersect = c("Seawater","Chlorophyta",
                                        "Rhodophyta", "Phaeophyta"),
      name = "Number of ", sort_sets=FALSE
)


##############################__REGION__##################################

## Upset a nivel de Phylum en regiones

phylum_rs <- merge_samples(algas_phylum, "Country", fun = sum)
phylum_rs_t <-as.data.frame(t(phylum_rs@otu_table))
phylum_rs_t[phylum_rs_t >= 4] <- 1

upset(data = phylum_rs_t, intersect = c("USA" ,"United Kingdom",
                                        "Canada", "China","Chile", "Australia"),
      name = "Number of Phyla",  sort_sets=FALSE
)

## Upset a nivel de Familia en regiones

family_rs <- merge_samples(algas_family, "Country", fun = sum)
family_rs_t <-as.data.frame(t(family_rs@otu_table))
family_rs_t[family_rs_t >= 4] <- 1

upset(data = family_rs_t, intersect = c("USA" ,"United Kingdom",
                                        "Canada", "China","Chile", "Australia"),
      name = "Number of Families",  sort_sets=FALSE
)


## Upset a nivel de ASV en regiones

asv_rs <- merge_samples(algas_asv, "Country", fun = sum)
asv_rs_t <-as.data.frame(t(asv_rs@otu_table))
asv_rs_t[asv_rs_t >= 4] <- 1

upset(data = asv_rs_t, intersect = c("USA" ,"United Kingdom",
                                        "Canada", "China","Chile", "Australia"),
      name = "Number of ASVs",  sort_sets=FALSE
)

