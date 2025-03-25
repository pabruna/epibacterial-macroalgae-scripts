# Script: Macroalgal Microbiome Diversity Analysis
# Updated: March 2025
# Description: Performs alpha diversity analysis of macroalgal-associated microbial communities across taxonomic levels (ASV, Family, Phylum). 
# Includes data filtering, diversity metric calculation, visualization, and statistical comparisons between macroalgal groups and host species.

# Key Features:
# 1. Data preprocessing:
#    - Filters low-abundance taxa (sum < 10)
#    - Removes low-variance features (bottom 10% quantile)
#    - Aggregates data at Phylum and Family levels

# 2. Alpha diversity analysis:
#    - Calculates Chao1 (richness), Shannon and inverse Simpson (evenness) indices
#    - Filters samples with Chao1 < 10
#    - Comparative analysis across:
#       * Macroalgal groups (Chlorophyta/Phaeophyta/Rhodophyta/Seawater)
#       * Individual host species (22 macroalgal species)

# 3. Visualization:
#    - Boxplot+jitter plots for group comparisons
#    - Custom color schemes for algal groups
#    - Publication-ready figure panels

# 4. Statistical analysis:
#    - Kruskal-Wallis tests with agricolae package
#    - Post-hoc grouping of significant differences

# Output: Comparative diversity plots and statistical results for:
# - ASV-level microbial communities
# - Family-level aggregations  
# - Phylum-level aggregations


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








algas_asv <- algas_ps
algas_phylum <- tax_glom(algas_ps, "Phylum")
algas_family <- tax_glom(algas_ps, "Family")


alfa_diversidad<-function(objeto_phyloseq){
  
  temp1<-richness(objeto_phyloseq, index = c("chao1"))
  temp6<-evenness(objeto_phyloseq, index = c("pielou", "simpson") ) 
  df<-cbind(temp1, temp6)
  df %>%
    mutate(sample_id=rownames(.))%>%
    select(sample_id, everything())-> df
  rownames(df)<-1:length(rownames(df))
  colnames(df) <- c("sample_id","Chao1", "Shannon", "inv_Simpson")
  
  metadata1 <- data.frame(sample_data(objeto_phyloseq))
  metadata1$sample_id <- rownames(metadata1)
  rownames(metadata1) <- 1:length(rownames(metadata1))
  
  df2 <- merge(metadata1, df, by.x = "sample_id",
              by.y = "sample_id")
  return(df2) 
  
}


# metricas de diversidad como tal.
alfa_phylum<-alfa_diversidad(algas_phylum)
alfa_family <- alfa_diversidad(algas_family)
alfa_asv <- alfa_diversidad(algas_ps)


# hay muestras pobres, por lo que es mejor excluirlas

alfa_asv %>% filter(Chao1 >= 10) -> alfa_asv
alfa_phylum %>% filter(Chao1 >=10) -> alfa_phylum
alfa_family %>% filter(Chao1 >= 10) -> alfa_family

# a continuacion son las comparacion de macroalgas


colores <- c("#008000", "#654321","#8B0000", "#ADD8E6")
indices <- c("Chao1", "Shannon", "inv_Simpson")


# a nivel de ASV
plots <- list()

for (index in indices){
  plot_name <- paste0("ASV_", index,"_Macroalgal-group", ".png")
  alfa_asv %>%
    ggplot(aes_string(x = "Macroalgal_group",
                      y = index,
                      color = "Macroalgal_group",
                      fill = "Macroalgal_group")) +
    geom_boxplot(alpha=0.5, show.legend =TRUE,
                 outlier.colour = NULL,
                 outlier.fill = NULL,
                 outlier.shape = NULL,) +
    geom_jitter(width = 0.2, show.legend = FALSE, shape=21, colour= "black", stroke=0.5)+
    theme_bw() +
    xlab("") +
    scale_color_manual(values = colores) +
    scale_fill_manual(values = colores) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.title = element_text(size=18, colour = "black"),
          legend.title = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size=13)) -> plots[[index]]
}

ggarrange(plots[[1]], plots[[2]], plots[[3]], labels = c("A", "B", "C"), ncol = 3,
          common.legend = TRUE, legend = "right",
          font.label = list(size = 20, color = "black", face = "bold", family = NULL))


# a nivel de Phylum

plots <- list()

for (index in indices){
  plot_name <- paste0("ASV_", index,"_Macroalgal-group", ".png")
  alfa_phylum %>%
    ggplot(aes_string(x = "Macroalgal_group",
                      y = index,
                      color = "Macroalgal_group",
                      fill = "Macroalgal_group")) +
    geom_boxplot(alpha=0.5, show.legend =TRUE,
                 outlier.colour = NULL,
                 outlier.fill = NULL,
                 outlier.shape = NULL,) +
    geom_jitter(width = 0.2, show.legend = FALSE, shape=21, colour= "black", stroke=0.5)+
    theme_bw() +
    xlab("") +
    scale_color_manual(values = colores) +
    scale_fill_manual(values = colores) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.title = element_text(size=18, colour = "black"),
          legend.title = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size=13)) -> plots[[index]]
}

ggarrange(plots[[1]], plots[[2]], plots[[3]], labels = c("A", "B", "C"), ncol = 3,
          common.legend = TRUE, legend = "right",
          font.label = list(size = 20, color = "black", face = "bold", family = NULL))



# a nivel de Familia
plots <- list()

for (index in indices){
  plot_name <- paste0("ASV_", index,"_Macroalgal-group", ".png")
  alfa_family %>%
    ggplot(aes_string(x = "Macroalgal_group",
                      y = index,
                      color = "Macroalgal_group",
                      fill = "Macroalgal_group")) +
    geom_boxplot(alpha=0.5, show.legend =TRUE,
                 outlier.colour = NULL,
                 outlier.fill = NULL,
                 outlier.shape = NULL,) +
    geom_jitter(width = 0.2, show.legend = FALSE, shape=21, colour= "black", stroke=0.5)+
    theme_bw() +
    xlab("") +
    scale_color_manual(values = colores) +
    scale_fill_manual(values = colores) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.title = element_text(size=18, colour = "black"),
          legend.title = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size=13)) -> plots[[index]]
}

ggarrange(plots[[1]], plots[[2]], plots[[3]], labels = c("A", "B", "C"), ncol = 3,
          common.legend = TRUE, legend = "right",
          font.label = list(size = 20, color = "black", face = "bold", family = NULL))


#Estadistica Macroalgas group

library(agricolae)

kruskal(alfa_asv$Chao1, alfa_asv$Macroalgal_group)$groups
kruskal(alfa_asv$Shannon, alfa_asv$Macroalgal_group)$groups
kruskal(alfa_asv$inv_Simpson, alfa_asv$Macroalgal_group)$groups

kruskal(alfa_family$Chao1, alfa_family$Macroalgal_group)$groups
kruskal(alfa_family$Shannon, alfa_family$Macroalgal_group)$groups
kruskal(alfa_family$inv_Simpson, alfa_family$Macroalgal_group)$groups

kruskal(alfa_phylum$Chao1, alfa_phylum$Macroalgal_group)$groups
kruskal(alfa_phylum$Shannon, alfa_phylum$Macroalgal_group)$groups
kruskal(alfa_phylum$inv_Simpson, alfa_phylum$Macroalgal_group)$groups


#############################################################################################################

resumen <- data.frame(table(alfa_asv$Source, alfa_asv$Macroalgal_group)) %>%
  filter(Freq != 0)

resumen2<-data.frame(table(resumen$Var2)) #tenemos 22 especies con 17 pardas

colores2 <- c("#008000", "#008000", "#654321", "#654321", "#654321", "#654321", 
             "#654321", "#654321", "#654321", "#654321", "#654321", "#654321", 
             "#654321", "#654321", "#654321", "#654321", "#8B0000", "#8B0000", 
             "#8B0000", "#ADD8E6")


# a nivel de ASV
plots <- list()
alfa_asv$Source <- trimws(alfa_asv$Source)


for (index in indices){
  #plot_name <- paste0("ASV_", index,"_source", ".png")
  alfa_asv %>%
    mutate(Source = case_when(
      Source == "Nereocystis leutkeana laboratory" ~ "Nereocystis luetkeana",
      Source == "Nereocystis leutkeana laboratory natural" ~ "Nereocystis luetkeana",
      TRUE ~ Source
    )) %>%
    mutate(Macroalgal_group = factor(Macroalgal_group, levels = c("Chlorophyta", "Phaeophyta", "Rhodophyta", "Seawater")),
           Source_Ordered = interaction(Macroalgal_group, Source, lex.order = TRUE)) %>%
    
    ggplot(aes(x = Source_Ordered, y = .data[[index]],
               color = Source_Ordered, fill = Source_Ordered)) +
    geom_boxplot(alpha=0.5, show.legend =FALSE,
                 outlier.colour = NULL,
                 outlier.fill = NULL,
                 outlier.shape = NULL,
                 ) +
    theme_bw() +
    xlab("") +
    scale_color_manual(values = colores2) +
    scale_fill_manual(values = colores2) +
    theme(axis.text.x = element_text(face = c(rep("italic", length.out = 19), "plain"),
                                     angle = 90,
                                     hjust = 1, vjust = 0.5, size = 11,
                                     colour = "black"),
          axis.text.y = element_text(size = 11, colour = "black"),
          axis.title = element_text(size=15, colour = "black"),
          legend.title = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size=13)) -> plots[[index]]
}

ggarrange(plots[[1]], plots[[2]], plots[[3]], labels = c("A", "B", "C"), ncol = 3,
          common.legend = TRUE, legend = "right",
          font.label = list(size = 20, color = "black", face = "bold", family = NULL))



# Family level
plots <- list()
alfa_family$Source <- trimws(alfa_family$Source)


for (index in indices){
  #plot_name <- paste0("ASV_", index,"_source", ".png")
  alfa_family %>%
    mutate(Source = case_when(
      Source == "Nereocystis leutkeana laboratory" ~ "Nereocystis luetkeana",
      Source == "Nereocystis leutkeana laboratory natural" ~ "Nereocystis luetkeana",
      TRUE ~ Source
    )) %>%
    mutate(Macroalgal_group = factor(Macroalgal_group, levels = c("Chlorophyta", "Phaeophyta", "Rhodophyta", "Seawater")),
           Source_Ordered = interaction(Macroalgal_group, Source, lex.order = TRUE)) %>%
    
    ggplot(aes(x = Source_Ordered, y = .data[[index]],
               color = Source_Ordered, fill = Source_Ordered)) +
    geom_boxplot(alpha=0.5, show.legend =FALSE,
                 outlier.colour = NULL,
                 outlier.fill = NULL,
                 outlier.shape = NULL,
    ) +
    theme_bw() +
    xlab("") +
    scale_color_manual(values = colores2) +
    scale_fill_manual(values = colores2) +
    theme(axis.text.x = element_text(face = c(rep("italic", length.out = 19), "plain"),
                                     angle = 90,
                                     hjust = 1, vjust = 0.5, size = 11,
                                     colour = "black"),
          axis.text.y = element_text(size = 11, colour = "black"),
          axis.title = element_text(size=15, colour = "black"),
          legend.title = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size=13)) -> plots[[index]]
}

ggarrange(plots[[1]], plots[[2]], plots[[3]], labels = c("A", "B", "C"), ncol = 3,
          common.legend = TRUE, legend = "right",
          font.label = list(size = 20, color = "black", face = "bold", family = NULL))


# Phylum level


colores2 <- c("#008000", "#008000", "#654321", "#654321", "#654321", "#654321", 
              "#654321", "#654321", "#654321", "#654321", "#654321", "#654321", 
              "#654321", "#654321", "#654321", "#654321", "#8B0000", "#8B0000", 
               "#ADD8E6")
plots <- list()
alfa_phylum$Source <- trimws(alfa_phylum$Source)


for (index in indices){
  #plot_name <- paste0("ASV_", index,"_source", ".png")
  alfa_phylum %>%
    mutate(Source = case_when(
      Source == "Nereocystis leutkeana laboratory" ~ "Nereocystis luetkeana",
      Source == "Nereocystis leutkeana laboratory natural" ~ "Nereocystis luetkeana",
      TRUE ~ Source
    )) %>%
    mutate(Macroalgal_group = factor(Macroalgal_group, levels = c("Chlorophyta", "Phaeophyta", "Rhodophyta", "Seawater")),
           Source_Ordered = interaction(Macroalgal_group, Source, lex.order = TRUE)) %>%
    
    ggplot(aes(x = Source_Ordered, y = .data[[index]],
               color = Source_Ordered, fill = Source_Ordered)) +
    geom_boxplot(alpha=0.5, show.legend =FALSE,
                 outlier.colour = NULL,
                 outlier.fill = NULL,
                 outlier.shape = NULL,
    ) +
    theme_bw() +
    xlab("") +
    scale_color_manual(values = colores2) +
    scale_fill_manual(values = colores2) +
    theme(axis.text.x = element_text(face = c(rep("italic", length.out = 18), "plain"),
                                     angle = 90,
                                     hjust = 1, vjust = 0.5, size = 11,
                                     colour = "black"),
          axis.text.y = element_text(size = 11, colour = "black"),
          axis.title = element_text(size=15, colour = "black"),
          legend.title = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size=13)) -> plots[[index]]
}

ggarrange(plots[[1]], plots[[2]], plots[[3]], labels = c("A", "B", "C"), ncol = 3,
          common.legend = TRUE, legend = "right",
          font.label = list(size = 20, color = "black", face = "bold", family = NULL))

# Stats
library(agricolae) #v1.3-7

kruskal(alfa_asv$Chao1, alfa_asv$Source)$groups
kruskal(alfa_asv$Shannon, alfa_asv$Source)$groups
kruskal(alfa_asv$inv_Simpson, alfa_asv$Source)$groups

kruskal(alfa_family$Chao1, alfa_family$Source)$groups
kruskal(alfa_family$Shannon, alfa_family$Source)$groups
kruskal(alfa_family$inv_Simpson, alfa_family$Source)$groups

kruskal(alfa_phylum$Chao1, alfa_phylum$Source)$groups
kruskal(alfa_phylum$Shannon, alfa_phylum$Source)$groups
kruskal(alfa_phylum$inv_Simpson, alfa_phylum$Source)$groups