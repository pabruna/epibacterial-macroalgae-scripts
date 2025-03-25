# Script: Microbial Community Composition Bubble Plots  
# Updated: March 2025  
# Description: Visualizes relative abundance of microbial taxa across macroalgal groups and geographic regions using bubble plots  

# Key Features:  
# 1. Data preprocessing:  
#    - Filters low-abundance taxa (sum < 10)  
#    - Removes low-variance features (bottom 10% quantile)  
#    - Aggregates data at three taxonomic levels:  
#       * Phylum  
#       * Class  
#       * Family  

# 2. Comparative analyses:  
#    A) By Macroalgal Group:  
#       - Chlorophyta (green algae)  
#       - Phaeophyta (brown algae)  
#       - Rhodophyta (red algae)  
#       - Seawater controls  
#    B) By Country:  
#       - USA, UK, Canada, China, Chile, Australia  

# 3. Visualization:  
#    - Bubble plots showing top 15 most abundant taxa  
#    - Bubble size represents relative abundance (%)  
#    - Color-coded by sample group  
#    - Clean taxonomic labels (removes prefixes like "p__")  
#    - Consistent sorting by total relative abundance  

# Output:  
# - Comparative bubble plots at different taxonomic resolutions:  
#   1. Phylum-level composition  
#   2. Class-level composition  
#   3. Family-level composition  
# - Reveals dominant microbial groups in each sample type  



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


algas_phylum <- tax_glom(algas_ps, "Phylum")
algas_family <- tax_glom(algas_ps, "Family")

################################# Regiones 16S ##############################


#Plot Phylum
set.seed(2024)
p_ord<-ordinate(algas_phylum, "PCoA", "bray")


plot_ordination(algas_phylum, p_ord, type = "samples", color = "Region") +
  stat_ellipse(size=0.7, show.legend = FALSE) +
  geom_point(show.legend = TRUE, shape=21, colour= "black", stroke=1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size=14, colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=13)) + 
  scale_color_npg() +
  guides(color = guide_legend(override.aes = list(size = 5))) -> p1


#Plot  Family
set.seed(2024)
p_ord2<-ordinate(algas_family, "PCoA", "bray")


plot_ordination(algas_family, p_ord2, type = "samples", color = "Region") +
  stat_ellipse(size=0.7, show.legend = FALSE) +
  geom_point(show.legend = TRUE, shape=21, colour= "black", stroke=1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size=14, colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=13)) + 
  scale_color_npg() +
  guides(color = guide_legend(override.aes = list(size = 5))) -> p2

ggarrange(p1, p2, labels = c("Phylum level", "Family level"),
          common.legend = TRUE, legend = "right",
          font.label = list(size = 15, color = "black", face = "bold", family = NULL))



############################      Paises     ################################


#Plot Phylum

plot_ordination(algas_phylum, p_ord, type = "samples", color = "Country") +
  stat_ellipse(size=0.7, show.legend = FALSE) +
  geom_point(show.legend = TRUE, shape=21, colour= "black", stroke=1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size=14, colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=13)) + 
  scale_color_jco() +
  guides(color = guide_legend(override.aes = list(size = 5))) -> a1


#Plot  Family

plot_ordination(algas_family, p_ord2, type = "samples", color = "Country") +
  stat_ellipse(size=0.7, show.legend = FALSE) +
  geom_point(show.legend = TRUE, shape=21, colour= "black", stroke=1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size=14, colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=13)) + 
  scale_color_jco() +
  guides(color = guide_legend(override.aes = list(size = 5))) -> a2

ggarrange(a1, a2, labels = c("Phylum level", "Family level"),
          common.legend = TRUE, legend = "right",
          font.label = list(size = 15, color = "black", face = "bold", family = NULL))

############################ Macroalgal_Group ################################


colores <- c("#008000", "#654321","#8B0000", "#ADD8E6")

#Plot Phylum

plot_ordination(algas_phylum, p_ord, type = "samples", color = "Macroalgal_group") +
  stat_ellipse(size=0.7, show.legend = FALSE) +
  geom_point(show.legend = TRUE, shape=21, colour= "black", stroke=1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size=14, colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=13)) + 
  scale_color_manual(values = colores) +
  guides(color = guide_legend(override.aes = list(size = 5))) -> m1


#Plot  Family

plot_ordination(algas_family, p_ord, type = "samples", color = "Macroalgal_group") +
  stat_ellipse(size=0.7, show.legend = FALSE) +
  geom_point(show.legend = TRUE, shape=21, colour= "black", stroke=1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size=14, colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=13)) + 
  scale_color_manual(values = colores) +
  guides(color = guide_legend(override.aes = list(size = 5))) -> m2


ggarrange(m1, m2, labels = c("Phylum level", "Family level"),
          common.legend = TRUE, legend = "right",
          font.label = list(size = 15, color = "black", face = "bold", family = NULL))



############################ Macroalgal_Group2 ################################

## Pardas vs Seawater

ps_phylum_sub <- subset_samples(algas_phylum,
                                Macroalgal_group %in% c("Phaeophyta", "Seawater"))


ps_family_sub <- subset_samples(algas_family,
                                Macroalgal_group %in% c("Phaeophyta", "Seawater"))


ord1<- ordinate(ps_phylum_sub, "PCoA", "bray")
ord2<- ordinate(ps_family_sub, "PCoA", "bray")



colores <- c( "#654321", "#ADD8E6")

#Plot Phylum

plot_ordination(ps_phylum_sub, ord1, type = "samples", color = "Macroalgal_group") +
  stat_ellipse(size=0.7, show.legend = FALSE) +
  geom_point(show.legend = TRUE, shape=21, colour= "black", stroke=1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size=14, colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=13)) + 
  scale_color_manual(values = colores) +
  guides(color = guide_legend(override.aes = list(size = 5))) -> b1


#Plot  Family

plot_ordination(ps_family_sub, ord2, type = "samples", color = "Macroalgal_group") +
  stat_ellipse(size=0.7, show.legend = FALSE) +
  geom_point(show.legend = TRUE, shape=21, colour= "black", stroke=1) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size=14, colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=13)) + 
  scale_color_manual(values = colores) +
  guides(color = guide_legend(override.aes = list(size = 5))) -> b2


ggarrange(b1, b2, labels = c("Phylum level", "Family level"),
          common.legend = TRUE, legend = "right",
          font.label = list(size = 15, color = "black", face = "bold", family = NULL))


###### ESTADISTICA ######



lista_phyloseq <- list(algas_phylum = algas_phylum,
                       algas_family = algas_family,
                       ps_phylum_sub = ps_phylum_sub,
                       ps_family_sub = ps_family_sub)

resultados <- data.frame(Objeto = character(), PValue = numeric(), stringsAsFactors = FALSE)

for (nombre in names(lista_phyloseq)) {
  ps_actual <- lista_phyloseq[[nombre]]
  distancias <- phyloseq::distance(ps_actual, method = "bray")
  datos_muestra <- as(sample_data(ps_actual), "data.frame")
  permanova_resultados <- vegan::adonis2(distancias ~ Country, data = datos_muestra, permutations = 999)
  p_value <- permanova_resultados$`Pr(>F)`[1]
  resultados <- rbind(resultados, data.frame(Objeto = nombre, PValue = p_value))
}

print(resultados)
