#BubblePlot

library(phyloseq) #v1.50.0
library(ggpubr) #v0.6.0
library(ggsci) #v3.2.0 
library(ggplot2) #v3.5.1
library(dplyr) #v1.1.4
library(fantaxtic) 


load("./algas_ps-gg2nb.RData")

algas_ps<- prune_taxa(taxa_sums(algas_ps) > 10, algas_ps)

#Filtro por baja varianza
otu_mat <- otu_table(algas_ps)
sd_taxa <- apply(otu_mat, 1, sd)
threshold_varianza <- quantile(sd_taxa, 0.1)
algas_ps <- prune_taxa(sd_taxa > threshold_varianza, algas_ps)



algas_phylum <- tax_glom(algas_ps, "Phylum")
algas_family <- tax_glom(algas_ps, "Family")
algas_class <- tax_glom(algas_ps, "Class")


phylum_ms <- merge_samples(algas_phylum, "Macroalgal_group", fun = sum)
family_ms <- merge_samples(algas_family, "Macroalgal_group", fun = sum)
class_ms <- merge_samples(algas_class, "Macroalgal_group", fun = sum)

############################# __ Macroalgal Group __ ##########################


#Phylum

top15 <- get_top_taxa(physeq_obj = phylum_ms, n = 15, relative = T,
                      discard_other = F, other_label = "Non-abundant Phylum")



topsito15 <- psmelt(top15)


colores_pasteles <- c("#E57373", "#81C784", "#64B5F6", "#FFB74D", "#9575CD",
                      "#4DB6AC", "#FF8A65", "#A1887F", "#FFD54F", "#90CAF9",
                      "#FFCC80", "#BA68C8", "#81D4FA", "#FFF176", "#80CBC4", "red")


topsito15 <- topsito15 %>%
  group_by(Sample) %>%
  mutate(Abundancia_total_muestra = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Abundancia_relativa = (Abundance / Abundancia_total_muestra) * 100) %>%
  group_by(Phylum) %>%
  mutate(Abundancia_relativa_total = sum(Abundancia_relativa)) %>%
  ungroup() %>%
  arrange(desc(Abundancia_relativa_total)) %>%
  mutate(Phylum = factor(Phylum, levels = unique(Phylum)))


ggplot(topsito15, aes(x = Sample, y = Phylum)) +
  geom_point(aes(size = Abundancia_relativa, fill = Sample), alpha = 0.75, shape = 21) +
  scale_fill_manual(values = rev(colores_pasteles)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size=13),
        axis.title.y= element_text(size = 15)) +
  ylab("Relative abundance (%)") +
  xlab("")


## Class


top15 <- get_top_taxa(physeq_obj = class_ms, n = 15, relative = T,
                      discard_other = F, other_label = "Non-abundant Class")



topsito15 <- psmelt(top15)


colores_pasteles <- c("#E57373", "#81C784", "#64B5F6", "#FFB74D", "#9575CD",
                      "#4DB6AC", "#FF8A65", "#A1887F", "#FFD54F", "#90CAF9",
                      "#FFCC80", "#BA68C8", "#81D4FA", "#FFF176", "#80CBC4", "red")


topsito15 <- topsito15 %>%
  group_by(Sample) %>%
  mutate(Abundancia_total_muestra = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Abundancia_relativa = (Abundance / Abundancia_total_muestra) * 100) %>%
  group_by(Class) %>%
  mutate(Abundancia_relativa_total = sum(Abundancia_relativa)) %>%
  ungroup() %>%
  arrange(desc(Abundancia_relativa_total)) %>%
  mutate(Class = factor(Class, levels = unique(Class)))


ggplot(topsito15, aes(x = Sample, y = Class)) +
  geom_point(aes(size = Abundancia_relativa, fill = Sample), alpha = 0.75, shape = 21) +
  scale_fill_manual(values = rev(colores_pasteles)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size=13),
        axis.title.y= element_text(size = 15)) +
  ylab("Relative abundance (%)") +
  xlab("")


## Family

top15 <- get_top_taxa(physeq_obj = family_ms, n = 15, relative = T,
                      discard_other = F, other_label = "Non-abundant Families")



topsito15 <- psmelt(top15)


colores_pasteles <- c("#E57373", "#81C784", "#64B5F6", "#FFB74D", "#9575CD",
                      "#4DB6AC", "#FF8A65", "#A1887F", "#FFD54F", "#90CAF9",
                      "#FFCC80", "#BA68C8", "#81D4FA", "#FFF176", "#80CBC4", "red")


topsito15 <- topsito15 %>%
  group_by(Sample) %>%
  mutate(Abundancia_total_muestra = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Abundancia_relativa = (Abundance / Abundancia_total_muestra) * 100) %>%
  group_by(Family) %>%
  mutate(Abundancia_relativa_total = sum(Abundancia_relativa)) %>%
  ungroup() %>%
  arrange(desc(Abundancia_relativa_total)) %>%
  mutate(Family = factor(Family, levels = unique(Family)))


ggplot(topsito15, aes(x = Sample, y = Family)) +
  geom_point(aes(size = Abundancia_relativa, fill = Sample), alpha = 0.75, shape = 21) +
  scale_fill_manual(values = rev(colores_pasteles)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size=13),
        axis.title.y= element_text(size = 15)) +
  ylab("Relative abundance (%)") +
  xlab("")


############################# __ Country __ ##########################


phylum_rs <- merge_samples(algas_phylum, "Country", fun = sum)
family_rs <- merge_samples(algas_family, "Country", fun = sum)
class_rs <- merge_samples(algas_class, "Country", fun = sum)


#Phylum

top15 <- get_top_taxa(physeq_obj = phylum_rs, n = 15, relative = T,
                      discard_other = F, other_label = "Z-Non-abundant Phylum")



topsito15 <- psmelt(top15)


topsito15 <- topsito15 %>%
  group_by(Sample) %>%
  mutate(Abundancia_total_muestra = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Abundancia_relativa = (Abundance / Abundancia_total_muestra) * 100) %>%
  group_by(Phylum) %>%
  mutate(Abundancia_relativa_total = sum(Abundancia_relativa)) %>%
  ungroup() %>%
  arrange(desc(Abundancia_relativa_total)) %>%
  mutate(Phylum = gsub("p__", "", Phylum), # Eliminar el prefijo "p__"
         Phylum = gsub("d__Bacteria", " Z-Unclassified bacteria", Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels = sort(unique(Phylum))))



ggplot(topsito15, aes(x = Sample, y = Phylum)) +
  geom_point(aes(size = Abundancia_relativa), alpha = 0.75, shape = 21, fill= "#28B463")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size=10),
        axis.title.y= element_text(size = 15)) +
  ylab("") +
  xlab("")


#Class

top15 <- get_top_taxa(physeq_obj = class_rs, n = 15, relative = T,
                      discard_other = F, other_label = "Z-Non-abundant Class")


topsito15 <- psmelt(top15)

topsito15 <- topsito15 %>%
  group_by(Sample) %>%
  mutate(Abundancia_total_muestra = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Abundancia_relativa = (Abundance / Abundancia_total_muestra) * 100) %>%
  group_by(Class) %>%
  mutate(Abundancia_relativa_total = sum(Abundancia_relativa)) %>%
  ungroup() %>%
  arrange(desc(Abundancia_relativa_total)) %>%
  mutate(Class = gsub("c__", "", Class),
         Class = gsub("d__Bacteria", " Z-Unclassified bacteria", Class)) %>%
  mutate(Class = factor(Class, levels = sort(unique(Class))))


ggplot(topsito15, aes(x = Sample, y = Class)) +
  geom_point(aes(size = Abundancia_relativa), alpha = 0.75, shape = 21, fill="red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size=10),
        axis.title.y= element_text(size = 15)) +
  ylab("") +
  xlab("")


## Family
algas_family <- tax_glom(algas_ps, "Family")
family_rs <- merge_samples(algas_family, "Country", fun = sum)

top15 <- get_top_taxa(physeq_obj = family_rs, n = 15, relative = T,
                      discard_other = F, other_label = "Non-abundant Families")



topsito15 <- psmelt(top15)


topsito15 <- topsito15 %>%
  group_by(Sample) %>%
  mutate(Abundancia_total_muestra = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Abundancia_relativa = (Abundance / Abundancia_total_muestra) * 100) %>%
  group_by(Family) %>%
  mutate(Abundancia_relativa_total = sum(Abundancia_relativa)) %>%
  ungroup() %>%
  arrange(desc(Abundancia_relativa_total)) %>%
  mutate(Family = gsub("f__", "", Family),
         Family = gsub("d__Bacteria", " Z-Unclassified bacteria", Family)) %>%
  mutate(Family = factor(Family, levels = sort(unique(Family))))


ggplot(topsito15, aes(x = Sample, y = Family)) +
  geom_point(aes(size = Abundancia_relativa), alpha = 0.75, shape = 21, fill="blue4") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size=10),
        axis.title.y= element_text(size = 15)) +
  ylab("") +
  xlab("")
