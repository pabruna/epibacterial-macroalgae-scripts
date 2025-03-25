# Script: Microbial Co-occurrence Network Analysis  
# Updated: March 2025  
# Description: Calculates pairwise microbial associations for network construction 
# using correlation methods  

# Key Features:  
# 1. Data Preparation:  
#    - Loads processed phyloseq object  
#    - Extracts and transposes OTU table for correlation analysis  

# 2. Network Analysis:  
#    - Implements pairwise correlation calculations between all ASVs  
#    - Supports both Pearson (linear) and Spearman (rank-based) methods  
#    - Progress tracking during computation  
#    - Stores complete correlation matrices with p-values  

# 3. Output:  
#    - Pearson correlation network (pearson_edges.csv)  
#    - Spearman correlation network (spearman_edges.csv)  
#    - Files contain: source node, target node, correlation coefficient, p-value  

# Note:  
# - Large datasets may require significant computation time (O(n^2) complexity)  
# - Results suitable for network visualization in tools like Gephi or Cytoscape  
# - Includes progress messages to track calculation status  

library(phyloseq)
library(dplyr)
library(tidyr)
library(SpiecEasi)
library(Hmisc)

load("algas_ps-gg2nb.RData")
ps<-algas_ps

otu_data <- t(as.matrix(otu_table(ps)))

# Función para calcular correlaciones y p-valores
correlation_network <- function(data, method) {
  n_taxa <- ncol(data)
  source <- character()
  target <- character()
  correlation <- numeric()
  p_value <- numeric()
  
  for (i in 1:(n_taxa - 1)) {
    # Mensaje de progreso
    cat(sprintf("Procesando taxón %d de %d\n", i, n_taxa))
    
    for (j in (i+1):n_taxa) {
      test_result <- cor.test(data[,i], data[,j], method=method, conf.level=0)
      ASV1 <- c(source, colnames(data)[i])
      ASV2 <- c(target, colnames(data)[j])
      correlation <- c(correlation, test_result$estimate)
      p_value <- c(p_value, test_result$p.value)
    }
  }
  return(data.frame(source, target, correlation, p_value))
}

# Red de correlación Pearson


pearson_network <- correlation_network(otu_data, "pearson")
warwrite.csv(pearson_network, "pearson_edges.csv", row.names=FALSE)

# Red de correlación Spearman
spearman_network <- correlation_network(otu_data, "spearman")
write.csv(spearman_network, "spearman_edges.csv", row.names=FALSE)
