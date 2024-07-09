library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)

general_celltype_gene_exp <- 
  read.csv("/Users/f006f9w/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/epithelial_celltype_pseudobulk_counts.csv", 
           row.names = 1)
metadata <- 
  read.csv("/Users/f006f9w/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/hmo_metadata_final.csv", row.names = 1)

## Remove the Non-secretor donors BM06 and BM13
metadata <- subset(metadata, !(donor %in% c('BM06', 'BM13')))

celltype_gene_exp <- na.omit(general_celltype_gene_exp)
celltype_gene_exp <- t(celltype_gene_exp)
celltype_gene_exp <- celltype_gene_exp[,rownames(metadata)]

# remove small samples
celltype_gene_exp <- celltype_gene_exp[,metadata$Cell_Number > 10]
metadata = metadata[colnames(celltype_gene_exp),]

all_hmos <- colnames(metadata)
all_hmos <- all_hmos[43:44]
# Determine the total number of columns
total_columns <- length(all_hmos)
#all_hmos <- c(all_hmos[1:19], all_hmos[(total_columns-2):total_columns])

scale(metadata$time_post_partum_days)

for( hmo in all_hmos){ 
  metadata[,hmo]<-scale(metadata[,hmo], center = TRUE, scale = TRUE)
}

for( i in c('LC1')){#unique(metadata$Epithelial.Cell.Subclusters)){
  for(hmo in all_hmos) {
    epi <- celltype_gene_exp[,metadata$Epithelial.Cell.Subclusters == 'LC1']
    
    epi <- na.omit(epi)
    
    # remove very late samples and low cell number samples
    epi_meta = metadata[colnames(epi),]
    epi <- epi[,epi_meta$time_post_partum_days < 400]
    epi_meta = epi_meta[colnames(epi),]
    epi <- epi[,epi_meta$Cell_Number > 10]
    epi_meta = epi_meta[colnames(epi),]
    
    dim(epi)
    
    formula_str <- paste("~donor+time_post_partum_days+",hmo)
    reduced_formula_str <- paste("~donor+time_post_partum_days")
    dds <- DESeqDataSetFromMatrix(countData = epi, colData = epi_meta, design = as.formula(formula_str))
    dds_lrt_time <- DESeq(dds, test="LRT", reduced = as.formula(reduced_formula_str))
    res<-results(dds_lrt_time)
    path <- paste("/Users/f006f9w/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/TypeIvsTypeII_final/",sub(" ", "_", i), sep = "")
    dir.create(path)
    write.csv(res,paste(path, "/", hmo, ".csv", sep = ""),quote=FALSE,row.names =TRUE)
  }
}


norm_mat <- estimateSizeFactors(dds_lrt_time)
norm_mat <- counts(norm_mat, normalized=TRUE)
saveRDS(norm_mat, "norm_mat_Cycling_Lactocytes.rds")

### scatter plots ###
library(patchwork)
hmo_genes <- c("B3GNT3")
#transporters <- c("SLC2A9", "SLC6A14","SLC02B1", "SLC30A2", "SLC30A8")
### Plotting the genes

# Create an empty list to store plots
plot_list <- list()

# Loop through each gene
for (gene in hmo_genes) {
  p <- plotCounts(dds_lrt_time, gene = gene, intgroup = "Type_I", returnData = TRUE)
  
  # Filter time_PP_weeks based on conditions
  time_PP_weeks <- metadata$time_post_partum_weeks[metadata$Epithelial.Cell.Subclusters == 'LC1']
  
  current_plot <- ggplot(p, aes(x = Type_I, y = count, color = time_PP_weeks)) +
    geom_point() +
    geom_smooth(method = "loess", color = "black", se = TRUE) +  # Add regression line with confidence interval
    labs(x = "Type_I HMO con.", y = "Normalized Counts") +
    ggtitle(paste("Gene", gene)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold.italic", size = 8),  # Bold and italicize the title
          axis.title = element_text(size = 8)) +  # Reduce the font size of axis labels
    theme(legend.position = "none")  # Remove legend
  
  # Store the plot in the list
  plot_list[[gene]] <- current_plot
}

# Combine all plots into a single plot
combined_plot <- wrap_plots(plot_list, ncol = 1)
combined_plot <- combined_plot + plot_layout(guides = "collect")
ggsave("combined_plot.pdf", combined_plot, width = 6, height = 7, units = "cm")


# Display the combined plot
print(combined_plot)

ggsave("hmogenes.pdf", combined_plot, device = 'pdf')

################################################################################

#mean_DSLNH = mean(metadata$DSLNH)

##Plot the HMOs

# Assuming 'metadata' is your data frame containing HMO data
hmo_column_names <- c("DSLNH",  "FDSLNH"  , "DFLNH"   , "FLNH"   ,  "DSLNT" ,   "LNH"   ,   "DFLNT"  ,  "LSTc",  
"LSTb","LNFP.III", "LNFP.II" , "LNFP.I" ,"LNnT" ,"LNT" , "X6.SL" ,"X3.SL"  ,  "DFLac"   , "X3FL"   ,  "X2.FL" )

pdf("combined_plots.pdf", width = 8, height = 12)  # Adjust width and height as needed
dev.off()

hist(metadata$DSLNT)


for( i in c('LC1')){ #unique(metadata$Epithelial.Cell.Subclusters)){
  for(hmo in all_hmos) {
    epi <- celltype_gene_exp[,metadata$Epithelial.Cell.Subclusters == i]
    
    epi <- na.omit(epi)
    
    # remove very late samples and low cell number samples
    epi_meta = metadata[colnames(epi),]
    epi <- epi[,epi_meta$time_post_partum_days < 400]
    epi_meta = epi_meta[colnames(epi),]
    epi <- epi[,epi_meta$Cell_Number > 10]
    epi_meta = epi_meta[colnames(epi),]
    
    dim(epi)
    
    pca_data <- prcomp(t(epi))
    pca_covariates <- as.data.frame(pca_data$x)
    formula_str <- "~ X6.SL + PC1 + PC2 + ..."
    reduced_formula_str <- "~ PC1 + PC2 + ..."
    dds <- DESeqDataSetFromMatrix(countData = epi, colData = epi_meta, design = as.formula(formula_str))
    dds_lrt_time <- DESeq(dds, test = "LRT", reduced = as.formula(reduced_formula_str))
    res <- results(dds_lrt_time)
    path <- paste("/Users/f006f9w/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/Test/",sub(" ", "_", i), sep = "")
    dir.create(path)
    write.csv(res,paste(path, "/", hmo, ".csv", sep = ""),quote=FALSE,row.names =TRUE)
  }
}


#Plotting HMO concentrations versus time post partum 
hmo_columns <- c(colnames(metadata[1:19]))
melted_df <- reshape2::melt(metadata, id.vars = "time_post_partum_days", measure.vars = hmo_columns)

ggplot(melted_df, aes(x = time_post_partum_days, y = value)) +
  geom_point() +
  facet_wrap(~ variable, scales = "free_y") +
  geom_vline(xintercept = 19, linetype = "dashed", color = "red") +
  labs(x = "Days Postpartum", y = "Concentration", title = "Concentration of HMOs against Days Postpartum")

