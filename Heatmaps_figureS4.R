##### HeatMaps for HMO and genes ###
library(ComplexHeatmap)
library(tidyverse)
library("RColorBrewer")
library(readxl)

metadata <- 
  read.csv("/Users/f006f9w/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/hmo_metadata_final.csv", 
           row.names = 1)
norm_mat <- readRDS("/Users/f006f9w/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/Data/Secretory_Lactocytes/norm_mat_secretory_lactocytes.rds")

all_hmos <- colnames(metadata)
all_hmos <- all_hmos[43:44]
#all_hmos <- list(all_hmos[1:19], all_hmos[43:44])
#all_hmos <- list(all_hmos[5],all_hmos[12])

IPA <- read_excel("/Users/f006f9w/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/IPA.xls", 
                        col_names = TRUE, skip = 1)

for( hmo in all_hmos){ 
  metadata[,hmo]<-scale(metadata[,hmo], center = TRUE, scale = TRUE)
}

# Get the column names of normmat
column_names <- colnames(norm_mat)
metadata_LC1 <- metadata[column_names, , drop = FALSE]

directory_path <- "~/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/Data/Secretory_Lactocytes/Filtered_GeneLists/"
setwd(directory_path)
#files_in_directory = list(paste0(directory_path,"filter_",all_hmos[[1]],'.csv'), paste0(directory_path,'filter_',all_hmos[[2]],'.csv'))
files_in_directory <- list.files(directory_path, pattern = "^filter_.*\\.csv$")

Type <- IPA$`Type(s)`
left_ann_df <- data.frame(Type, row.names = IPA$ID)
unique_levels <- unique(left_ann_df$Type)
colors_list <- c("lightgray", "yellow2", "yellow4", "black", "skyblue", "blue4", "cyan2", "cyan4", "green2", "green4", "brown", "coral2", "purple4","orange4", "pink")
level_colors <- (length(unique_levels))
level_colors <- setNames(colors_list, unique_levels)

# Loop through each HMO file
for (hmo_file_name in files_in_directory) {
  
  hmo_name <- gsub(".*_(.*?)csv.*", "\\1", hmo_file_name)
  hmo_name <- substr(hmo_name, 1, nchar(hmo_name) - 1)
  hmo_name <- paste("Type_", hmo_name, sep = '')
  
  # Read the HMO data
  hmo <- read.csv(hmo_file_name)
  rownames(hmo) <- hmo['X'][,]
  
  # Calculate normalized matrix and Z-scores
  result <- tryCatch(
    {
      if (nrow(hmo) < 40) {
        next  # Skip this iteration if length is less than 40
      }
      mat <- norm_mat[rownames(hmo),]
      if (!is.null(mat)) {
        mat.z <- t(apply(mat, 1, scale))
        colnames(mat.z) <- colnames(mat)
        hmo_concentration <- as.vector(metadata_LC1[, colnames(metadata_LC1) == hmo_name])
        time_post_partum_days <- metadata_LC1$time_post_partum_days
        ann_df <- data.frame(time_post_partum_days, hmo_concentration)
        rownames(ann_df) <- rownames(metadata_LC1)
        colnames(ann_df) <- c("time_post_partum_days", "hmo_concentration")

        mat.z <- mat.z[, order(as.numeric(ann_df$hmo_concentration))]
        
        n = nrow(ann_df)
        row_font <- 5
        if (length(rownames(hmo)) < 60) {
          row_font <- 9
        }
        ## time post partum
        min_v_tpp = 0
        max_v_tpp = 400
        col_seq_tpp <- seq(min_v_tpp, max_v_tpp, 5)
        Var_tpp = circlize::colorRamp2(col_seq_tpp, rev(hcl.colors(length(col_seq_tpp),"Blues")))
        ## hmo_conc
     
        min_v_conc = min(ann_df$hmo_concentration)
        max_v_conc = max(ann_df$hmo_concentration)
        col_seq_conc <- seq(min_v_conc, max_v_conc, 0.5)
        Var_conc <- circlize::colorRamp2(col_seq_conc, rev(hcl.colors(length(col_seq_conc), "Rocket")))
        # Create the heatmap with dynamically adapted annotation column
        
        Type <- IPA$`Type(s)`
        left_ann_df <- data.frame(Type, row.names = IPA$ID)
        
        #common_genes <- intersect(rownames(left_ann_df), rownames(hmo))
        common_genes <- hmo[hmo$X %in% rownames(left_ann_df), ]
        common_genes <- common_genes[order(-common_genes$log2FoldChange), ]
        common_genes <- rbind(
          common_genes[1:25, ], 
        common_genes[(nrow(common_genes) - 24):nrow(common_genes), ])
        common_genes = common_genes$X
        
        left_ann_df <- left_ann_df[common_genes, , drop = FALSE]
        colnames(left_ann_df) <- c("Type_of_Molecule")
        temp_mat.z <- mat.z[row.names(left_ann_df), ]
        
        left_ann_df$Type_of_Molecule <- factor(left_ann_df$Type_of_Molecule, levels = unique_levels)
        Var_molecule <- level_colors[left_ann_df$Type_of_Molecule]
        ha <- rowAnnotation(Type_of_Molecule = left_ann_df$Type_of_Molecule,
                            col = list(Type_of_Molecule = Var_molecule),
                            show_annotation_name = FALSE)
        
        ht <- Heatmap(temp_mat.z,
                      name = "Z score",
                      cluster_columns = FALSE,
                      clustering_distance_rows = "manhattan",
                      clustering_method_rows = "average",
                      show_row_dend = TRUE,
                      show_column_dend = TRUE,
                      show_row_names = TRUE,
                      show_column_names = FALSE,
                      row_names_side = "left",
                      top_annotation = HeatmapAnnotation(Time_Post_Partum_Days = ann_df[colnames(mat.z), ]$time_post_partum_days,
                                                         HMO_concentration = ann_df[colnames(mat.z), ]$hmo_concentration,
                                                         col = list(Time_Post_Partum_Days = Var_tpp, HMO_concentration = Var_conc), 
                                                         show_annotation_name = FALSE),
                      right_annotation = ha,
                      row_names_gp = gpar(fontsize = 7, fontfamily = "sans", fontface = "italic"),
                      column_names_gp = grid::gpar(fontsize = 3))
        #pdf(file = paste0("Heatmap_", gsub(".csv", ".pdf", hmo_name)), height = 6, width = 8)
        file_path <- paste0("Heatmap_", hmo_name, ".pdf")
        pdf(file = file_path, height = 6, width = 8)
        print(ht)
        dev.off()
      }
  },
    error = function(err) {
      print("Error during hmo:")
      print(hmo_name)
      print(errorCondition(err))
    },
    warning = function(warn) {
      print("Warning during hmo:")
      print(hmo_name)
    },
    finally = {
    }
  )
}
