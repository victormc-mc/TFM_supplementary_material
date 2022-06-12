library("ggplot2")
library("RColorBrewer")
data <- read.table("matrix_15.txt",  header = TRUE)
data_num <- data
data_num$Corner <- NULL
matrix <- as.matrix(data_num)
rownames(matrix) <- data$Corner

col_order <- data.frame("Pipolin_ID" = data$Corner)
col_data <- read.table("genus_col.txt", header=TRUE)
join_df <- merge(col_order, col_data, sort=FALSE)

col_code <- as.numeric(as.factor(join_df$Genus))
ramp_man <- c("#ff0000", "#944d15", "#7d04d2", "#808080", "#ff9200", 
              "#EFF3FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", 
              "#08519C", "#ec7fb4", "#ff009e", "#6e3f0a", "#93c47d", 
              "#869120") #
colSide <- ramp_man[col_code]


svg("heatmap_clustering_10.svg", width = 15, height = 10)
heatmap(matrix, scale = "none",
        cexRow = 0.15,cexCol=1,
        RowSideColors = colSide,
        Rowv = NULL, Colv = NULL)


dev.off()

#BiocManager::install("ComplexHeatmap")
