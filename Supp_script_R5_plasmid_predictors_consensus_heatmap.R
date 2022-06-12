library('gtools')
library(RColorBrewer)
data <- read.table("consensus_matrix.txt",  header = TRUE)
data <- data[mixedorder(data$corner),]
data_num <- data
data_num$corner <- NULL
data_num$genus <-NULL
matrix <- as.matrix(data_num)
rownames(matrix) <- data$corner



coul <- colorRampPalette(brewer.pal(8, "Oranges"))(25)
coul[1] <- "#FFFFFF"

col_order <- data.frame("Pipolin_ID" = data$corner)
col_data <- read.table("genus_col.txt", header=TRUE)
join_df <- merge(col_order, col_data, sort=FALSE)

col_code <- as.numeric(as.factor(join_df$Genus))
ramp_man <- c("#ff0000", "#944d15", "#7d04d2", "#808080", "#ff9200", 
              "#EFF3FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", 
              "#08519C", "#ec7fb4", "#ff009e", "#6e3f0a", "#93c47d", 
              "#869120") #
colSide <- ramp_man[col_code]


svg("heatmap_consensus.svg", width = 10, height = 10)
heatmap(matrix, scale = "none", col = coul,
        cexRow = 0.15,cexCol=0.15,
        RowSideColors = colSide, 
        Rowv = NA, Colv = NA,
        margins = c(0.1,20))
#legend(x="right", legend=c("min", "med", "max"),fill=heat.colors(3))
dev.off()
