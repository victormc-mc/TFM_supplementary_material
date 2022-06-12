library('gtools')
library(RColorBrewer)
data <- read.table("relaxases_MOBP_matrix_id.txt", header = TRUE)
data <- data[mixedorder(data$corner),]
data_num <- data
data_num$corner <- NULL
matrix_rel <- as.matrix(data_num)
rownames(matrix_rel) <- data$corner
svg("heatmap_MOBP_relaxases_id.svg", width = 30, height = 30)
heatmap(matrix_rel, cexRow = 0.3,cexCol=0.3,
        #Rowv = NA, Colv = NA
        )
dev.off()


library(pheatmap)
pheatmap(matrix_rel, display_numbers = T, fontsize_number = 5,
         fontsize_row = 2, fontsize_col = 2)




library(reshape2)
library(ggplot2)
df_ggplot <- setNames(melt(matrix_rel), c('rows', 'vars', 'values'))

p <- ggplot(df_ggplot, aes(rows, vars)) +
  geom_tile(aes(fill = values)) + 
  #geom_text(aes(label = round(values, 1))) +
  scale_fill_gradient(low = "white", high = "red") 
p <- p + theme(text = element_text(size = 0.001)) 
p
