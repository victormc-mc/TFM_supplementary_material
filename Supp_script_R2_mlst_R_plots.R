# library
library(ggplot2)

data_mlst <- read.table("mlst_firmicutes_feb_2022_updated.txt", header=TRUE, sep = "\t")
filtered_data_mlst <- data_mlst[data_mlst$ST!="NA",]
freq_mlst <- as.data.frame(table(filtered_data_mlst$Specie, filtered_data_mlst$ST))

ggplot(freq_mlst, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="stack", stat="identity") +
    coord_flip() +
    xlab("Genomes") + ylab("Species") + ggtitle("In silico MLST according to pubMLST schemes")


