# library
library(ggplot2)

data <- read.table("classN_vs_genus_table.txt", header=TRUE)
ggplot(data, aes(fill=Genus, y=Freq, x=AMR_class_Number)) + 
    geom_bar(position="stack", stat="identity") + scale_x_continuous(breaks=0:12) +
    xlab("(Resistant to) Number of AMR classes") + ggtitle("Number of AMR class resistances predicted for each genome")
    

data <- read.table("classID_vs_genus_table.txt", header=TRUE, sep="\t")
ggplot(data, aes(fill=Genus, y=AMR_class, x=Freq)) + 
  geom_bar(position="stack", stat="identity") +
  xlab("Genomes") + ylab("AMR class") +
  ggtitle("Frequency of AMR class resistance by genus")

#Genomes with no AMR
sum(data[data$AMR_class_Number==0,]$Freq)
225-54

#MDR
sum(data[data$AMR_class_Number>2,]$Freq)
