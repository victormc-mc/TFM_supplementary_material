library(ggplot2)
library(tibble)
library(scales)
library(ggrepel)
library(forcats)
library(RColorBrewer)
library(randomcoloR)



genome_data <- read.table("genome_summary.txt", header = TRUE, sep="\t", quote = "\"")
pipolin_data <- read.delim("pipolin_summary_updated.txt", header = TRUE, sep="\t", quote = "$")
#assembly <- read.table("assembly_genebank_genus_freq.txt", header = TRUE)


#Pie chart
myPalette <- brewer.pal(12, "Paired") 
pie(table(genome_data$Genus), border="white", col=myPalette) 

#Bar chart
data_freq <- tibble(as.data.frame(table(genome_data$Genus)))

pipolins_genus_freq <- ggplot(data = data_freq, aes(x = reorder(Var1,Freq),Freq), y = Freq) +
  geom_col(aes(fill = Var1) , show.legend = FALSE) +
  ggtitle(paste("Genus frequency in host genomes")) +
  coord_flip() + 
  geom_label(aes(label = paste(Freq, percent(Freq/sum(Freq)), sep = "\n"),
                 y = -5, fill = Var1),
             show.legend = FALSE,
             size = 4, label.padding = unit(0.01, "lines")) 

pipolins_genus_freq

# Stacked
n <- 37
palette <- distinctColorPalette(n)
specie_genus_freq = as.data.frame(table(genome_data$Genus, genome_data$Specie))
ggplot(specie_genus_freq, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=palette, name = "Species") +
  coord_flip() +
  scale_x_discrete(limits = levels(reorder(specie_genus_freq$Var1,specie_genus_freq$Freq))) +
  xlab("Genus") + ylab("Genomes") + ggtitle("Species frequency in genomes with pipolins")


  
host_genus_freq = as.data.frame(table(genome_data$Genus, genome_data$Host))
ggplot(host_genus_freq, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  xlab("Genomes") + ylab("Host") + ggtitle("Genome host")


isor_genus_freq = as.data.frame(table(genome_data$Genus, genome_data$Isolation_source))
ggplot(isor_genus_freq, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  xlab("Genomes") + ylab("Isolation source") + ggtitle("Genome isolation source")


loc_genus_freq = as.data.frame(table(genome_data$Genus, genome_data$Geo_location))
ggplot(loc_genus_freq, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip()


#Assembly stats
assembly <- read.delim("assembly_stats.txt", header = TRUE, sep = "\t")
data_assembly <- tibble(assembly)

assembly_genus_freq <- ggplot(data = assembly, aes(x = reorder(Genus,Genomes),Genomes), y = Genomes) +
  geom_col(aes(fill = Genus) , show.legend = FALSE) +
  ggtitle(paste("Genus frequency in assembly databse")) +
  coord_flip() + 
  geom_label(aes(label = paste(Genomes, percent(Genomes/sum(Genomes)), sep = "\n"),
                 y = -2000, fill = Genus),
             show.legend = FALSE,
             size = 4, label.padding = unit(0.01, "lines")) 

assembly_genus_freq

### Pipolin and contig length histograms
ggplot(data=pipolin_data, aes(x=Pipolin_length, fill=Genus)) + geom_histogram(binwidth=1000) + 
  geom_vline(data=pipolin_data, aes(xintercept = 100000), colour="black", linetype="dotted") +
  xlab("Pipolin length (bp)") + ylab("Number of pipolins") + ggtitle("Pipolin length histogram")

ggplot(data=pipolin_data, aes(x=Contig_length, fill=Genus)) + geom_histogram(binwidth=10000) + 
  geom_vline(data=pipolin_data, aes(xintercept = 100000), colour="black", linetype="dotted") +
  xlab("Contig length (bp)") + ylab("Number of contigs") + ggtitle("Contig length histogram (all contigs)")

ggplot(data=pipolin_data, aes(x=Contig_length, fill=Genus)) + geom_histogram(binwidth=1000) + 
  geom_vline(data=pipolin_data, aes(xintercept = 100000), colour="black", linetype="dotted") +
  xlim(0, 100000) +
  xlab("Contig length (bp) Bacillus Salmonella Campylobacter") + ylab("Number of contigs") + ggtitle("Contig length histogram (contigs <100 kbp)")
