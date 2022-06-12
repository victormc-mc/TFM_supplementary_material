pipolin_genus = ["Blautia","Butyricicoccus","Clostridium","Contaminated","Hespellia","Lacticaseibacillus",
                "Lactobacillus","Lentilactobacillus","Ligilactobacillus","Limosilactobacillus","Loigolactobacillus"
                "Porcincola","Roseburia","Ruminococcus","Staphylococcus","Trichococcus"]
genus_stats = {}
with open("assembly_summary_genbank.txt", "r") as f:
    for line in f:
        if "#" != line[0]:
            fields = line.split("\t")
            genus = fields[7].split()[0]
            if genus in pipolin_genus:
                if genus not in list(genus_stats.keys()):
                    genus_stats[genus] = 1
                else:
                    genus_stats[genus] += 1


output_table = "Genus\tGenomes\n"
for key in list(genus_stats.keys()):
    output_table += key+"\t"+str(genus_stats[key])+"\n"

with open("assembly_stats.txt", "w") as f:
    f.write(output_table)