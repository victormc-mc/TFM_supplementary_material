genus_name = []
G_info = {}


with open("genome_summary.txt") as f:
    for line in f:
        fields = line.split("\t")
        G_info[fields[0]] = [fields[3], fields[4]]

        if fields[3] not in genus_name:
            genus_name.append(fields[3])

total_genes = {}
total_unique_genes = {}
total_unique_genes_num = {}
total_subclass = {}
total_subclass_num = {}
total_class = {}
total_class_num = {}

uniq_G = []

with open("firmicutes_feb_2022_pipolin_amrfinder.txt") as f:
    for line in f:
        fields = line.split("\t")
        G = fields[0]

        if G not in uniq_G:
            uniq_G.append(G)

        gene = fields[6]
        AMR_class = fields[11].split("/")
        AMR_subclass = fields[12].split("/")
        element_type = fields[9]

        if element_type == "AMR":
            if G not in list(total_genes.keys()):
                total_genes[G] = 1
                total_unique_genes[G] = [gene]
                total_unique_genes_num[G] = 1
                total_subclass[G] = AMR_subclass
                total_subclass_num[G] = 1
                total_class[G] = AMR_class
                total_class_num[G] = 1

            else:
                total_genes[G] += 1
                if gene not in total_unique_genes[G]:
                    total_unique_genes[G].append(gene)
                    total_unique_genes_num[G] += 1
                for subclass in AMR_subclass:
                    if subclass not in total_subclass[G]:
                        total_subclass[G].append(subclass)
                        total_subclass_num[G] += 1
                for classs in AMR_class:
                    if classs not in total_class[G]:
                        total_class[G].append(classs)
                        total_class_num[G] += 1


full_table = ""
with open("genome_summary.txt") as f:
    for line in f:
        fields = line.split("\t")
        G = fields[0]

        if "G_" in G:
            if G in list(total_class.keys()):
                full_table += line.replace("\n", "\t")+";".join(total_class[G])+"\t"+str(total_class_num[G])+"\t"+";".join(total_subclass[G])+"\t"+str(total_subclass_num[G])+"\n"
            else:
                full_table += line.replace("\n", "\t")+"none\t0\tnone\t0\n"
        else:
            full_table += line.replace("\n", "\t")+"Total_class\tTotal_num\tTotal_subclass\tTotal_subclass_num\n"
with open("genome_summary_amr.txt", "w") as f:
    f.write(full_table)



gene_subclass_freq = {}
gene_class_freq = {}

for G in total_class:
    for AMR_class in total_class[G]:
        if AMR_class not in list(gene_class_freq.keys()):
            gene_class_freq[AMR_class] = 1
        else:
            gene_class_freq[AMR_class] += 1

for G in total_subclass:
    for AMR_subclass in total_subclass[G]:
        if AMR_subclass not in list(gene_subclass_freq.keys()):
            gene_subclass_freq[AMR_subclass] = 1
        else:
            gene_subclass_freq[AMR_subclass] += 1


classN_vs_genus = {}
for G in total_class_num:
    if G != "Genome_ID":
        N = total_class_num[G]
        genus = G_info[G][0]
        if genus not in list(classN_vs_genus.keys()):
            classN_vs_genus[genus] = {N:1}
        else:
            if N not in list(classN_vs_genus[genus].keys()):
                classN_vs_genus[genus][N] = 1
            else:
                classN_vs_genus[genus][N] += 1

for G in G_info:
    if G not in list(total_class_num.keys()) and "G_" in G:
        genus = G_info[G][0]
        if genus not in list(classN_vs_genus.keys()):
            classN_vs_genus[genus] = {0:1}
        else:
            if 0 not in list(classN_vs_genus[genus].keys()):
                classN_vs_genus[genus][0] = 1
            else:
                classN_vs_genus[genus][0] += 1

classN_vs_genus_table = "Genus\tAMR_class_Number\tFreq\n"
for genus in classN_vs_genus:
    for num in classN_vs_genus[genus]:
        classN_vs_genus_table += genus+"\t"+str(num)+"\t"+str(classN_vs_genus[genus][num])+"\n"

with open("classN_vs_genus_table.txt", "w") as f:
    f.write(classN_vs_genus_table)


gene_class_freq_table = "AMR_Class\tFreq\n"
for key in gene_class_freq:
    gene_class_freq_table += key+"\t"+str(gene_class_freq[key])+"\n"
with open("gene_class_freq_table.txt", "w") as f:
    f.write(gene_class_freq_table)



classID_vs_genus = {}
for G in total_class:
    if G != "Protein identifier":
        class_list = total_class[G]
        genus = G_info[G][0]
        for ID in class_list:
            if genus not in list(classID_vs_genus.keys()):
                classID_vs_genus[genus] = {ID:1}
            else:
                if ID not in list(classID_vs_genus[genus].keys()):
                    classID_vs_genus[genus][ID] = 1
                else:
                    classID_vs_genus[genus][ID] += 1
classID_vs_genus_table = "Genus\tAMR_class\tFreq\n"
for genus in classID_vs_genus:
    for num in classID_vs_genus[genus]:
        classID_vs_genus_table += genus+"\t"+str(num)+"\t"+str(classID_vs_genus[genus][num])+"\n"

with open("classID_vs_genus_table.txt", "w") as f:
    f.write(classID_vs_genus_table)



