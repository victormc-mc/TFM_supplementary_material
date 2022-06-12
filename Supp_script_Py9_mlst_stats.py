from dataclasses import field

pubmlst_species  = ["Staphylococcus epidermidis", "Staphylococcus aureus", "Staphylococcus haemolyticus", "Staphylococcus hominis", "Ligilactobacillus salivarius"]
contaminated_genomes = ["G_132856", "G_136197", "G_149386", "G_136176", "G_49870"]
G_info = {}
with open("genome_summary.txt") as f:
    for line in f:
        fields = line.split("\t")
        G_info[fields[0]] = fields[4]
ST_info = {}
with open("mlst_firmicutes_feb_2022.txt") as f:
    for line in f:
        fields = line.split("\t")
        G = fields[0].split("/")[-1].replace(".fa","")
        if G not in contaminated_genomes and G_info[G] in pubmlst_species:
            ST_info[G] = G_info[G]+"\t"+"ST"+fields[2]
        else:
            ST_info[G] = G_info[G]+"\t""NA"

mlst_updated = "Genome\tSpecie\tST\n"
for key in ST_info:
    mlst_updated += key +"\t"+ST_info[key]+"\n"

with open("mlst_firmicutes_feb_2022_updated.txt", "w") as f:
    f.write(mlst_updated)


"""
mlst_spec_freq = {}
with open("mlst_firmicutes_feb_2022.txt.txt") as f:
    for line in f:
        fields = line.split("\t")
        G = fields[0].split("/")[1]
        specie = fields[1]
        ST = fields[2]
        if G not in contaminated_genomes and len(fields)>4:
            if specie not in list(mlst_spec_freq.keys()):
                mlst_spec_freq[specie] = {ST:1}
            else:
                if ST not in list(mlst_spec_freq[specie].keys()):
                    mlst_spec_freq[specie][ST] = 1
                else:
                    mlst_spec_freq[specie][ST] += 1

out_table = "Specie\tST\tFreq\n"
for specie in mlst_spec_freq:
    for ST in mlst_spec_freq[specie]:
        out_table += specie+"\t"+ST+"\t"+str(mlst_spec_freq[specie][ST])+"\n"

with open("mlst_freq_table.txt", "w") as f:
    f.write(out_table)


"""



"""
        if len(fields) > 4:
            print(fields)

genus_name = []
G_info = {}

with open("genome_summary.txt") as f:
    for line in f:
        fields = line.split("\t")
        G_info[fields[0]] = fields[3]

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

with open("firmicutes_genomes_amrfinder.txt") as f:
    for line in f:
        fields = line.split("\t")
        G = fields[0]

        if G not in uniq_G:
            uniq_G.append(G)

        gene = fields[2]
        AMR_class = fields[7]
        AMR_subclass = fields[8]

        if G not in list(total_genes.keys()):
            total_genes[G] = 1
            total_unique_genes[G] = [gene]
            total_unique_genes_num[G] = 1
            total_subclass[G] = [AMR_subclass]
            total_subclass_num[G] = 1
            total_class[G] = [AMR_subclass]
            total_class_num[G] = 1

        else:
            total_genes[G] += 1
            if gene not in total_unique_genes[G]:
                total_unique_genes[G].append(gene)
                total_unique_genes_num[G] += 1
            if AMR_subclass not in total_subclass[G]:
                total_subclass[G].append(AMR_subclass)
                total_subclass_num[G] += 1
            if AMR_class not in total_class[G]:
                total_class[G].append(AMR_subclass)
                total_class_num[G] += 1


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
    if G != "Protein identifier":
        N = total_class_num[G]
        genus = G_info[G]
        if genus not in list(classN_vs_genus.keys()):
            classN_vs_genus[genus] = {N:1}
        else:
            if N not in list(classN_vs_genus[genus].keys()):
                classN_vs_genus[genus][N] = 1
            else:
                classN_vs_genus[genus][N] += 1

for G in G_info:
    if G not in list(total_class_num.keys()) and "G_" in G:
        genus = G_info[G]
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
        genus = G_info[G]
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


#Pipolins statistics:
with open("firmicutes_pipolin_amrfinder.txt") as f:
    for line in f:
        print(line)
"""
