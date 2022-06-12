world_count = {}
genomes_dict = {}
with open("pipolin_summary_updated.txt", "r") as f:
    for line in f:
        fields = line.replace("\n", "").split("\t")
        if "G_" in line:
            G = fields[1]
            if G not in list(genomes_dict.keys()):
                genomes_dict[G] = [1, fields[14], fields[17], fields[18], fields[20], fields[24], fields[25], fields[26], fields[27].split(":")[0]]
            else:
                genomes_dict[G][0] += 1
        
            country = fields[27].split(":")[0]
            if country in list(world_count.keys()):
                world_count[country] += 1
            else:
                world_count[country] = 1

table = "Genome\tPipolins\tAssemblyAN\tGenus\tSpecie\tBioSampleAN\tHost\tIsolation_source\tSRA\tGeo_location\n"
for key in genomes_dict:
    table += key + "\t"+str(genomes_dict[key][0])+ "\t" + "\t".join(genomes_dict[key][1:]) + "\n"

with open("genome_summary.txt", "w") as f:
    f.write(table)


print(world_count)
