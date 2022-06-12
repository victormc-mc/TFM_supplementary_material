from hashlib import new


contaminated_genomes = ["G_132856", "G_136197", "G_149386", "G_136176", "G_49870"]

new_summary = ""
with open("pipolin_summary.txt", "r") as f:
    for line in f:
        fields = line.split("\t")
        if "G_" in line:
            G = fields[1]
            if G in contaminated_genomes:
                updated_line = "\t".join(fields[:16])+"\tContaminated\tContaminated\tContaminated\t"+"\t".join(fields[19:])
                new_summary += updated_line.replace("\n", "\tContaminated\n")
            else:
                new_summary += line.replace("\n", "\tNot_contaminated\n")
        else:
            new_summary += line.replace("\n", "\tContamination_status\n")

with open("pipolin_summary_updated.txt", "w") as f:
    f.write(new_summary)