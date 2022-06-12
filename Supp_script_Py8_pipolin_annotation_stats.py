import os
from Bio import SeqIO

n_pipolin = 0
n_cds = 0
n_hypo = 0


for G in os.listdir("Results_EP"):
    for file in os.listdir("Results_EP/"+G+"/pipolins"):
        if "single_record" in file and "v0" in file:
            n_pipolin += 1
            with open("Results_EP/"+G+"/pipolins/"+file) as f:
                for record in SeqIO.parse(f, "genbank"):
                    for feat in record.features:
                        if feat.type == "CDS":
                            n_cds += 1
                            if feat.qualifiers["product"][0] == "hypothetical protein":
                                n_hypo += 1

print(n_pipolin)
print(n_cds)
print(n_hypo)