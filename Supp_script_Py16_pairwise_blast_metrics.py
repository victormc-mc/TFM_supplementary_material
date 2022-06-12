import os
import subprocess
from Bio import SeqIO
import math


subprocess.run("blastp -query T4SS_MOBP1_Firmicutes_relaxases_clean.faa  -subject T4SS_MOBP1_Firmicutes_relaxases_clean.faa -out T4SS_MOBP1_Firmicutes_relaxases_clean_pairwiseBLAST_rel.txt -outfmt 7 -num_threads 30 -max_target_seqs 10000 -evalue 0.1", shell = True)

max_eval = 0

blast_qs_stats = {}
with open("T4SS_MOBP1_Firmicutes_relaxases_clean_pairwiseBLAST_rel.txt", "r") as f:
    for line in f:
        if "#" not in line:
            fields = line.split("\t")
            query = fields[0]
            subject = fields[1]
            e_val = eval(fields[10])
            if e_val != 0.0:
                e_val = round(-math.log(e_val), 3)
            else:
                e_val = 500

            if e_val > max_eval:
                max_eval = e_val

            e_val = str(round(math.log(e_val), 3)) #doble transformacion para q se vea algo en el heatmap
            id_per = fields[2]
            score = fields[11].replace("\n", "")
            if query not in list(blast_qs_stats.keys()):
                blast_qs_stats[query] = {subject:[[e_val, id_per, score]]}
            else:
                if subject not in list(blast_qs_stats[query].keys()):
                    blast_qs_stats[query][subject] = [[e_val, id_per, score]]
                else:
                    blast_qs_stats[query][subject].append([e_val, id_per, score])

### CON lim e-value < 0.1 no hay mas de 1 hit per subject. Si quitas el parametro sí
for q in blast_qs_stats:
    for s in blast_qs_stats[q]:
        if len(blast_qs_stats[q][s]) > 1:
            print(q, s, blast_qs_stats[q][s])

e_val_heatmap = "corner\t"+"\t".join(list(blast_qs_stats.keys()))+"\n"
for q in blast_qs_stats:
    e_val_heatmap += q+"\t"
    for s in blast_qs_stats:
        if s in list(blast_qs_stats[q].keys()):
            e_val_heatmap += blast_qs_stats[q][s][0][0]+"\t"
        else:
            e_val_heatmap += "0\t"

    e_val_heatmap += "\n"

with open("relaxases_MOBP_matrix_e_val.txt", "w") as f:
    f.write(e_val_heatmap)



id_heatmap = "corner\t"+"\t".join(list(blast_qs_stats.keys()))+"\n"
for q in blast_qs_stats:
    id_heatmap += q+"\t"
    for s in blast_qs_stats:
        if s in list(blast_qs_stats[q].keys()):
            id_heatmap += blast_qs_stats[q][s][0][1]+"\t"
        else:
            id_heatmap += "0\t"

    id_heatmap += "\n"

with open("relaxases_MOBP_matrix_id.txt", "w") as f:
    f.write(id_heatmap)