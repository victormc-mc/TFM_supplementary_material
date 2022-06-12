import os
from re import sub
import subprocess
import math


if "tmp_blast_results" not in os.listdir():
    subprocess.run("mkdir tmp_blast_results", shell=True)
else:
    subprocess.run("rm -r tmp_blast_results/", shell=True)
    subprocess.run("mkdir tmp_blast_results", shell=True)

if "tmp_blast_results" not in os.listdir():
    subprocess.run("mkdir tmp_blast_db", shell=True)
else:
    subprocess.run("rm -r tmp_blast_db/", shell=True)
    subprocess.run("mkdir tmp_blast_db", shell=True)

subprocess.run("makeblastdb -in firmicutes_prodigal_partial.faa -dbtype prot -out tmp_blast_db/db", shell=True) #proteins_subset.faa

graph = "target\tsource\te-value\t-log10_e-value_200min\t-log10_e-value_50min\n" 

subprocess.run("ls all_proteins/* | parallel --verbose 'blastp -query {} -db tmp_blast_db/db -outfmt 7 -num_threads 30 -max_target_seqs 10000 -out tmp_blast_results/{/.}_blsatp.txt'", shell = True)
n_q = 0
for result in os.listdir("tmp_blast_results"):
    n_q += 1
    print("Processing output : ", result, "number", n_q)
    with open("tmp_blast_results/"+result, "r") as f:
        for line in f:
            if "#" != line[0]:
                fields = line.split()
                if fields[0] == fields[1]:
                    graph += fields[0]+"\n"    
                else:
                    if eval(fields[10]) <= 1E-5:
                        if eval(fields[10]) == 0:
                            log_eval = 200
                            log_eval50 = 50
                        elif eval(fields[10]) <= 1E-50:
                            log_eval = -math.log(eval(fields[10]),10)
                            log_eval50 = 50
                        else:
                            log_eval = -math.log(eval(fields[10]),10)
                            log_eval50 = -math.log(eval(fields[10]),10)
                        graph += fields[0]+"\t"+fields[1]+"\t"+fields[10]+"\t"+str(log_eval)+"\t"+str(log_eval50)+"\n"

with open("graph_eval_E-5.txt", "w") as f:
    f.write(graph)



"""
file_n = 0
for file in os.listdir("tmp_blast_results"):
    file_n += 1
    if "000" in str(file_n)[-3:]:
        print(file_n)
    with open("tmp_blast_results/"+file, "r") as f:
        n = 0
        for line in f:
            if "#" not in line:
                n += 1
                if n == 1:
                    query = line.split()[0]
                    subject = line.split()[1]
                    e_value = line.split()[10]
                    if eval(e_value) < 1E-10:
                        if e_value == "0.0":
                            e_value_log = -250
                        else:
                            e_value_log = math.log(eval(e_value),10)
                        graph += query+"\t"+subject+"\t"+str(e_value_log)+"\n"

"""
