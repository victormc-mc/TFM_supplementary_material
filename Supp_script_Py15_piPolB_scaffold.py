import os
from Bio import SeqIO
import subprocess

### piPolB len check

piPolB_len_dic = {}

for result in os.listdir("hmmsearch_pipolb_results_positive"):
    faa_name = result.replace("_piPolB_positive.hmm.tbl", ".fa")
    
    with open("all_proteins_fixed/"+faa_name, "r") as f:
        pipolb_len = 0
        for record in SeqIO.parse(f, "fasta"):
            pipolb_len += len(record.seq)
    
    #piPolB_name =  "_".join(result.split("_")[:3])+"_piPolB_"+result.split("_")[10]
    piPolB_len_dic[faa_name] = pipolb_len


### pipolin-piPolB genome dict

pipolin_piPolB_dict = {}

for piPolB in piPolB_len_dic:
    pipolin = "_".join(piPolB.split("_")[:3])
    if pipolin not in list(pipolin_piPolB_dict.keys()):
        pipolin_piPolB_dict[pipolin] = {piPolB:piPolB_len_dic[piPolB]}
    else:
        pipolin_piPolB_dict[pipolin][piPolB] = piPolB_len_dic[piPolB]

### Ubicate fragmetns with blast
if "firmicutes_partial_pipolbs" not in os.listdir():
    subprocess.run("mkdir firmicutes_partial_pipolbs", shell=True)
else:
    subprocess.run("rm -r firmicutes_partial_pipolbs", shell=True)
    subprocess.run("mkdir firmicutes_partial_pipolbs", shell=True)

if "firmicutes_blastRef_pipolbs" not in os.listdir():
    subprocess.run("mkdir firmicutes_blastRef_pipolbs", shell=True)
else:
    subprocess.run("rm -r firmicutes_blastRef_pipolbs", shell=True)
    subprocess.run("mkdir firmicutes_blastRef_pipolbs", shell=True)

if "firmicutes_scaffolded_pipolbs" not in os.listdir():
    subprocess.run("mkdir firmicutes_scaffolded_pipolbs", shell=True)
else:
    subprocess.run("rm -r firmicutes_scaffolded_pipolbs", shell=True)
    subprocess.run("mkdir firmicutes_scaffolded_pipolbs", shell=True)

for pipolin in pipolin_piPolB_dict:
    if len(pipolin_piPolB_dict[pipolin]) > 1:
        for piPolB in pipolin_piPolB_dict[pipolin]:
            subprocess.run("cp all_proteins_fixed/"+piPolB+" firmicutes_partial_pipolbs/"+piPolB,shell=True)
    else:
        piPolB = list(pipolin_piPolB_dict[pipolin].keys())[0] 
        subprocess.run("cp all_proteins_fixed/"+piPolB+" firmicutes_scaffolded_pipolbs/"+piPolB,shell=True)

subprocess.run("ls firmicutes_partial_pipolbs/* | parallel --verbose 'blastp -num_alignments 1 -outfmt 7 -max_hsps 1 -query {} -subject piPolB_ref.fasta -out firmicutes_blastRef_pipolbs/{/.}.blast.out'", shell=True)


pipolb_fragments_range_dict = {}

for blast_result in os.listdir("firmicutes_blastRef_pipolbs"):
    piPolB = blast_result.replace(".blast.out", ".fa")
    with open("firmicutes_blastRef_pipolbs/"+blast_result, "r") as f:
        for line in f:
            if "#" != line[0]:
                pipolb_fragments_range_dict[piPolB] = [eval(line.split()[8]), eval(line.split()[9])]
                


# ordering and scaffolding piPolBs
pipolin_scaffolded_piPolBs = {}
for pipolin in pipolin_piPolB_dict:
    if len(pipolin_piPolB_dict[pipolin]) > 1:

        piPolB_len_dic_subset = {}
        for piPolB in pipolin_piPolB_dict[pipolin]:
            piPolB_len_dic_subset[piPolB] = piPolB_len_dic[piPolB]

        initial_order = dict(sorted(piPolB_len_dic_subset.items(), key=lambda item: -item[1]))

        scaffolded_pipolbs = {}

        piPolB_range = "Empty"
        for piPolB_fragment in initial_order:
            piPolB_fragment_range = range(pipolb_fragments_range_dict[piPolB_fragment][0]+5, pipolb_fragments_range_dict[piPolB_fragment][-1]-5)
            if piPolB_range == "Empty":
                piPolB_range = [piPolB_fragment_range]
                scaffolded_pipolbs[piPolB_fragment] = [(piPolB_fragment_range[0]+piPolB_fragment_range[-1])/2, piPolB_fragment_range[0], piPolB_fragment_range[-1]]
            else:
                check = 0
                for subrange in piPolB_range:
                    if piPolB_fragment_range[0] in subrange or piPolB_fragment_range[-1] in subrange:
                        check += 1
                
                if check == 0:
                    scaffolded_pipolbs[piPolB_fragment] = [(piPolB_fragment_range[0]+piPolB_fragment_range[-1])/2, piPolB_fragment_range[0], piPolB_fragment_range[-1]]
                    piPolB_range.append(piPolB_fragment_range)
        
        pipolin_scaffolded_piPolBs[pipolin] = dict(sorted(scaffolded_pipolbs.items(), key=lambda item: item[1][0]))


###Write new piPolBs
for pipolin in pipolin_scaffolded_piPolBs:
    fragmented_piPolBs = list(pipolin_scaffolded_piPolBs[pipolin].keys())
    print(pipolin_scaffolded_piPolBs[pipolin])
    file_info = fragmented_piPolBs[0].replace(".fa", ".scaffolded.fa")
    file_name = "firmicutes_scaffolded_pipolbs/"+file_info

    piPolB_paths = []
    for fragment in fragmented_piPolBs:
        piPolB_paths.append("all_proteins_fixed/"+fragment)

    print(piPolB_paths)
    subprocess.run("head -n 1 all_proteins_fixed/"+fragmented_piPolBs[0]+" > "+file_name, shell=True)
    subprocess.run("cat "+" ".join(piPolB_paths)+ "| grep -v '>' | tr -d '*' >> "+file_name, shell=True)



### Filter piPolBs scaf

filter_list = ['G_131247', 'G_48777', 'G_131451', 'G_131572', 'G_123932', 'G_50173', 'G_48873', 'G_132944', 'G_124008', 'G_169491', 'G_124016', 'G_124119', #no RNA
               'G_149386', 'G_136197', 'G_132856', 'G_49870', 'G_136176'] #Contaminated

if "firmicutes_scaffolded_pipolbs_clean" not in os.listdir():
    subprocess.run("mkdir firmicutes_scaffolded_pipolbs_clean", shell=True)
else:
    subprocess.run("rm -r firmicutes_scaffolded_pipolbs_clean", shell=True)
    subprocess.run("mkdir firmicutes_scaffolded_pipolbs_clean", shell=True)

for piPolB_scaf in os.listdir("firmicutes_scaffolded_pipolbs"):
    genome = "G_"+piPolB_scaf.split("_")[1]
    if genome not in filter_list:
        print("cp firmicutes_scaffolded_pipolbs/"+piPolB_scaf+" firmicutes_scaffolded_pipolbs_clean/"+piPolB_scaf)
        subprocess.run("cp firmicutes_scaffolded_pipolbs/"+piPolB_scaf+" firmicutes_scaffolded_pipolbs_clean/"+piPolB_scaf, shell=True)






























### longest piPolB selection
"""
pipolin_best_piPolB = {}
for piPolB in piPolB_len_dic:
    pipolin = piPolB.split("_piPolB_")[0]
    if pipolin not in list(pipolin_best_piPolB.keys()):
        pipolin_best_piPolB[pipolin] = piPolB
    else:
        if piPolB_len_dic[piPolB] > piPolB_len_dic[pipolin_best_piPolB[pipolin]]:
            pipolin_best_piPolB[pipolin] = piPolB


if "firmicutes_longest_pipolbs" not in os.listdir():
    subprocess.run("mkdir firmicutes_longest_pipolbs", shell=True)
else:
    subprocess.run("rm -r firmicutes_longest_pipolbs", shell=True)
    subprocess.run("mkdir firmicutes_longest_pipolbs", shell=True)

for piPolB in pipolin_best_piPolB:
    pipolin = pipolin_best_piPolB[piPolB].split("_piPolB_")[0]
    num = "_"+pipolin_best_piPolB[piPolB].split("_piPolB_")[1]+".faa"
    for file in os.listdir("all_proteins_fixed"):
        if pipolin in file and num in file:
            subprocess.run("cp all_proteins_fixed/"+file+" firmicutes_longest_pipolbs/"+file, shell=True)
"""




"""
### NO FUNCIONA PORQUE NO TIENE EN CUENTA LA STRAND

### semi-scafolding (just adjacent piPolBs)
piPolB_len_dic_plus_adj = {}
top = 0
for piPolB in piPolB_len_dic:
    num = piPolB.split("_piPolB_")[1]
    prev = piPolB.split("_piPolB_")[0]+"_piPolB_"+str(eval(num)-1)
    next = piPolB.split("_piPolB_")[0]+"_piPolB_"+str(eval(num)+1)
    if prev in list(piPolB_len_dic_plus_adj.keys()):
        piPolB_len_dic_plus_adj[prev+"+"+piPolB] = piPolB_len_dic_plus_adj[prev] + piPolB_len_dic[piPolB]
        piPolB_len_dic_plus_adj.pop(prev)
    elif next in list(piPolB_len_dic_plus_adj.keys()):
        piPolB_len_dic_plus_adj[piPolB+"+"+next] = piPolB_len_dic_plus_adj[next] + piPolB_len_dic[piPolB]
        piPolB_len_dic_plus_adj.pop(next)
    else:
        piPolB_len_dic_plus_adj[piPolB] = piPolB_len_dic[piPolB]

### Select longest from the scaffolded
pipolin_best_adj_piPolB = {}
for piPolB in piPolB_len_dic_plus_adj:
    pipolin = piPolB.split("_piPolB_")[0]
    if pipolin not in list(pipolin_best_adj_piPolB.keys()):
        pipolin_best_adj_piPolB[pipolin] = piPolB
    else:
        if piPolB_len_dic_plus_adj[piPolB] > piPolB_len_dic_plus_adj[pipolin_best_adj_piPolB[pipolin]]:
            pipolin_best_adj_piPolB[pipolin] = piPolB

if "firmicutes_longest_adj_pipolbs" not in os.listdir():
    subprocess.run("mkdir firmicutes_longest_adj_pipolbs", shell=True)
else:
    subprocess.run("rm -r firmicutes_longest_adj_pipolbs", shell=True)
    subprocess.run("mkdir firmicutes_longest_adj_pipolbs", shell=True)

for pipolin in pipolin_best_adj_piPolB:
    if "+" in pipolin_best_adj_piPolB[pipolin]:
        fragments = pipolin_best_adj_piPolB[pipolin].split("+")
        piPolB_list = []
        num_list = []
        for fragment in fragments:
            num = fragment.split("_piPolB_")[1]
            num_list.append(num)
            piPolB_list.append("all_proteins_fixed/firmicutes_prodigal_partial.id_"+pipolin+"*"+"_"+num+".faa") #pipolb file

        file_name = "firmicutes_longest_adj_pipolbs/firmicutes_prodigal_partial.id_"+pipolin+"_"+num_list[0]+"_plus_adj.faa"
        subprocess.run("head -n 1 "+piPolB_list[0]+" > "+file_name, shell=True)
        subprocess.run("cat "+" ".join(piPolB_list)+ "| grep -v '>' | tr -d '*' >> "+file_name, shell=True)
    else:
        num = "_"+pipolin_best_adj_piPolB[pipolin].split("_piPolB_")[1]+".faa"
        for file in os.listdir("all_proteins_fixed"):
            if pipolin in file and num in file:
                subprocess.run("cp all_proteins_fixed/"+file+" firmicutes_longest_adj_pipolbs/"+file, shell=True)

"""