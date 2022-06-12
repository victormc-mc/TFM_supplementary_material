pipolin_cluster = {}
clust_list = []
clust_freq = {}

with open("graph_eval_E-5_mcl_clustered.csv", "r") as f:
    for line in f:
        if "__mclCluster" not in line:
            fields = line.replace('"','').split(',')

            if "G_" in fields[-1]:
                pipolin = "_".join(fields[-1].split("_")[:3])
            else:
                pipolin = "_".join(fields[-1].split("_")[:2])
            cluster = fields[0]
            
            singleton = False
            if cluster == '':
                singleton = True
            
            if not singleton:
                if pipolin not in pipolin_cluster.keys():
                    pipolin_cluster[pipolin] = [cluster]
                else:
                    pipolin_cluster[pipolin].append(cluster)
            
                if cluster not in list(clust_freq.keys()):
                    clust_list.append(cluster)
                    clust_freq[cluster] = 1
                else:
                    clust_freq[cluster] += 1


print(len(list(pipolin_cluster.keys())))
print(clust_freq)

matrix = "Corner\t"+"\t".join(clust_list)+"\n"

for pipolin_key in pipolin_cluster:
    matrix += pipolin_key+"\t"
    for cluster in clust_list:
        if cluster in pipolin_cluster[pipolin_key]:
            matrix += "1\t"
        else:
            matrix += "0\t"
    matrix += "\n"

with open("matrix.txt", "w") as f:
    f.write(matrix)


###5 or more matrix
matrix_15 = ""
clust_15 = []
for clu in clust_freq:
    if clust_freq[clu] >= 15:
        clust_15.append(clu)
        matrix_15 += "\t"+clu
matrix_15 = "Corner"+matrix_15+"\n"

for pipolin_key in pipolin_cluster:
    matrix_15 += pipolin_key+"\t"
    for cluster in clust_15:
        if cluster in pipolin_cluster[pipolin_key]:
            matrix_15 += "1\t"
        else:
            matrix_15 += "0\t"
    matrix_15 += "\n"

with open("matrix_15.txt", "w") as f:
    f.write(matrix_15)
