filter_pipolins = ["Pipolin_ID", "G_49878_0", "G_49885_1", "G_100479_0", "G_129016_0", "G_131542_0", "G_131949_1", "G_131987_0", "G_136322_0"]

G_pred = {}
with open("pipolin_summary_updated.txt", "r") as f:
    for line in f:
        if line.split()[0] not in filter_pipolins:
            G_pred[line.split()[0]] = [line.split("\t")[17]]

#1 plasmid finder
plasmid_finder_R = {}
with open("plasmidfinder_pipolins_firmicutes_feb2022_filtered.tsv", "r") as f:
    for line in f:
        if "G_" in line:
            plasmid_finder_R[line.split("\t")[4].split()[0]] = "plasmid"

#2 MOB-suite
MOBSuite_R = {}
with open("MOBSuit_recon_firmicutes_feb2022.txt", "r") as f:
    for line in f:
        if "G_" in line:
            pipolin_id = line.split("\t")[4].split("_G_")[0]
            prediction = line.split("\t")[1]
            MOBSuite_R[pipolin_id] = prediction


#3 MOB-scan
MOBScan_R = {}
with open("MOBScan_web.txt", "r") as f:
    for line in f:
        if "G_" in line:
            pipolin_id = "_".join(line.split("_")[0:3])
            MOBScan_R[pipolin_id] = "plasmid"


#4 plasmidVerify
plasmidVerify_R = {}
with open("pipolins_firmicutes_feb2022_renamed_result_table.csv", "r") as f:
    for line in f:
        pipolin_id = line.split(",")[0]
        prediction = line.split(",")[1].lower()
        plasmidVerify_R[pipolin_id] = prediction

for G in G_pred:
    if G in list(plasmid_finder_R.keys()):
        G_pred[G].append("1")
    else:
        G_pred[G].append("0")
    if G in list(MOBSuite_R.keys()):
        if MOBSuite_R[G] == "plasmid":
            G_pred[G].append("1")
        else:
            G_pred[G].append("0")
    else:
        G_pred[G].append("0")
    if G in list(MOBScan_R.keys()):
        G_pred[G].append("1")
    else:
        G_pred[G].append("0")
    if G in list(plasmidVerify_R.keys()):
        if plasmidVerify_R[G] == "plasmid":
            G_pred[G].append("1")
        else:
            G_pred[G].append("0")
    else:
        G_pred[G].append("0")

matrix = "corner\tgenus\tplasmid_finder\tMOB-suite\tMOBScan\tplasmidVerify\n"
for G in G_pred:
    matrix += G+"\t"+"\t".join(G_pred[G])+"\n"

with open("consensus_matrix.txt", "w") as f:
    f.write(matrix)
