import os
import subprocess
from Bio import SeqIO

if "Results_EP_filtered" not in os.listdir():
    subprocess.run("mkdir Results_EP_filtered", shell = True)
    subprocess.run("mkdir Results_EP_filtered/pipolins_gbk_single_record", shell = True)
    subprocess.run("mkdir Results_EP_filtered/pipolins_fa", shell = True)
    subprocess.run("mkdir Results_EP_filtered/genomes", shell = True)
else:
    subprocess.run("rm -r Results_EP_filtered", shell = True)
    subprocess.run("mkdir Results_EP_filtered", shell = True)
    subprocess.run("mkdir Results_EP_filtered/pipolins_gbk_single_record", shell = True)
    subprocess.run("mkdir Results_EP_filtered/pipolins_fa", shell = True)
    subprocess.run("mkdir Results_EP_filtered/genomes", shell = True)

if "Results_EP_discarded" not in os.listdir():
    subprocess.run("mkdir Results_EP_discarded", shell = True)
    subprocess.run("mkdir Results_EP_discarded/pipolins_gbk_single_record", shell = True)
    subprocess.run("mkdir Results_EP_discarded/pipolins_fa", shell = True)
    subprocess.run("mkdir Results_EP_discarded/genomes", shell = True)
else:
    subprocess.run("rm -r Results_EP_discarded", shell = True)
    subprocess.run("mkdir Results_EP_discarded", shell = True)
    subprocess.run("mkdir Results_EP_discarded/pipolins_gbk_single_record", shell = True)
    subprocess.run("mkdir Results_EP_discarded/pipolins_fa", shell = True)
    subprocess.run("mkdir Results_EP_discarded/genomes", shell = True)

len_list = []
piPolB_num = 0
piPolB_cds = 0
pipolin_num = 0
n = 0
for G in os.listdir("Results_EP"):
    n += 1
    if "000" in str(n):
        print(n)
    for pipolin in os.listdir("Results_EP/"+G+"/pipolins"):
        if "v0" in pipolin and "single_record" in pipolin:
            pipolin_num += 1
            piPolB_check = False
            pipolb_len_sum = 0
            cds_number = 0
            len_record = 0
            with open("Results_EP/"+G+"/pipolins/"+pipolin, "r") as f:
                for record in SeqIO.parse(f, "genbank"):
                    len_record += len(record.seq)
                    for feature in record.features:
                        if feature.type == "CDS":
                            cds_number += 1
                            for prod in feature.qualifiers["product"]:
                                if "Primer-independent DNA polymerase PolB" in prod:
                                    piPolB_check = True
                                    piPolB_cds += 1
                                    pipolb_len_sum += len(feature.qualifiers["translation"][0])
            len_list.append(pipolb_len_sum)

            pipolin_len_check = True
            if len_record > 1000000:
                pipolin_len_check = False

            pipolb_len_check = True
            if pipolb_len_sum < 400:
                pipolb_len_check = False
            
            cds_check = True
            if cds_number < 5:
                cds_check = False

            if pipolin_len_check and pipolb_len_check and cds_check:
                subprocess.run("cp Results_EP/"+G+"/pipolins/"+pipolin+" Results_EP_filtered/pipolins_gbk_single_record/"+pipolin, shell=True)
                subprocess.run("cp Results_EP/"+G+"/pipolins/"+pipolin.replace(".single_record.gbk", ".fa") +" Results_EP_filtered/pipolins_fa/"+pipolin.replace(".single_record.gbk", ".fa"), shell=True)
                subprocess.run("cp Results_EP/"+G+"/"+G+".fa"+" Results_EP_filtered/genomes/"+G+".fa", shell=True)
            else:
                subprocess.run("cp Results_EP/"+G+"/pipolins/"+pipolin+" Results_EP_discarded/pipolins_gbk_single_record/"+pipolin, shell=True)
                subprocess.run("cp Results_EP/"+G+"/pipolins/"+pipolin.replace(".single_record.gbk", ".fa") +" Results_EP_discarded/pipolins_fa/"+pipolin.replace(".single_record.gbk", ".fa"), shell=True)
                subprocess.run("cp Results_EP/"+G+"/"+G+".fa"+" Results_EP_discarded/genomes/"+G+".fa", shell=True)

            if piPolB_check:
                piPolB_num += 1
            else:
                print(pipolin, "doesn't meet the criteria")

print("Expected piPolBs:", pipolin_num)
print("Pipolins with piPolBs:", piPolB_num)
print("Total piPolBs CDSs:", piPolB_cds)
