'''
Estado actual: parece que funciona. Mejorar report de tareas y stats pipolb
tal vez estaria bien guardar los pipolb positive genomes en tfp friendly format (para easy download)
probar ungroup en parallel

'''

'''
Part 1. Preparations must be made
'''
import multiprocessing
import os
import sys
import subprocess
import time
from datetime import datetime
from Bio import SeqIO

start_time = time.time()
current_day_h = str(datetime.now()).replace(" ","_").replace(":",".")[:-7]
resume = input("Are you resuming a previous screening? (Yes or No) ").upper()

if resume[0] == "Y":
    resume = True
    print("Resuming latest screening process.")
    rec_dir = input("Please, introduce the output folder of the screening to be resumed: ")
    os.chdir(rec_dir)
    for rec_folder in os.listdir():
        if "dataset" in rec_folder:
            fetch_folder_name = rec_folder
            fetch_path = fetch_folder_name + "/ncbi_dataset/"
    print(os.listdir())
    print("fetchpath: ", fetch_path)
     
elif resume[0] == "N":
    resume = False
    term = input("Enter taxon name or id: ").lower()
    filename = term.replace(" ","_")+".zip"
    out_prefix = input("Input prefix of results directory: ")
    os.mkdir("results_"+out_prefix+"_"+current_day_h)
    os.chdir("results_"+out_prefix+"_"+current_day_h)
    os.mkdir("tmp_genomes")
    os.mkdir("Results_EP")
    os.mkdir("EP_pool")
    os.mkdir("EP_error")
else:
    print("Your answer was taken as a no.")
    resume = False
    term = input("Enter taxon name or id: ").lower()
    filename = term.replace(" ","_")+".zip"
    out_prefix = input("Input prefix of results directory: ")
    os.mkdir("results_"+out_prefix+"_"+current_day_h)
    os.chdir("results_"+out_prefix+"_"+current_day_h)
    os.mkdir("tmp_genomes")
    os.mkdir("Results_EP")
    os.mkdir("EP_pool")
    os.mkdir("EP_error")


'''
Part 2. NCBI Datasets search
'''
if not resume:
    start_time_datasets = time.time()
    filter_answer = input("Do you want to rehydrating genomes from GenBank, RefSeq or Both (GenBank default)? ")
    subprocess.run(["datasets", "download", "genome", "taxon", term, "--dehydrated", "--filename", filename,
                    "--exclude-genomic-cds", "--exclude-gff3", "--exclude-protein", "--exclude-rna"])
    fetch_folder_name = filename[:-4]+"_dataset"
    subprocess.run(["unzip", filename, "-d", fetch_folder_name])
    subprocess.run(["rm", filename])

    fetch_path = fetch_folder_name+"/ncbi_dataset/"
    subprocess.run("cat "+fetch_path+"fetch.txt | grep -v json > "+fetch_path+"fetch_new.txt", shell=True)
    subprocess.run("rm "+fetch_path+"fetch.txt", shell = True)
    print("Splitting process:")
    subprocess.run("wc -l "+fetch_path+"fetch_new.txt", shell = True)
    #subprocess.run("split --verbose -d -a 10 -l 50 --additional-suffix '.txt' "+fetch_path+"fetch_new.txt "+fetch_path+"split_fetch_", shell = True) #just for linux (not mac) 
    
    if filter_answer == "Both":
        print("Rehydrating genomes from both GenBank and RefSeq.")
        filter_string = "/"
    elif filter_answer == "RefSeq":
        print("Rehydrating genomes from RefSeq.")
        filter_string = "/GCF_"
    else:
        print("Rehydrating genomes from GenBank.")
        filter_string = "/GCA_"

    with open(fetch_path+"fetch_new.txt", "r") as f:
        content = f.readlines()
    split_n = 0
    assembly_n = 0
    #Filter -> keep GCA (Genbank assemblies)
    Filtered_assemblies = []
    for GC_line in content:
        if filter_string in GC_line.split("\t")[2]:
            Filtered_assemblies.append(GC_line)

    ### Latest version filter ###
    cod_G = {}
    for line in Filtered_assemblies:
        assembly_version = line.split("\t")[-1].split("/")[1]
        assembly_AN = assembly_version.split(".")[0]
        if assembly_AN not in list(cod_G.keys()):
            cod_G[assembly_AN] = {assembly_version:[line]}
        else:
            if assembly_version not in list(cod_G[assembly_AN].keys()):
                cod_G[assembly_AN][assembly_version] = [line]
            else:
                cod_G[assembly_AN][assembly_version].append(line)
    
    cod_G_latest = {}
    for key in cod_G:
        latest_num = 0

        for version in cod_G[key]:
            num = eval(version.split(".")[1])
            if num > latest_num:
                latest_ver = cod_G[key][version]
                latest_AN = version
        cod_G_latest[latest_AN] = latest_ver
        
        if len(cod_G[key]) > 1:
            print("More than one version for", key, ":", list(cod_G[key].keys()))
            print(latest_AN, "selected")
       

    filtered_latest_assemblies = []
    for key in cod_G_latest:
        for ver in cod_G_latest[key]:
            filtered_latest_assemblies.append(ver)


    ### Incluir assembly split check
    split_n = 0
    while len(filtered_latest_assemblies)!=0:
        split_n += 1
        if len(filtered_latest_assemblies) > 100:
            cut_site = 100
            for i in range(100, len(filtered_latest_assemblies)):
                current_AN = filtered_latest_assemblies[i-1].split("\t")[-1].split("/")[1]
                next_AN = filtered_latest_assemblies[i].split("\t")[-1].split("/")[1]
                if current_AN != next_AN:
                    cut_site = i
                    break

        with open(fetch_path+"split_fetch_"+str(split_n)+".txt", "w") as f_out:
            f_out.write("".join(filtered_latest_assemblies[:cut_site]))
        filtered_latest_assemblies = filtered_latest_assemblies[cut_site:]

    print("Splitting completed. Length")
    subprocess.run("wc -l "+fetch_path+"split_fetch*", shell = True)
    start_split = 1
    print(time.time()-start_time_datasets, "seconds")
else:
    start_time_datasets = time.time()
    print("Removing temporal and incomplete files.")
    subprocess.run("rm -r tmp_genomes", shell = True)
    subprocess.run("mkdir tmp_genomes", shell = True)
    subprocess.run("rm -r tmp_results_EP", shell = True)
    subprocess.run("mkdir tmp_results_EP", shell = True)
    split_counter = 0
    for file in os.listdir(fetch_path):
        if "split" in file:
            split_counter += 1
        if ".recovery." in file:
            split_counter -= 1
            n_recovery = file.split(".")[0]
            n_recovery = n_recovery.split("_")[-1]
            RT_rec_file = file.replace(".recovery.txt", "")
            subprocess.run(["mv", fetch_path+file, fetch_path+RT_rec_file])
    split_n = split_counter+int(n_recovery)
    start_split = int(n_recovery)
    print("Recovery completed in: ", end="")
    print(time.time()-start_time_datasets, "seconds")
    print("Starting from split", start_split)

'''
Part 3. Genome processing
'''
n_tries = "5"
for i in range (start_split, split_n+2): #from split 1 to split + 2 more iterations to finish delayed tasks #Ej: for 5 splits: 1 (A) 2(AB) 3(AB) 4(AB) 5(AB) 6(B)
    print("\n###Iteration "+str(i)+".\n")
    start_genome_set=time.time()
    file = "split_fetch_"+str(i)+".txt"

    ### 1) Rehydartion
    if i == 1: 
        print("Starting rehydration process for split:", str(i), "out of", str(split_n))
        p_rehydration = subprocess.Popen('python3 ../genome_downloading.py '+str(i)+' '+fetch_path+" "+fetch_folder_name+" "+file+" "+n_tries+" results_"+current_day_h, shell=True)
        print("Rehydration process finished in:", str(time.time() - start_genome_set), "seconds.")
    elif 1 < i < split_n+1:
        print("Starting rehydration process for split:", str(i), "out of", str(split_n))
        if not resume:
            subprocess.run("mv tmp_genomes/* EP_pool/", shell=True) #mv genomes from earlier split
        else:
            print("Analyzing genomes from previous split")
            resume = False #From this point the screening will be as it was when it was stopped
        p_rehydration = subprocess.Popen('python3 ../genome_downloading.py '+str(i)+' '+fetch_path+" "+fetch_folder_name+" "+file+" "+n_tries+" results_"+current_day_h, shell=True) #, stdout=subprocess.DEVNULL
    elif i == split_n+1:
        print("Moving genomes from last split.")
        subprocess.run("mv tmp_genomes/* EP_pool/", shell=True) #mv genomes from last split

    ### 3) ExplorePipolin analysis (write better)
    if 1 < i < split_n+2:
        #Run EP
        print("Starting ExplorePipolin analysis for split: ", str(i-1))
        p_EP = subprocess.Popen('ls EP_pool/* | parallel --verbose "explore_pipolin --keep-tmp --out-dir tmp_results_EP --cpus 4 {}"', shell=True) #, stdout=subprocess.DEVNULL

    ### 4) Waiting time and check results
    if i > 1:
        # a) EP wait 
        p_EP.wait()
        print("ExplorePipolin process finished in:", str(time.time() - start_genome_set), "seconds.")
        print("Checking EP results:")
        #Check results (rm pipolb negative)
        for folder in os.listdir("tmp_results_EP"):
            EP_ERROR = False
            EP_piPolB = True
            with open("tmp_results_EP/"+folder+"/"+folder+".log") as f:
                for line in f:
                    if "ERROR" in line.upper():
                        EP_ERROR = True
                    if "No piPolBs were found!" in line:
                        EP_piPolB = False
            if EP_ERROR:
                subprocess.run("mv tmp_results_EP/"+folder+" EP_error/"+folder, shell=True)
                print("[E] EP raised error for:", folder)
            else:
                if EP_piPolB:
                    print("[+] Pipolins found in:", folder)
                    subprocess.run("mv tmp_results_EP/"+folder+" Results_EP/"+folder, shell=True)
                else:
                    print("[-] Pipolins not found in:", folder)
                    subprocess.run("rm -r tmp_results_EP/"+folder, shell=True)
            

    # b) Rehydration wait
    p_rehydration.wait()

    #Delete analised genomes
    if i > 1:
        subprocess.run("rm -r EP_pool/*", shell=True)

        
    ### 5) Info report
    print("Printing lists content")
    subprocess.run("wc -l *.txt", shell=True) #change by report.txt?
    subprocess.run("ls Results_EP/ | wc -l", shell = True)

    print("Processing this split took:", str((time.time() - start_genome_set)/60), "minutes.")
    print("Current running time is: ", str((time.time() - start_time)/60), "minutes.")
    print("\n###Iteration "+str(i)+" ended.\n\n")

'''
Part 4. Cleaning
'''
subprocess.run("rm -r tmp_genomes", shell = True)
subprocess.run("rm -r tmp_results_EP", shell = True)
subprocess.run("rm -r EP_pool", shell = True)

print("End of script")
print((time.time() - start_time)/60)
print("Minutes")





