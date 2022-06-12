import os
from subprocess import check_call
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Entrez

### 1) Read gbks from Results_EP/genome_id/pipolins/*single_record.gbk

pipolin_dictionary = {}
#form: Genome_id : [] 


for genome_id in os.listdir("Results_EP/"):
    for pipolin in os.listdir("Results_EP/"+genome_id+"/pipolins"):
        if "single_record" in pipolin and "v0" in pipolin:
            pipolin_id = pipolin.split(".")[0].replace("_v0", "")

            #Record genome_id [Field 1]
            pipolin_dictionary[pipolin_id] = [genome_id]

            #Record pipolin file [Field 2]
            pipolin_dictionary[pipolin_id].append(pipolin)

            #Record pipolin length [Field 3]
            rec = [re for re in SeqIO.parse("Results_EP/"+genome_id+"/pipolins/"+pipolin, "genbank")][0]
            pipolin_dictionary[pipolin_id].append(str(len(rec))) 

            #get pipolin identifier
            pipolin_fa = pipolin.replace("single_record.gbk", "fa")
            rec_fa = [re for re in SeqIO.parse("Results_EP/"+genome_id+"/pipolins/"+pipolin_fa, "fasta")][0]

            #Record contig length [Field 4], genome length [Field 5], contig_name [Field 6], contig_info [Field 7] and accesion_number [Field 8]
            full_genome_length = 0
            pipolin_contig = "NA"
            pipolin_contig_info = "NA"
            pipolin_contig_accesion_number = "NA"
            pipolin_contig_lenght = ""
            for fa in SeqIO.parse("Results_EP/"+genome_id+"/"+genome_id+".fa", "fasta"):
                full_genome_length += len(fa)
                if fa.id == rec_fa.id:
                    pipolin_contig = str(fa.id)
                    pipolin_contig_info = str(fa.description).split(" ")[2:]
                    pipolin_contig_accesion_number = str(fa.description).split(" ")[1]
                    pipolin_contig_lenght = str(len(fa))
            pipolin_dictionary[pipolin_id].append(pipolin_contig_lenght)
            pipolin_dictionary[pipolin_id].append(str(full_genome_length))
            pipolin_dictionary[pipolin_id].append(pipolin_contig)
            pipolin_dictionary[pipolin_id].append(" ".join(pipolin_contig_info))
            pipolin_dictionary[pipolin_id].append(pipolin_contig_accesion_number)

            #Record att presence [Field 9], number [Field 10], coordinates [Field 11] and length [Field 12]
            repeats = [feat for feat in rec.features if feat.type == "repeat_region"]
            if repeats == []:
                pipolin_dictionary[pipolin_id].append("False")
                pipolin_dictionary[pipolin_id].append("0")
                pipolin_dictionary[pipolin_id].append("NA")
                pipolin_dictionary[pipolin_id].append("NA")
                pipolin_dictionary[pipolin_id].append("NA")
            else: 
                pipolin_dictionary[pipolin_id].append("True")
                repeat_coord_list = []
                repeat_length_list = []
                functional_cord_list = []
                for repeat in repeats:
                    repeat_coord_list.append(str(repeat.location))
                    repeat_length_list.append(int(abs(repeat.location.end-repeat.location.start)))
                    functional_cord_list.append(range(int(repeat.location.start), int(repeat.location.end)))    
                pipolin_dictionary[pipolin_id].append(str(len(repeats)))
                pipolin_dictionary[pipolin_id].append(str(repeat_coord_list))
                pipolin_dictionary[pipolin_id].append(str(repeat_length_list))
            
                #Record integration site (features overlapping with repeat_regions)
                for feat in rec.features:
                    feat_overlap = False
                    for att_coord in functional_cord_list:
                        if feat.type != "repeat_region" and feat.type != "source":
                            if feat.location.start in att_coord or feat.location.end in att_coord:
                                feat_overlap = True
                    if feat_overlap:
                        if "product" in list(feat.qualifiers.keys()):
                            pipolin_dictionary[pipolin_id].append(str(feat.qualifiers["product"]))
                        else: #In case it has no product
                            pipolin_dictionary[pipolin_id].append(str(feat.type))

conversion_dict = {}
with open("conversion_table.txt", "r") as f:
    for line in f:
        conversion_dict[line.split("\t")[1].split(".")[0]] = line.strip("\n").split("\t")

#Record Genome AN [Field 13] and Genome File [Field 14]
for key in pipolin_dictionary:
    pipolin_dictionary[key].append(conversion_dict[pipolin_dictionary[key][0]][2])
    pipolin_dictionary[key].append(conversion_dict[pipolin_dictionary[key][0]][3])

n_count = 0
#Add organism, species name, assembly status
for key in pipolin_dictionary:
    n_count += 1
    Entrez.email = "vmc11298@gmail.com"
    AN = pipolin_dictionary[key][13]
    search_id = Entrez.esearch(db="assembly", term=AN+"[Assembly Accession]")
    
    check = False
    for i in range(5):
        try:
            content = Entrez.read(search_id)
            check = True
            break
        except:
            print(AN, "There was an error accesing Assembly databse. Retrying")
    
    if not check:
        print("Couldn't access, skipping.")
        continue

    id_num = content['IdList'][0]
    fetch_id = Entrez.efetch(db="assembly", id=id_num, rettype="docsum", retmode="xml")
    content2 = Entrez.read(fetch_id)
    #Record species name [Field 15], organism [Fiel 16], genus [Field 17], AssemblyStatus [Field 18]
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["Organism"])
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["Organism"].split()[0])
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["SpeciesName"])
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["AssemblyStatus"])
    print(key, AN)


    #BiosampleAN [Field 19], #BiosampleID [Field 20], #Biosource [Field 21], #Coverage [Field 22]
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["BioSampleAccn"])
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["BioSampleId"])
    pipolin_dictionary[key].append(str(content2['DocumentSummarySet']['DocumentSummary'][0]["Biosource"]))
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["Coverage"])
    

    #Host and isolation Source, SRA and geo location
    biosample_id_num = content2['DocumentSummarySet']['DocumentSummary'][0]["BioSampleId"]
    fetch_id_biosample = Entrez.efetch(db="biosample", id=biosample_id_num, rettype="full", retmode="text")
    tmp_dict = {"host":"NA", "iso_source":"NA", "SRA":"NA", "geo_loc":"NA"}
    for line in fetch_id_biosample:
        if "/host=" in line:
            tmp_dict["host"] = line.split('"')[1]
        if "/isolation source=" in line:
            tmp_dict["iso_source"] = line.split('"')[1]
        if "SRA:" in line:
            tmp_dict["SRA"] = line.split("SRA:")[1].replace("\n", "")
        if "/geographic location=" in line:
            tmp_dict["geo_loc"] = line.split('"')[1]
        if "Error" in line:
            print(line)
            print("Error in", biosample_id_num)
    
    for tmp_key in tmp_dict:
        pipolin_dictionary[key].append(tmp_dict[tmp_key])
    
    print(n_count)
    print(pipolin_dictionary[key])


#print or write dict #add haeder 
table_output = "Pipolin_ID\tGenome_ID\tPipolin_file\tPipolin_length\tContig_length\tGenome_length\tContig_name\tContig_description\tContig_AN\tAtts\tNumber_atts\tCoordinates_atts\tLength_atts\tIntegration_site\tGenome_AN\tGenome_file\tOrganism\tGenus\tSpecie\tAssemblyStatus\tBiosampleAN\tBiosampleID\tBiosource\tCoverage\tHost\tisolation_source\tSRA\tgeo_location\n"
for key in pipolin_dictionary:
    table_output += str(key)+"\t"+"\t".join(pipolin_dictionary[key])+"\n"

with open("pipolin_summary.txt", "w") as f:
    f.write(table_output)
