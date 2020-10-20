import pandas as pd
import aa_matching as aa_mm
import re
from collections import defaultdict
import json
import requests

# generate list of IMGT/HLA alleles that correspond to DPB1 HVRs

DPB1_HVR_AlleleList = defaultdict(list)

# get hypervariable region array
url = "http://rest.hlatools.org/hladpb1/resources/hypervariableRegions"
response = requests.get(url)
# print(response.text)

hypervariableRegions = json.loads(response.text)

# pretty print JSON response with indenting
print(json.dumps(hypervariableRegions, indent = 4, sort_keys=True))

hvr_index = 0
hvr_name_list = []
hvr_codon_list = []
hvr_poly_to_name = {}  # 84G_85G_86P_87M or 84V_85G_86P_87M to GGPM|VGPM 
for hvr in hypervariableRegions:
    # print (hvr)
    hvr_name = hypervariableRegions[hvr_index]['hypervariableRegionName']
    codon_list = hypervariableRegions[hvr_index]['codonNumberList']
    variant_map = hypervariableRegions[hvr_index]['variantMap']

    hvr_index = hvr_index + 1
    print ("HVR_Name: " + hvr_name.capitalize())
    print ("Codon_List: " + str(codon_list))
    # print (variant_map)

    hvr_name_list.append(hvr_name.capitalize())
    hvr_codon_list.append(codon_list)

    # print each variant ID

    # for variant in variant_map:
    #     print (variant_map)
        # print ('Variant ID = ' + key, ':', variant_map[key])

    # print protein sequence for each variant ID

    for protein_sequence in variant_map.values():
        print ("Protein Sequence: " + str(protein_sequence['proteinSequenceList']))
        
        # send protein strings onto an array
        hvr_motif_list = protein_sequence['proteinSequenceList']

        hvr_name_motif = "|".join(hvr_motif_list)

        # hvr_motif_array = str(hvr_motif_array).replace("['","").replace("']","")
        hvr_name_string = 'DPB1_HVR_' + hvr_name.capitalize() + '_' + hvr_name_motif

        hvr_poly_list = []
        for motif in hvr_motif_list:
            print ("HVR Motif: " + motif)
            # for i in range(0,len(motif)):
            i = 0
            aa_poly_list = []
            for aa in motif:
                aa_pos = codon_list[i]
                i = i + 1
                print ("AA: " + str(aa_pos) + aa)
                aa_poly_list.append(str(aa_pos) + aa)
            hvr_poly_string = 'DPB1_HVR_' + hvr_name.capitalize() + "_" + '_'.join(aa_poly_list) 
            hvr_poly_to_name[hvr_poly_string] = hvr_name_string
         # print abbreviated name (DPB1_HVR_hypervariable region name_proteinsequence)


# print alleles and HVRs
HLA_DPB1_HVR_AlleleList = defaultdict(list)

for allele_loctyp in aa_mm.HLA_seq:
    (allele_loc, allele_typ) = allele_loctyp.split('*')

    if (allele_loc != "DPB1"):
        continue

    # position_list_int = ([8,9,11],[35,36],[55,56,57],[65,69],[76],[84,85,86,87],[96])
    for hvr_name, position_list in zip(hvr_name_list, hvr_codon_list):
        # print (hvr_name, position_list_int)

        hvr_motif = aa_mm.getEpitope(allele_loctyp,position_list)

        DPB1_HVR_poly = 'DPB1_HVR_' + hvr_name + '_' + hvr_motif

        if DPB1_HVR_poly in hvr_poly_to_name: # some polymorphisms don't have variant IDs
            DPB1_HVR_name = hvr_poly_to_name[DPB1_HVR_poly]
            HLA_DPB1_HVR_AlleleList[DPB1_HVR_name].append(allele_loctyp)
        # print (allele_loctyp +  " : DPB1_HVR_" + (hvr_name) + "_" + allele_epitope_string)
        
        print ("Allele " + allele_loctyp + " has HVR " + hvr_name + " with polymorphism(s) " + hvr_motif)

        HLA_DPB1_HVR_AlleleList[DPB1_HVR_poly].append(allele_loctyp)

# output allele list to file
DPB1_HVR_allelelist_filename = "AA_DPB1_HVR_to_alleles.csv"
DPB1_HVR_allelelist_file = open(DPB1_HVR_allelelist_filename, 'w')

DPB1_HVR_allelelist_file.write("AA_Poly,Alleles\n")

for AA_poly in HLA_DPB1_HVR_AlleleList:
    allelelist = HLA_DPB1_HVR_AlleleList[AA_poly]
    allele_string = ','.join(allelelist)
    DPB1_HVR_allelelist_file.write(AA_poly + "," + allele_string + "\n")

DPB1_HVR_allelelist_file.close()
# print (HLA_AA_AlleleList[])