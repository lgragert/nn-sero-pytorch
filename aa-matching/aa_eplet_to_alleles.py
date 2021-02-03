import pandas as pd
import aa_matching as aa_mm
import re
from collections import defaultdict


def complex_poly_to_eplet_poly(complex_poly):

    poly_list = complex_poly.split('|')
    eplet_poly_list = []

    for poly in poly_list:
        aa_array = re.split('([A-Z])',poly)

        eplet_list = []
    
        # parsing the polymorphic AA list
        for j in range (0,len(aa_array)):
            if (aa_array[j] == ""): # skip blank
                continue

            if ((j % 2) == 1): # skip every other array index
                continue

            eplet = aa_array[j] + aa_array[j+1]
            eplet_list.append(eplet)

            j = j + 2
        
        eplet_list_string = '_'.join(eplet_list)
        eplet_poly_list.append(eplet_list_string)
    
    eplet_poly_string = "|".join(eplet_poly_list)

    return eplet_poly_string

# generate list of IMGT/HLA alleles that correspond to eplets in the Eplet Registry

HLA_Eplet_AlleleList = defaultdict(list)

# for loc in ["DRB1"]:
# for loc in ["DRB1","DRB3","DRB4","DRB5"]:
# for loc in ["A","B","C"]:
# for loc in ["DQA1","DQB1"]:
for loc in ["A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1"]:
# for loc in aa_mm.ard_start_pos:

    print (loc)
    # print (hlaProteinOffset[hla])

    # extract data from file Epitope Registry files

    if loc in ['A','B','C']:
        epitope_filename = "Epitope_Registry_ABC.txt"
    elif loc in ['DRB1','DRB3','DRB4','DRB5']:
        epitope_filename = "Epitope_Registry_DRB.txt"
    elif loc in ["DQA1","DQB1"]:
        epitope_filename = "Epitope_Registry_DQ.txt"
    elif loc in ["DPA1","DPB1"]:
        epitope_filename = "Epitope_Registry_DP.txt"
    else:
        print ("Invalid Locus: " + loc)

    epitope = pd.read_csv(epitope_filename,sep='\t')

    # print (epitope)

    (epitope_db_prefix,epitope_db_suffix) = epitope_filename.split('.')
    (ep,reg,epitope_db_name) = epitope_db_prefix.split('_') 

    for allele_loctyp in aa_mm.HLA_seq:
        (allele_loc, allele_typ) = allele_loctyp.split('*')

        if (allele_loc != loc):
            continue

        # print (loc + '*' + allele_typ)

        for i in range (0,len(epitope)):  # eplet loop 
 
            complex_poly = epitope['Polymorphic'][i]

            # eplet name cleanup
            complex_poly = complex_poly.replace(" ","") # remove spaces to concatenate polymorphisms
            complex_poly = complex_poly.replace("pairwith","") # 'pair with' means the same thing as spaces
            complex_poly = complex_poly.replace("163L-167G/S","163L167G|163L167S") # logical OR
            complex_poly = complex_poly.replace("(","") # remove opened parentheses
            complex_poly = complex_poly.replace(")","") # remove closed parantheses 
            if ('or' in complex_poly):  # remove eplet definitions for DP from DRB
                (complex_poly,post_or) = complex_poly.split('or')

            poly_list = complex_poly.split('|')

            # print ("Poly List: " + str(poly_list))

            eplet_name = "Eplet_" + epitope_db_name + '_' + epitope['Eplet'][i]
            eplet_poly_string = complex_poly_to_eplet_poly(complex_poly)

            # loop through all AA motifs
            for poly in poly_list:

                aa_array = re.split('([A-Z])',poly)

                eplet_list = []
                
                # parsing the polymorphic AA list
                for j in range (0,len(aa_array)):
                    if (aa_array[j] == ""): # skip blank
                        continue

                    if ((j % 2) == 1): # skip every other array index
                        continue

                    eplet = aa_array[j] + aa_array[j+1]
                    eplet_list.append(eplet)

                    j = j + 2

                # parsing eplet_list to create position_list

                position_string = re.sub('([A-Z])',' ', poly)
                position_string = position_string.replace(" ",",")
                position_string = position_string[:-1]
                
                position_list = position_string.split(",")

                # print ("Poly: " + poly)
                # print ("AA_array: " + str(aa_array))
                # print ("Eplet_list: " + str(eplet_list))
                # print ("Position_String: " + position_string)
                # print ("Position_List: " + str(position_list))

                position_list_int = [int(i) for i in position_list]

                eplet_list_string = '_'.join(eplet_list)
                # eplet_poly = "Eplet_" + epitope_db_name + '_' + '_'.join(eplet_list)
                eplet_poly = "Eplet_" + epitope_db_name + '_' + eplet_poly_string
            
                # print ("Position_List: " + str(position_list_int))
                # print ("Eplet_Name: " + eplet_name)
                # print ("Eplet_Poly: " + eplet_poly)
                # print ("Eplet_List_String: " + eplet_list_string)

                allele_epitope_string = aa_mm.getEpitope(allele_loctyp,position_list_int)

                # print ("Allele_Epitope_String: " + allele_epitope_string)

                if (eplet_list_string == allele_epitope_string):
                    HLA_Eplet_AlleleList[eplet_name].append(allele_loctyp)
                    HLA_Eplet_AlleleList[eplet_poly].append(allele_loctyp)
                    # print ("Allele " + allele_loctyp + " has eplet " + eplet_name + " with polymorphisms " + eplet_list_string)
                else:
                    continue
                    # print ("Allele " + allele_loctyp + " MISSING EPLET")   

            
# output allele list to file
eplet_allelelist_filename = "AA_eplet_to_alleles.csv"
eplet_allelelist_file = open(eplet_allelelist_filename, 'w')

eplet_allelelist_file.write("AA_Poly,Alleles\n")

for AA_poly in HLA_Eplet_AlleleList:
    allelelist = HLA_Eplet_AlleleList[AA_poly]
    allele_string = ','.join(allelelist)
    eplet_allelelist_file.write(AA_poly + "," + allele_string + "\n")

eplet_allelelist_file.close()
# print (HLA_AA_AlleleList[])
