import pandas as pd
import aa_matching as aa_mm
import re
from collections import defaultdict

# extract data from file Epitope_Registry_DRB.txt

epitope_filename = "Epitope_Registry_ABC.txt"
# epitope_filename = "Epitope_Registry_DRB.txt"
epitope = pd.read_csv(epitope_filename,sep='\t')

print (epitope)

(epitope_db_prefix,epitope_db_suffix) = epitope_filename.split('.')
(ep,reg,epitope_db_name) = epitope_db_prefix.split('_') 

# generate list of IMGT/HLA alleles that correspond to eplets in the Eplet Registry

HLA_Eplet_AlleleList = defaultdict(list)

# for loc in ["DRB1"]:
# for loc in ["DRB1","DRB3","DRB4","DRB5"]:
for loc in ["A","B","C"]:
# for loc in aa_mm.ard_start_pos:
    print (loc)
    # print (hlaProteinOffset[hla])
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

            poly_array = complex_poly.split('|')

            print ("Poly Array: " + str(poly_array))

            aa_array = re.split('([A-Z])',complex_poly)

            eplet_name = "Eplet_" + epitope_db_name + '_' + epitope['Eplet'][i]
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

            position_string = re.sub('([A-Z])',' ', complex_poly)
            position_string = position_string.replace(" ",",")
            position_string = position_string[:-1]
            
            position_list = position_string.split(",")

            print ("Complex_Poly: " + complex_poly)
            print ("AA_array: " + str(aa_array))
            print ("Eplet_list: " + str(eplet_list))
            print ("Position_String: " + position_string)
            print ("Position_List: " + str(position_list))

            position_list_int = [int(i) for i in position_list]

            eplet_list_string = '_'.join(eplet_list)
            eplet_poly = "Eplet_" + epitope_db_name + '_' + '_'.join(eplet_list)
        
            print ("Position_List: " + str(position_list_int))
            print ("Eplet_Name: " + eplet_name)
            print ("Eplet_Poly: " + eplet_poly)
            print ("Eplet_List_String: " + eplet_list_string)

            allele_epitope_string = aa_mm.getEpitope(allele_loctyp,position_list_int)


            print ("Allele_Epitope_String: " + allele_epitope_string)

            if (eplet_list_string == allele_epitope_string):
                HLA_Eplet_AlleleList[eplet_name].append(allele_loctyp)
                print ("Allele " + allele_loctyp + " has eplet " + eplet_name + " with polymorphisms " + eplet_list_string)
            else:
                print ("Allele " + allele_loctyp + " MISSING EPLET")   

            
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

# testing the functions

allele1 = "A*02:01"
position1 = 15

aa_position_list = []
for pos in position_list:
    aa_position_list.append(pos)
    print (aa_position_list)

print (allele1)

AA = aa_mm.getAAposition(allele1,position1)
print ("AA at Position: " + str(position1) + ": " + AA)

epitope_string = aa_mm.getEpitope(allele1,[1,3,5,8,10])
print ("Epitope for AA Position List [1,3,5,8,10]: " + epitope_string)


