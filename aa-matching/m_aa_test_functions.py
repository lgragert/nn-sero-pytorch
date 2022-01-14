import pandas as pd
import aa_matching as aa_mm
import re
from collections import defaultdict

epitope_filename = "Epitope_Registry_ABC.txt"
epitope = pd.read_csv(epitope_filename,sep='\t')

print (epitope)

(epitope_db_prefix,epitope_db_suffix) = epitope_filename.split('.')
(ep,reg,epitope_db_name) = epitope_db_prefix.split('_') 

# each row is an eplet
for i in range (0,len(epitope)): 

    complex_poly = epitope['Polymorphic'][i]

    print (complex_poly)

    # eplet name cleanup
    complex_poly = complex_poly.replace(" ","") # remove spaces to concatenate polymorphisms
    complex_poly = complex_poly.replace("pairwith","") # 'pair with' means the same thing as spaces
    complex_poly = complex_poly.replace("(","") # remove opened parentheses
    complex_poly = complex_poly.replace(")","") # remove closed parantheses 
    if ('or' in complex_poly):  # remove eplet definitions for DP from DRB
        (complex_poly,post_or) = complex_poly.split('or')

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

        # parsing eplet_list to create position_list

        position_list = re.sub('([A-Z])',' ', complex_poly)
        position_list = position_list.replace(" ",",")
        position_list = position_list[:-1]

        eplet_poly = "Eplet_" + epitope_db_name + '_' + '_'.join(eplet_list)

        j = j + 2

    print (complex_poly)
    print (aa_array)
    print (eplet_list)
    print (position_list)
    print (eplet_name)
    print (eplet_poly)

# testing the functions

allele1 = "A*02:01"
allele2 = "A*01:01"

position1 = 15

print (allele1)

AA = aa_mm.getAAposition(allele1,position1)
print ("AA at Position: " + str(position1) + ": " + AA)

AA_substring = aa_mm.getAAsubstring(allele1,1,10)
print ("AAs from Position 1 to 10: " + AA_substring)

epitope_string = aa_mm.getEpitope(allele1,[1,3,5,8,10])
print ("Epitope for AA Position List [1,3,5,8,10]: " + epitope_string)

AA = aa_mm.getAAposition(allele1,44)
print ("A*02:01 position 44: " + AA)
AA = aa_mm.getAAposition(allele2,44)
print ("A*01:01 position 44: " + AA)

mismatched = aa_mm.isPositionMismatched(allele1,allele2,44)
print ("Are A*02:01 and A*01:01 mismatched at position 44?: " + str(mismatched))


AA = aa_mm.getAAposition(allele1,45)
print ("A*02:01 position 45: " + AA)
AA = aa_mm.getAAposition(allele2,45)
print ("A*01:01 position 45: " + AA)

mismatched = aa_mm.isPositionMismatched(allele1,allele2,45)
print ("Are A*02:01 and A*01:01 mismatched at position 45?: " + str(mismatched))

mm_count = aa_mm.count_AA_Mismatches_Allele(allele1,allele1,allele2,allele2,44)
print ("# of Mismatches between A*02:01 homozygous donor and A*01:01 homozygous recip at position 44: " + str(mm_count))

AA1_donor = "Y"
AA2_donor = "Y"
AA1_recip = "D"
AA2_recip = "D"
mm_count = aa_mm.count_AA_Mismatches(AA1_donor,AA2_donor,AA1_recip,AA2_recip)
print ("# of Mismatches between YY AAs in donor and DD in recip: " + str(mm_count))
AA1_recip = "Y"
mm_count = aa_mm.count_AA_Mismatches(AA1_donor,AA2_donor,AA1_recip,AA2_recip)
print ("# of Mismatches between YY AAs in donor and YD in recip: " + str(mm_count))

allele1_donor = "A*02:01"
allele2_donor = "A*01:01"
allele1_recip = "A*03:01"
allele2_recip = "A*66:01"

