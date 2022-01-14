import pandas as pd
from aa_matching_msf import *
import re

aa_mm = AAMatch(dbversion=3420)

# testing the functions

allele1 = "A*02:01"
allele2 = "A*01:01"

position1 = 8

print (allele1)

AA = aa_mm.getAAposition(allele1,position1)
print ("AA at Position: " + str(position1) + ": " + AA)

AA_substring = aa_mm.getAAsubstring(allele1,1,11)
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

epitope_filename = "Epitope_Registry_DRB.txt"
epitope = pd.read_csv(epitope_filename,sep='\t')

print (epitope)

# ABC - Loci A, B, C,
# DRB - Loci DRB1, DRB3, DRB4, DRB5 (DP is dropped from definition)
# DQ - Loci DQA1, DQB1
# DP - Loci DPA1, DPB1

(epitope_db_prefix,epitope_db_suffix) = epitope_filename.split('.')
(ep,reg,epitope_db_name) = epitope_db_prefix.split('_') 

# each row is an eplet
for i in range (0,len(epitope)): 

    complex_poly = epitope['Polymorphic'][i]

    print (complex_poly)

    # eplet name cleanup
    complex_poly = complex_poly.replace(" ","") # remove spaces to concatenate polymorphisms
    complex_poly = complex_poly.replace("pairwith","") # 'pair with' means the same thing as spaces
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

        eplet_poly = "Eplet_" + epitope_db_name + '_' + '_'.join(eplet_list)

        j = j + 2


    print (complex_poly)
    print (aa_array)
    print (eplet_list)
    print (eplet_name)
    print (eplet_poly)
