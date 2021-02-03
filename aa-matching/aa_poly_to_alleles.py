import pandas as pd
import aa_matching as aa_mm
from collections import defaultdict

# generate list of IMGT/HLA alleles that have each single AA polymorphism

HLA_AA_AlleleList = defaultdict(list)

for loc in aa_mm.ard_start_pos:
    print (loc)
    # print (hlaProteinOffset[hla])
    for allele_loctyp in aa_mm.HLA_seq:
        (allele_loc, allele_typ) = allele_loctyp.split('*')
        if (allele_loc != loc):
            continue
        # print (loc + '*' + allele_typ)

        for AA_pos in range (aa_mm.ard_start_pos[loc],aa_mm.ard_end_pos[loc]):
            # print (AA_pos)
            side_chain = aa_mm.getAAposition(allele_loctyp,AA_pos)
            # if (side_chain == "-"):

            AA_poly = (str(AA_pos) + side_chain)
            loc_AA_poly = loc + "_" + AA_poly
            # print (allele_loctyp + " " + AA_poly)
            HLA_AA_AlleleList[loc_AA_poly].append(allele_loctyp)
            

# output allele list to file
AA_allelelist_filename = "AA_poly_to_alleles.csv"
AA_allelelist_file = open(AA_allelelist_filename, 'w')

AA_allelelist_file.write("AA_Poly,Alleles\n")

for AA_poly in HLA_AA_AlleleList:
    allelelist = HLA_AA_AlleleList[AA_poly]
    allele_string = ','.join(allelelist)
    AA_allelelist_file.write(AA_poly + "," + allele_string + "\n")

# print (HLA_AA_AlleleList[])