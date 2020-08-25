#TODO - OOP approach to aa-matching and sequence information

import pandas as pd
import numpy as np
import aa_matching_msf as aa_mm
from collections import OrderedDict

def ungap(dataframe, refseq, loc):
    # the dashes will be put at the beginning of every set of possible polymorphisms per residue
    ## this is to prevent all of the '-' characters from being sent to front
    #output_frame = output_frame.sort_index(axis=1, key = lambda x: (int(x[1])))
    i = 0
    j = 0
    new_cols = {}

    for column in dataframe:
        num = int(column[1:])
        if (column[0] == '-'):
            if dataframe.loc[refseq[loc]][column] == 1:
                i += 1
                j += 1
                new_col = (str(num - i) + '_INS_' + str(j))
                new_cols[column] = new_col
            else: 
                j = 0
                new_col = (str(num - i) + '_DEL_' + str(j))
                new_cols[column] = new_col
        else:
            if j == 0:
                new_col = (str(column[0]) + str(num - i))
                new_cols[column] = new_col
            else:
                new_col = (str(column[0]) + str(num - i) + '_INS_' + str(j))
                new_cols[column] = new_col

    
    dataframe = dataframe.rename(columns=new_cols)
    return dataframe

# generate list of IMGT/HLA alleles that have each single AA polymorphism

for loc in aa_mm.ard_start_pos:
    print (loc)
    HLA_alleles = []
    AA_polys = {}
    # print (hlaProteinOffset[hla])
    for allele_loctyp in aa_mm.HLA_seq:
        (allele_loc, allele_typ) = allele_loctyp.split('*')
        if (allele_loc != loc):
            continue
        # print (loc + '*' + allele_typ)
<<<<<<< Updated upstream:aa_matching/aa_poly_to_alleles_msf.py

        for AA_pos in range (aa_mm.ard_start_pos[loc],aa_mm.ard_end_pos[loc]):
            # print (AA_pos)
=======
        for AA_pos in range(aa_mm.ard_start_pos[loc],aa_mm.ard_end_pos[loc]):
>>>>>>> Stashed changes:aa-matching/aa_poly_to_alleles_msf.py
            side_chain = aa_mm.getAAposition(allele_loctyp,AA_pos)
            # if (side_chain == "-"):
            AA_poly = (side_chain + str(AA_pos))
            #AA_poly = (str(AA_pos) + side_chain)
            # print (allele_loctyp + " " + AA_poly)
            if AA_poly not in AA_polys.keys():
                AA_polys[AA_poly] = 0
                
    for Xallele_loctyp in aa_mm.HLA_seq:
        (Xallele_loc, Xallele_typ) = Xallele_loctyp.split('*')
        if (Xallele_loc != loc):
            continue
        # print (loc + '*' + allele_typ)

        HLA_AA_polys = {}
        
        for key in AA_polys.keys():
            HLA_AA_polys[key] = 0

        for XAA_pos in range(aa_mm.ard_start_pos[loc],aa_mm.ard_end_pos[loc]):
            Xside_chain = aa_mm.getAAposition(Xallele_loctyp,XAA_pos)
            # if (side_chain == "-"):

            XAA_poly = (Xside_chain + str(XAA_pos))
            #XAA_poly = (str(XAA_pos) + Xside_chain)
            
            for poly in HLA_AA_polys.keys():
                if poly == XAA_poly:
                    HLA_AA_polys[poly] = 1

            
<<<<<<< Updated upstream:aa_matching/aa_poly_to_alleles_msf.py
        outie = OrderedDict(sorted(HLA_AA_polys.items(), key = lambda x: (int(x[0][1:]))))
=======
        ## lambda function to sort on the number following the AA residue (i.e. the residue position), then on the actual character for the residue 
        outie = OrderedDict(sorted(HLA_AA_polys.items(), key = lambda x: ((int(x[0][1:])), str(x[0][0]))))
>>>>>>> Stashed changes:aa-matching/aa_poly_to_alleles_msf.py
        outie['allele'] = Xallele_loctyp
        outie.move_to_end('allele', last=False)
        
        HLA_alleles.append(outie)


    
    output_frame = pd.DataFrame(HLA_alleles)
<<<<<<< Updated upstream:aa_matching/aa_poly_to_alleles_msf.py
    output_frame.to_csv('aa_matching/output/' + loc + '_AA_poly.csv', index=False)
=======
    output_frame = output_frame.set_index('allele')

    output_frame = ungap(output_frame, refseq, loc)

    output_frame.to_csv('./output/' + loc + '_AA_poly.csv', index=True)
>>>>>>> Stashed changes:aa-matching/aa_poly_to_alleles_msf.py
    # print (HLA_AA_AlleleList[])