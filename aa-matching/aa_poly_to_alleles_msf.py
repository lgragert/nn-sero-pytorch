#TODO - OOP approach to aa-matching and sequence information

import pandas as pd
#import aa_matching_msf as aa_mm
from aa_matching_msf import *
from collections import OrderedDict
from tqdm import tqdm

#aa_mm.main()
aa_mm = AAMatch(dbversion=3420, ungap=False)

def ungap(dataframe, refseq, loc):
    # the dashes will be put at the beginning of every set of possible
    # polymorphisms per residue
    ## this is to prevent all of the '-' characters from being sent to front
    #output_frame = output_frame.sort_index(axis=1, key = lambda x: (int(x[1])))
    i = 0
    j = 0
    new_cols = {}
    for column in dataframe:
        num = int(column[1:])
        if (column[0] == '-'):
            try:
                if dataframe.loc[refseq[loc]][column] == 1:
                    i += 1
                    j += 1
                    new_col = (str(num - i) + '_INS_' + str(j))
                    new_cols[column] = new_col
                else: 
                    j = 0
                    new_col = (str(num - i) + '_DEL_' + str(j))
                    new_cols[column] = new_col
            except KeyError:
                if dataframe.loc[refseq[loc][:-3]][column] == 1:
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


def getAApolys(loc, start, end):
    for allele_loctyp in aa_mm.HLA_full_allele:
        (allele_loc, allele_typ) = allele_loctyp.split('*')
        if (allele_loc != loc):
            continue
        # print (loc + '*' + allele_typ)

        for AA_pos in range(start,end):
            # print (AA_pos)
            
            side_chain = aa_mm.getAAposition(allele_loctyp,AA_pos)
            # if (side_chain == "-"):
            AA_poly = (side_chain + str(AA_pos+1))
            #AA_poly = (str(AA_pos) + side_chain)
            # print (allele_loctyp + " " + AA_poly)
            if AA_poly not in AA_polys.keys():
                AA_polys[AA_poly] = 0
                
    return AA_polys

# generate list of IMGT/HLA alleles that have each single AA polymorphism

for loc in aa_mm.ard_start_pos:
    print (loc)
    HLA_alleles = []
    AA_polys = {}
    # ARD starts at position 1 = Python index 0
    start = aa_mm.ard_start_pos[loc] - 1
    end = aa_mm.ard_end_pos[loc]

    AA_polys = getAApolys(loc, start, end)

    for allele_loctyp in tqdm(aa_mm.HLA_full_allele):
        (allele_loc, allele_typ) = allele_loctyp.split('*')
        if (allele_loc != loc):
            continue
        # print (loc + '*' + allele_typ)

        HLA_AA_polys = {}

        HLA_AA_polys = {key: 0 for key in AA_polys.keys()}
        
        XAA_polys = []

        # Generate list of polymorphisms in format AA# (# = position). Position adjusted to account for zero index.
        XAA_polys = [ aa_mm.getAAposition(allele_loctyp,AA_pos) + str(AA_pos + 1) for AA_pos in range(start,end) ] 

        HLA_AA_polys = {poly : 1 if poly in XAA_polys else 0 for poly in HLA_AA_polys.keys()}
            
        ## lambda function to sort on the number following the AA residue
        # (i.e. the residue position), then on the actual character for the
        # residue
        outie = OrderedDict(sorted(HLA_AA_polys.items(),
                                   key=lambda x: (int(x[0][1:]),str(x[0][0]))))
        outie['allele'] = allele_loctyp
        outie.move_to_end('allele', last=False)
        
        HLA_alleles.append(outie)


    
    output_frame = pd.DataFrame(HLA_alleles)
    output_frame = output_frame.set_index('allele')

    output_frame = ungap(output_frame, aa_mm.refseq_full, loc)

    output_frame.to_csv('./output/' + loc + '_AA_poly.csv', index=True)
    # print (HLA_AA_AlleleList[])