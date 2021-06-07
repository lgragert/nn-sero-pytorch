#!/usr/bin/env python3
import pandas as pd
import aa_matching_msf as aa

aa_mm = aa.AAMatch(dbversion=3420, ungap=False)

def ungap(dataframe, refseq, loc):
    # the dashes will be put at the beginning of every set of possible
    # polymorphisms per residue
    ## this is to prevent all of the '-' characters from being sent to front
    #output_frame = output_frame.sort_index(axis=1, key = lambda x: (int(x[1])))
    i = 0
    j = 0
    new_cols = {}
    for column in dataframe:
        num = int(column)
        if dataframe.loc[refseq[loc]][column] == '-':
            i += 1
            j += 1
            new_col = (str(num - i) + '_INS_' + str(j))
            new_cols[column] = new_col
        else: 
            j = 0
            new_col = (str(num - i))
            new_cols[column] = new_col
    
    dataframe = dataframe.rename(columns=new_cols)
    return dataframe


def getAApolys(loc, start, end):
    HLA_alleles = []
    for allele_loctyp in aa_mm.HLA_full_allele:
        AA_polys = {}
        (allele_loc, allele_typ) = allele_loctyp.split('*')
        if (allele_loc != loc):
            continue
        # print (loc + '*' + allele_typ)
        AA_polys['allele'] = allele_loctyp
        AA_polys.update({str(AA_pos): aa_mm.getAAposition(allele_loctyp,AA_pos) for AA_pos in range(start,end)})
        HLA_alleles.append(AA_polys)

    return HLA_alleles


for loc in aa_mm.ard_start_pos:
    print(loc)

    start = aa_mm.ard_start_pos[loc]
    end = aa_mm.ard_end_pos[loc]
    print(end)
    
    HLA_alleles = getAApolys(loc, start, end)
    output_frame = pd.DataFrame(HLA_alleles)
    output_frame = output_frame.set_index('allele')
    
    output_frame = ungap(output_frame, aa_mm.refseq_full, loc)
    output_frame = pd.get_dummies(output_frame)
    output_frame.to_csv('./output/{}_AA_poly.csv'.format(loc), index=True)
