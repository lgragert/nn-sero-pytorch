import pandas as pd
import aa_matching_msf as aa_mm
from collections import OrderedDict

# generate list of IMGT/HLA alleles that have each single AA polymorphism

refseq = {"A":"A*01:01", "B":"B*07:02", "C":"C*01:02", "DRB1":"DRB1*01:01", "DQB1":"DQB1*05:01", "DPB1":"DPB1*01:01"}

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

        for AA_pos in range (aa_mm.ard_start_pos[loc],aa_mm.ard_end_pos[loc]):
            # print (AA_pos)
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
            # print (AA_pos)
            Xside_chain = aa_mm.getAAposition(Xallele_loctyp,XAA_pos)
            # if (side_chain == "-"):

            XAA_poly = (Xside_chain + str(XAA_pos))
            #XAA_poly = (str(XAA_pos) + Xside_chain)
            
            for poly in HLA_AA_polys.keys():
                if poly == XAA_poly:
                    HLA_AA_polys[poly] = 1

            
        outie = OrderedDict(sorted(HLA_AA_polys.items(), key = lambda x: (int(x[0][1:]))))
        outie['allele'] = Xallele_loctyp
        outie.move_to_end('allele', last=False)
        
        HLA_alleles.append(outie)


    
    output_frame = pd.DataFrame(HLA_alleles)
    output_frame = output_frame.set_index('allele')

    i = aa_mm.ard_start_pos[loc]
    j = 0
    new_cols = {}

    for column in output_frame:
        if (column[0] == '-'):
            if output_frame.loc[refseq[loc]][column] == 1:
                new_col = (str(i) + '_insert_' + str(j))
                new_cols[column] = new_col
                j += 1
            else: 
                j = 1
                new_col = (str(i) + '_insert_' + str(j))
                new_cols[column] = new_col
        else:
            if output_frame.loc[refseq[loc]][column] == 1:
                new_col = (str(column[0]) + str(i))
                new_cols[column] = new_col
                i += 1
            else:
                new_col = (str(column[0]) + str(i))
                new_cols[column] = new_col

        
    output_frame = output_frame.rename(columns=new_cols)


    output_frame.to_csv('aa_matching/output/' + loc + '_AA_poly.csv', index=True)
    # print (HLA_AA_AlleleList[])