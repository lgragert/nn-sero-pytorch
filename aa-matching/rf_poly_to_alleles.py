import pandas as pd
import os
import os.path
import aa_matching_msf as aa
from collections import OrderedDict
from tqdm import tqdm

pathloc =  str(os.path.dirname(os.path.realpath(__file__)))
aa_mm = aa.AAMatch(dbversion=3420, ungap=False)

class Poly2Alleles:

    def __init__(self, one_hot=False):
        self.one_hot = one_hot
        self.main()

    def ungap(self, dataframe, refseq, loc):
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


    def getAApolys(self, loc, start, end):
        HLA_alleles = []
        for allele_loctyp in aa_mm.HLA_full_allele:
            AA_polys = {}
            (allele_loc, allele_typ) = allele_loctyp.split('*')
            if (allele_loc != loc):
                continue
            # print (loc + '*' + allele_typ)
            AA_polys['allele'] = allele_loctyp
            AA_polys.update({str(AA_pos): aa_mm.getAAposition(allele_loctyp,AA_pos) for AA_pos in range(start,end+1)})
            HLA_alleles.append(AA_polys)

        return HLA_alleles

    def main(self):
        for loc in aa_mm.ard_start_pos:
            print(loc)

            start = aa_mm.ard_start_pos[loc]
            end = aa_mm.ard_end_pos[loc]

            HLA_alleles = self.getAApolys(loc, start, end)
            output_frame = pd.DataFrame(HLA_alleles)
            output_frame = output_frame.set_index('allele')

            output_frame = self.ungap(output_frame, aa_mm.refseq_full, loc)

            if self.one_hot:

                features = list(HLA_alleles[0].keys())
                features = features.remove('allele')
                output_frame = pd.get_dummies(output_frame, columns=features)
                output_frame.to_csv('{}/OHoutput/{}_AA_poly.csv'.format(pathloc,loc), index=True)

            else:
                output_frame.to_csv('{}/output/{}_AA_poly.csv'.format(pathloc,loc), index=True)


Poly2Alleles(one_hot=True)