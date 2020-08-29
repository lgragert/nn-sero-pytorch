import re
import pandas as pd
import aa_matching_msf as aa_mm
from collections import OrderedDict


def ungap(dataframe, refseq, loc):
    # the dashes will be put at the beginning of every set of possible polymorphisms per residue
    ## this is to prevent all of the '-' characters from being sent to front
    # output_frame = output_frame.sort_index(axis=1, key = lambda x: (int(x[1])))
    i = 0
    j = 0
    new_cols = {}

    for column in dataframe:
        num = int(column[1:])
        if (column[0] == '-'):
            if dataframe.loc[refseq[loc]][column] == 1:
                i += 1
                j += 2
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

def toBinary(string):
    string = ''.join(format(ord(x), 'b')) for x in string
    return string

#from StackOverflow:  https://stackoverflow.com/questions/52452911/finding-all-positions-of-a-character-in-a-string
def findIns(sequence):
    seqIns = []
    idx = sequence.find('-')
    if (idx != -1):
        seqIns.append(idx)
    while (idx != -1):
        yield idx
        idx = sequence.find('-', idx+1)
        seqIns.append(idx)
    return seqIns

def impute(locDict, refseq):
    seqIns = findIns(locDict[refseq])
    for key in locDict.keys():
        locDict[key] = toBinary(key)
    binFrame = pd.DataFrame.from_dict(locDict)


refseq = {
    "A" : "A*01:01",
    "B" : "B*07:02",
    "C" : "C*01:02",
    "DRB1" : "DRB1*01:01",
    "DRB3" : "DRB3*01:01",
    "DRB4" : "DRB4*01:01",
    "DRB5" : "DRB5*01:01",
    "DQA1" : "DQA1*01:01",
    "DQB1" : "DQB1*05:01",
    "DPA1" : "DPA1*01:03",
    "DPB1" : "DPB1*01:01",
    }
#
# refseq = {
#     "A" : "A*01:01:01:01",
#     "B" : "B*07:02:01:01",
#     "C" : "C*01:02:01:01",
#     "DRB1" : "DRB1*01:01:01:01",
#     "DRB3" : "DRB3*01:01:02:01",
#     "DRB4" : "DRB4*01:01:01:01",
#     "DRB5" : "DRB5*01:01:01:01",
#     "DQA1" : "DQA1*01:01:01:01",
#     "DQB1" : "DQB1*05:01:01:01",
#     "DPA1" : "DPA1*01:03:01:01",
#     "DPB1" : "DPB1*01:01:01:01",
#     }
#

HLA_seq = aa_mm.HLA_seq
for loc in aa_mm.ard_start_pos:
    print("Processing locus " + loc + "...")
    locDict = { newKey: (HLA_seq[newKey])[aa_mm.ard_start_pos[loc]-1:aa_mm.ard_end_pos[loc]] for newKey in HLA_seq.keys() if (newKey.split('*')[0] == loc) }
    locFrame = pd.DataFrame.from_dict(locDict)
    locDict = {}
    locFrame = locFrame.transpose()
    locFrame = pd.get_dummies(locFrame, prefix_sep='')
    locFrame = locFrame.rename(mapper=(lambda x: (str(x[-1]) + str(int(x[:-1])+1))), axis=1)
    locFrame.index.names = ['allele']
    locFrame = ungap(locFrame, refseq, loc)
    locFrame.to_csv('./output/' + loc + '_AA_poly.csv', index=True)
    print("Done")