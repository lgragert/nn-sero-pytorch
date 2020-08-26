import re
import pandas as pd
import aa_matching_msf as aa_mm
from collections import OrderedDict

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

HLA_seq = aa_mm.HLA_seq
for loc in aa_mm.ard_start_pos:
    locDict = { newKey: HLA_seq[newKey] for newKey in HLA_seq.keys() if (newKey.split('*')[0] == loc) }
    locFrame = pd.DataFrame.from_dict(locDict)
    locDict = {}
    locFrame = locFrame.transpose()
    locFrame = pd.get_dummies(locFrame, prefix_sep='')
    locFrame = locFrame.rename(mapper=(lambda x: (str(x[-1]) + str(int(x[:-1])+1))), axis=1)
    locFrame.to_csv('./output/' + loc + '_AA_poly.csv', index=True)