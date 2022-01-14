import pandas as pd
import numpy as np

loci = ['A', 'B', 'C', 'DPB1', 'DQB1', 'DRB1']

base_dir = './randomforest/'

def clear(cell):
    #cell = [x for x in cell if x.isdigit()]
    cell = cell.replace('[','')
    cell = cell.replace(']','')
    cell = cell.replace("'",'')
    cell = cell.replace(" ",'')
    cell = cell.split(',')
    return cell

for loc in loci:
    compframe = pd.read_csv(base_dir + './comparison/' + loc + '_compfile.csv')
    locdict = {}

    olds = [clear(x) for x in compframe['Old Assignment']]
    for each in olds:
        for x in each:
            if x not in locdict.keys():
                locdict[x] = 1
            else:
                locdict[x] += 1
    

    resdict = {y : [locdict[y]] for y in locdict.keys()}

    resFrame = pd.DataFrame.from_dict(resdict)
    resFrame = resFrame.transpose()
    resFrame.to_csv(base_dir + "comparison/" + loc + "_compser.csv", index=True)