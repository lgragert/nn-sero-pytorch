import re
import pandas as pd
import numpy as np
import aa_matching_msf as aa_mm
from collections import OrderedDict

def ungap(dataframe, refseq, loc):
    # the dashes will be put at the beginning of every set of possible polymorphisms per residue
    ## this is to prevent all of the '-' characters from being sent to front
    # output_frame = output_frame.sort_index(axis=1, key = lambda x: (int(x[1])))
    i = 0
    j = 0
    new_cols = {}
    tNum = 0

    for column in dataframe:
        num = int(column[1:])
        if (column[0] == '-'):
            if dataframe.loc[refseq[loc]][column] == 1:
                tNum = num
                i += 1
                j += 1
                new_col = (str(num - i) + '_INS_' + str(j))
                new_cols[column] = new_col
        else:
            if tNum == num:
                new_col = (str(column[0]) + str(num - i) + '_INS_' + str(j))
                new_cols[column] = new_col
            else:
                new_col = (str(column[0]) + str(num - i))
                new_cols[column] = new_col

    dataframe = dataframe.rename(columns=new_cols)
    return dataframe

def toBinary(string):
    string = ''.join(format(ord(x), 'b') for x in string)
    return string

#from StackOverflow:  https://stackoverflow.com/questions/52452911/finding-all-positions-of-a-character-in-a-string
def findIns(sequence):
    seqIns = []
    idx = sequence.find('-')
    if (idx != -1):
        seqIns.append(idx)
    while (idx != -1):
        idx = sequence.find('-', idx+1)
        if (idx != -1):
            seqIns.append(idx)
    return seqIns

def checkComplete(sequence, seqIns):
    checkIns = []
    idx = sequence.find('-')
    if (idx != -1):
        checkIns.append(idx)
    while (idx != -1):
        idx = sequence.find('-', idx+1)
        if (idx != -1):
            checkIns.append(idx)
    checkIns = [x for x in checkIns if x not in seqIns]
    return checkIns

def sumHam(binNum):
    sum = 0
    binNum = str(binNum)
    for char in binNum:
        sum += int(char)
    return sum

def impute(locDict, refseq):
    seqIns = findIns(locDict[refseq])
    replacePos = {}
    binDict = {}
    for key in locDict.keys():
        replacePos[key] = checkComplete(locDict[key], seqIns)
        binDict[key] = toBinary(locDict[key])
    for rKey in replacePos.keys():
        rDist = {}
        if (len(replacePos[rKey]) != 0):
            print("Imputing peptide sequence for allele " + str(rKey))
            hDict = {hKey: binDict[hKey] for hKey in
                       binDict.keys() if (hKey.split(':')[0] == rKey.split(':')[0]) and len(replacePos[hKey]) == 0}
            for binKey in hDict.keys():
                if binKey != rKey:
                    summed = int(binDict[rKey], 2) ^ int(hDict[binKey], 2)
                    rDist[binKey] = bin(summed)[2:].zfill(len(locDict[rKey]))
                    rDist[binKey] = sumHam(rDist[binKey])
                else:
                    next
            nNearest = 100000
            nearest = "NA"
            for near in rDist.keys():
                nNear = int(rDist[near])
                if (nNear < nNearest):
                    nNearest = nNear
                    nearest = near
                else:
                    next
            for rVal in replacePos[rKey]:
                if nearest != "NA":
                    locDict[rKey] = locDict[rKey][:rVal] + locDict[nearest][rVal] + locDict[rKey][rVal+1:]
    return locDict

refseq = aa_mm.refseq
HLA_seq = aa_mm.HLA_seq
for loc in aa_mm.ard_start_pos:
    print("Processing locus " + loc + "...")
    locDict = { newKey: str(HLA_seq[newKey].seq) for newKey in HLA_seq.keys() }
    newDict = { locKey: locDict[locKey][(aa_mm.ard_start_pos[loc] - 1):(aa_mm.ard_end_pos[loc])] for locKey in
               locDict.keys() if (locKey.split('*')[0] == loc) }
    locDict = newDict
    del(newDict)
    imputed = impute(locDict, refseq[loc])
    repDict = { repKey: list(imputed[repKey]) for repKey in imputed.keys() }
    del(imputed)
    repFrame = pd.DataFrame.from_dict(repDict)
    repFrame = repFrame.transpose()
    repFrame = pd.get_dummies(repFrame, prefix_sep='')
    repFrame = repFrame.rename(mapper=(lambda x: (str(x[-1]) + str(int(x[:-1]) + 1))), axis=1)
    repFrame.index.names = ['allele']
    repFrame = ungap(repFrame, refseq, loc)
    repFrame.to_csv('./imputed/' + loc + '_imputed_poly.csv', index=True)
    print("Done with locus " + loc)
