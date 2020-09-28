#!/usr/bin/env python
###############################################################################
#   SCRIPT NAME:    impute_nuc.py
#   DESCRIPTION:    Module for inference of HLA nucleotide sequences
#   OUTPUT:
#   DATE:           September 01, 2020
#   AUTHOR:         Giovanni Biagini (dbiagini@tulane.edu ; GitHub: gbiagini)
#   PI:             Loren Gragert, Ph.D.
#   ORGANIZATION:   Tulane University School of Medicine
#   NOTES:          Written upon suggestions from Martin Maiers (NMDP)
###############################################################################

import pandas as pd
import nuc_matching_msf as nuc_mm
import aa_matching_msf as aa_mm
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from tqdm import tqdm

suffixes = ["L", "S", "C", "A", "Q", "N"]

def ungap(dataframe, refseq, loc):
    # the dashes will be put at the beginning of every set of possible
    # polymorphisms per residue
    ## this is to prevent all of the '-' characters from being sent to front
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


# converts sequence string to binary
# needed for XOR comparing binary sequences
def toBinary(string):
    string = ''.join(format(ord(x), 'b') for x in string)
    return string

# returns list of indexes for dash characters in all sequences
def findIns(sequence):
    seqIns = []
    idx = sequence.find('-')
    if (idx != -1):
        seqIns.append(idx)
    while (idx != -1):
        idx = sequence.find('-', idx + 1)
        if (idx != -1):
            seqIns.append(idx)
    return seqIns


# removes indexes if also in refseq[loc]
def checkComplete(sequence, seqIns):
    checkIns = []
    idx = sequence.find('-')
    if (idx != -1):
        checkIns.append(idx)
    while (idx != -1):
        idx = sequence.find('-', idx + 1)
        if (idx != -1):
            checkIns.append(idx)
    checkIns = [x for x in checkIns if x not in seqIns]
    return checkIns


# Sums result of XOR operation for nearest neighbor determination
def sumHam(binNum):
    sum = 0
    binNum = str(binNum)
    for char in binNum:
        sum += int(char)
    return sum

def translate_nuc(nuc_dict, seqIns):
    translated = {}
    incorrect = []
    for record in tqdm(nuc_dict.keys()):
        nuc_seq = SeqRecord(Seq(nuc_dict[record]))
        try:
            aa_seq = SeqRecord(seq=nuc_seq.seq.translate(cds=False,
                                                         to_stop=True))
        except TranslationError:
            incorrect.append(record)
            new_seq = str(nuc_seq.seq)
            new_sequel = []
            new_sequel[:] = new_seq
            new_sequel = [ new_sequel[x] for x in range(0,len(new_sequel)) if (x
                           not in seqIns) ]
            nuc_seq = SeqRecord(Seq(''.join(new_sequel)))
            aa_seq = SeqRecord(seq=nuc_seq.seq.translate(cds=False,
                                                         to_stop=True))
            print("Sequence for allele " + record + ": " + aa_seq.seq)
        translated[record] = aa_seq
    return translated, incorrect

def correction(incorrect, translated, aaIns):
    for each in incorrect:
        for insert in aaIns:
            translated[each] = translated[each][:insert] + '-' + \
                               translated[each][insert:]
    return translated

def finish_null(refseq, repDict):
    removal = []
    length = len(repDict[refseq])
    for entry in repDict.keys():
        if any(x in entry[1:] for x in suffixes):
            if (entry[-1] == "N"):
                trunk = len(repDict[entry])
                diff = length - trunk
                if (diff > 0):
                    filler = "*" * diff
                    fillist = list(filler)
                    fix = repDict[entry] + fillist
                    repDict[entry] = fix
            else:
                removal.append(entry)
    for bad in removal:
        del repDict[bad]
    return repDict

# testing function for validation of the translation approach
# def test_translate(sequence, rKey, seqIns):
#     nuc_seq = SeqRecord(Seq(sequence))
#     try:
#         aa_seq = SeqRecord(seq=nuc_seq.seq.translate(cds=True))
#         print(aa_seq)
#     except TranslationError:
#         new_seq = str(nuc_seq.seq)
#         new_sequel = []
#         new_sequel[:] = new_seq
#         new_sequel = [ new_sequel[x] for x in range(0,len(new_sequel)) if (x
#                        not in seqIns) ]
#         nuc_seq = SeqRecord(Seq(''.join(new_sequel)))
#         aa_seq = SeqRecord(seq=nuc_seq.seq.translate(cds=True))
#         print("Sequence for allele " + rKey + ": " + aa_seq.seq)
#     return

def impute(locDict, refseq, aaDict):
    seqIns = findIns(locDict[refseq])
    aaIns = findIns(str(aaDict[refseq].seq))
    replacePos = {}
    binDict = {}
    for key in locDict.keys():
        replacePos[key] = checkComplete(locDict[key], seqIns)
        binDict[key] = toBinary(locDict[key])
    for rKey in tqdm(replacePos.keys()):
        rDist = {}
        if (len(replacePos[rKey]) != 0):
            print("Imputing nucleotide sequence for allele " + str(rKey))
            # TODO (gbiagini) - difference accumulation - possible sorting
            hDict = {hKey: binDict[hKey] for hKey in
                     binDict.keys() if len(replacePos[hKey]) == 0}
            for binKey in hDict.keys():
                if binKey != rKey:
                    xoresult = int(binDict[rKey], 2) ^ int(hDict[binKey], 2)
                    rDist[binKey] = bin(xoresult)[2:].zfill(len(locDict[rKey]))
                    rDist[binKey] = sumHam(rDist[binKey])
            # arbitrary nearest value to be overwritten
            nNearest = 100000
            # placeholder for nearest allele
            nearest = "NA"
            # return allele closest to incomplete sequence
            for near in rDist.keys():
                nNear = int(rDist[near])
                if (nNear < nNearest):
                    nNearest = nNear
                    nearest = near
            if rKey in ["A*26:03", "C*03:23", "C*03:46", "DPB1*35:01",
                        "DRB1*04:20", "DRB1*05:13", "DQB1*06:06"]:
                print(rKey + " nearest neighbor: " + nearest)
            # infers sequence from nearest neighbor
            for rVal in replacePos[rKey]:
                if nearest != "NA":
                    locDict[rKey] = locDict[rKey][:rVal] + \
                                    locDict[nearest][rVal] + \
                                    locDict[rKey][rVal + 1:]
		# Here to test functionality of translation
        # elif (len(replacePos[rKey]) == 0):
        #     test_translate(locDict[rKey], rKey, seqIns)
    locDict, incorrect = translate_nuc(locDict, seqIns)
    translated = correction(incorrect, locDict, aaIns)
    return translated

# apply hlaProteinOffset and then limit to antigen recognition domain
def post_trans_mod(repDict, loc):
    for each in repDict.keys():
        repDict[each] = repDict[each][aa_mm.hlaProteinOffset[loc]:]
        repDict[each] = repDict[each][aa_mm.ard_start_pos[loc]:aa_mmard_end_pos[
            loc]]
    return repDict

aaDict = aa_mm.HLA_seq
refseq = nuc_mm.refseq
HLA_seq = nuc_mm.HLA_seq
#for loc in nuc_mm.refseq:
for loc in ["C", "DPB1", "DRB1", "DQB1"]:
    print("Processing locus " + loc + "...")
    locDict = {newKey: str(HLA_seq[newKey].seq) for newKey in HLA_seq.keys()}
    # TODO (gbiagini) - Removing the indexes that limit this to the antigen
    #  recognition domain. I'll need to reintroduce this in the future.
    newDict = {locKey: locDict[locKey] for locKey in locDict.keys() if (
                locKey.split('*')[0] == loc)}
    locDict = newDict
    del (newDict)
    imputed = impute(locDict, refseq[loc], aaDict)
    # creates list from sequence strings for Pandas dataframe
    repDict = {repKey: list(imputed[repKey]) for repKey in imputed.keys()}
    del (imputed)
    repDict = finish_null(refseq[loc], repDict)
    repDict = post_trans_mod(repDict, loc)
    repFrame = pd.DataFrame.from_dict(repDict)
    repFrame = repFrame.transpose()
    repFrame = pd.get_dummies(repFrame, prefix_sep='')
    # lambda function to rename columns with string character first then
    # position
    repFrame = repFrame.rename(
        mapper=(lambda x: (str(x[-1]) + str(int(x[:-1]) + 1))), axis=1)
    repFrame.index.names = ['allele']
    repFrame = ungap(repFrame, refseq, loc)
    if loc == "DRB1":
        repFrame.insert(1, "ZZ1", 0)
        repFrame.insert(2, "ZZ2", 0)
    repFrame.to_csv('./imputed/' + loc + '_imputed_poly.csv', index=True)
    print("Done with locus " + loc)
