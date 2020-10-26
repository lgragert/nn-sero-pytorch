#!/usr/bin/env python
###############################################################################
#   SCRIPT NAME:    impute_nuc.py
#   DESCRIPTION:    Module for inference of HLA nucleotide sequences
#   OUTPUT:
#   DATE:           September 01, 2020
#   LATEST UPDATE:  October 5, 2020
#   UPDATE REASON:  Distance matrix addition for computational efficiency
#   AUTHOR:         Giovanni Biagini (dbiagini@tulane.edu ; GitHub: gbiagini)
#   PI:             Loren Gragert, Ph.D.
#   ORGANIZATION:   Tulane University School of Medicine
#   NOTES:          Written upon suggestions from Martin Maiers (NMDP)
###############################################################################

import pandas as pd
import numpy as np
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
    # this is to prevent all of the '-' characters from being sent to front
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

# translate the nucleotide CDS into the amino acid sequence for all alleles


def translate_nuc(nuc_dict, seqIns):
    translated = {}
    incorrect = []
    print("Tranlsation in progress...")
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
            new_sequel = [new_sequel[x]
                          for x in range(0, len(new_sequel))
                          if (x not in seqIns)]
            nuc_seq = SeqRecord(Seq(''.join(new_sequel)))
            aa_seq = SeqRecord(seq=nuc_seq.seq.translate(cds=False,
                                                         to_stop=True))
        translated[record] = aa_seq
    return translated, incorrect

# add back gaps that are present in the reference (and other) alleles so that
# arrays are the same length


def correction(incorrect, translated, aaIns):
    for each in incorrect:
        for insert in aaIns:
            translated[each] = translated[each][:insert] + '-' + \
                translated[each][insert:]
    return translated

# complete null alleles so that arrays are the same length


def finish_null(refseq, repDict):
    removal = []
    length = len(repDict[refseq])
    for entry in repDict.keys():
        separator = entry.find("*")
        if any(x in entry[separator:] for x in suffixes):
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


# algorithm to generate distance matrices for AA sequences and output to csv
# files

# self is a Boolean to describe if the allele should be compared to itself.
# Set to False when performing imputation
def distmat(locDict, binDict, hDict, self=True):
    matDist = {}
    print("Generating distance matrix...")
    for bKey in tqdm(binDict.keys()):
        rDist = {}
        for binKey in hDict.keys():
            if self:
                xoresult = int(binDict[bKey], 2) ^ int(hDict[binKey], 2)
                rDist[binKey] = bin(xoresult)[2:].zfill(len(locDict[bKey]))
                rDist[binKey] = sumHam(rDist[binKey])
            else:
                if (bKey != binKey):
                    xoresult = int(binDict[bKey], 2) ^ int(hDict[binKey], 2)
                    rDist[binKey] = bin(xoresult)[2:].zfill(len(locDict[bKey]))
                    rDist[binKey] = sumHam(rDist[binKey])
        matDist[bKey] = rDist
    disFrame = pd.DataFrame.from_dict(matDist)
    return disFrame

# nearest 10 vote based on distance matrix
# TODO (gbiagini) - this function is not yet complete
def nearest10(loc, disFrame, rPos):
    with open("./data/nearest/" + loc + "_topten.txt", "w+") as handle:
        handle.write("NEAREST 10 NEIGHBORS FOR EACH IMPUTED ALLELE + \n")
        handle.write("Format for neighbors is $ALLELE (HamDist for NA Seq)\n")
        print("Nearest 10 imputation...")
        for rKey in tqdm(rPos.keys()):
            disFrame = disFrame.sort_values(by=rKey, axis=1, ignore_index=False)
            # pull top 10 closest alleles (skip index 0 --> "Unnamed")
            handle.write("Imputed allele:\t\t\t\t" + rKey + "\n")
            topten = {}
            i = 0
            for n in list(disFrame.columns[0:10]):
                i += 1
                topten[n] = disFrame.loc[rKey][n]
                handle.write("Neighbor #" + str(i) + ":\t\t\t\t" + n + " (" +
                             str(topten[n]) + ")\n")

            # infers sequence from nearest 10 neighbor vote
            # TODO (gbiagini) - Weight sequence contribution for each
            #  position by Hamming distance
            for rVal in rPos[rKey]:
                repseq = {}
                for neighbor in topten.keys():
                    n_acid = locDict[neighbor][rVal]
                    if n_acid not in repseq.keys():
                        repseq[n_acid] = -abs(float(topten[neighbor]))
                    else:
                        repseq[n_acid] -= abs(float(topten[neighbor]))
                new_val = max(repseq, key = lambda x: repseq[x])
                locDict[rKey] = locDict[rKey][:rVal] + \
                                new_val + locDict[rKey][rVal + 1:]
    return locDict

def impute(loc, locDict, refseq, aaDict):
    seqIns = findIns(locDict[refseq])
    aaIns = findIns(str(aaDict[refseq].seq))
    binFrame = pd.read_csv("./data/nadist/" + loc + ".csv", index_col=0)
    binDict = binFrame.to_dict()

    # check to see if any new alleles have been added to the dataset
    # if so, update distance matrix
    # if not, use previously computed distance matrix
    if (set(binDict.keys()) != set(locDict.keys())):
        # TODO (gbiagini) - modify this to simply update previous distance
        #  matrix, rather than re-do the entire computation
        print("New alleles detected - must generate distance matrix!")
        replacePos = {}
        binDict = {}
        with open("./data/nabin/" + loc + ".tsv", "w+") as handle:
            handle.write("Allele\t\tBinary Nucleotide Sequence\n")
            for key in locDict.keys():
                replacePos[key] = checkComplete(locDict[key], seqIns)
                binDict[key] = toBinary(locDict[key])
                handle.write(key + '\t\t' + binDict[key] + '\n')
        disFrame = distmat(locDict, binDict, binDict)
        disFrame.to_csv("./data/nadist/" + loc + ".csv")
        del(disFrame)

        hDict = {hKey: binDict[hKey] for hKey in binDict.keys() if len(
            replacePos[hKey]) == 0}
        rPos = {r: replacePos[r] for r in replacePos.keys() if
                (len(replacePos[r]) != 0)}
        del(replacePos)
        disFrame = distmat(locDict, binDict, hDict, self=False)
        disFrame.to_csv("./data/compdist/" + loc + ".csv")
    else:
        replacePos = {}
        for key in locDict.keys():
                replacePos[key] = checkComplete(locDict[key], seqIns)
        rPos = {r: replacePos[r] for r in replacePos.keys() if
                (len(replacePos[r]) != 0)}
        del (replacePos)
        # TODO (gbiagini) - Need to fix this for situations in which the
        #  distance matrix is regenerated (make index_col = 0)
        disFrame = pd.read_csv("./data/compdist/" + loc + ".csv", index_col=0)

    disFrame = disFrame.transpose()
    locDict = nearest10(loc, disFrame, rPos)
    locDict, incorrect = translate_nuc(locDict, seqIns)
    binDict = {}
    with open("./data/nabini/" + loc + ".tsv", "w+") as handle:
        handle.write("Allele\t\tBinary Nucleotide Sequence (with Imputation)\n")
        for key in locDict.keys():
            binDict[key] = toBinary(locDict[key])
            handle.write(key + '\t\t' + binDict[key] + '\n')

    del(disFrame)
    # generation of distance matrix for imputed sequences
    #disFrame = distmat(locDict, binDict, binDict)
    #disFrame.to_csv("./data/nadisti/" + loc + ".csv")
    #del(disFrame)
    translated = correction(incorrect, locDict, aaIns)
    return translated

# apply hlaProteinOffset and then limit to antigen recognition domain
# !!IMPORTANT!! to use nuc_mm and NOT aa_mm due to modified offsets

def post_trans_mod(repDict, loc):
    for each in repDict.keys():
        repDict[each] = repDict[each][nuc_mm.hlaProteinOffset[loc]:]
        repDict[each] = repDict[each][nuc_mm.ard_start_pos[
            loc]:nuc_mm.ard_end_pos[loc]]
    return repDict


aaDict = aa_mm.HLA_seq
refseq = nuc_mm.refseq
HLA_seq = nuc_mm.HLA_seq
# for loc in ["A", "B", "C", "DPB1", "DRB1", "DQB1"]:
for loc in nuc_mm.refseq:
    print("Processing locus " + loc + "...")
    locDict = {newKey: str(HLA_seq[newKey].seq) for newKey in HLA_seq.keys()}
    newDict = {locKey: locDict[locKey] for locKey in locDict.keys() if (
        locKey.split('*')[0] == loc)}
    locDict = newDict
    del (newDict)
    imputed = impute(loc, locDict, refseq[loc], aaDict)
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
    print(repFrame)
    repFrame = ungap(repFrame, refseq, loc)
    if loc == "DRB1":
        repFrame.insert(1, "ZZ1", 0)
        repFrame.insert(2, "ZZ2", 0)
    repFrame.to_csv('./imputed/' + loc + '_imputed_poly.csv', index=True)
    print("Done with locus " + loc)
