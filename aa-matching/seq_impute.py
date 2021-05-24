#!/usr/bin/env python3
###############################################################################
#   SCRIPT NAME:    seq_impute.py
#   DESCRIPTION:    Module for inference of HLA nucleotide sequences
#   OUTPUT:
#   DATE:           May 01, 2021
#   LATEST UPDATE:  
#   UPDATE REASON:  
#   AUTHOR:         Giovanni Biagini (dbiagini@tulane.edu ; GitHub: gbiagini)
#   PI:             Loren Gragert, Ph.D.
#   ORGANIZATION:   Tulane University School of Medicine
#   NOTES:          Formerly nuc_impute.py
###############################################################################


import math
import re
import os
import os.path
from os import wait
import pandas as pd
import numpy as np
import aa_matching_msf as aa
import nuc_matching_msf as na
from pathlib import Path
from collections import defaultdict
from functools import reduce
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

pathloc = str(os.getcwd()) + '/'

aa_mm = aa.AAMatch(dbversion=3400)
nuc_mm = na.NucMatch(dbversion=3400)


class Impute():

    suffixes = ["L", "S", "C", "A", "Q", "N"]

    def __init__(self, suffixes=suffixes):
        self.suffixes = suffixes    
        self.main()
        
    # removes indexes if also in refseq[loc]
    def checkComplete(self, sequence, seqIns):
        checkIns = [m.start() for m in re.finditer('-', sequence)]
        checkIns = [x for x in checkIns if x not in seqIns]
        return checkIns


    # translate the nucleotide CDS into the amino acid sequence for all alleles
    def translate_nuc(self, nuc_dict, seqIns):
        translated = {}
        incorrect = []
        print("Tranlsation in progress...")
        for record in nuc_dict.keys():
            nuc_seq = SeqRecord(Seq(nuc_dict[record]))
            try:
                aa_seq = SeqRecord(seq=nuc_seq.seq.translate(cds=False,
                                                            to_stop=True))
            except TranslationError:
                incorrect.append(record)
                new_seq = []
                new_seq[:] = str(nuc_seq.seq)
                new_seq = [new_seq[x]
                            for x in range(0, len(new_seq))
                            if (x not in seqIns)]
                new_sequel = ['N' for y in new_sequel if y == '-']
                nuc_seq = SeqRecord(Seq(''.join(new_sequel)))
                aa_seq = SeqRecord(seq=nuc_seq.seq.translate(cds=False,
                                                            to_stop=True))
            translated[record] = aa_seq
        return translated, incorrect

    # add back gaps that are present in the reference (and other) alleles so that
    # arrays are the same length
    def correction(self, incorrect, translated, aaIns):
        for each in incorrect:
            translated = {each: "{}{}{}".format(translated[each][:insert],'-',translated[each][insert:]) for insert in aaIns}
        return translated

    # complete null alleles so that arrays are the same length
    def finish_null(self, refseq, repDict):
        suffixes = self.suffixes
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


    # Normalized Hamming distance is based on the number of sequence differences
    #   divided by the number of sequence positions that are defined in both
    #   sequences.
    def distmatrix(self, locDict):
        distance = {}
        print("Generating distance matrix...")
        for allele1 in locDict.keys():
            seq_dist = {}
            seq1 = locDict[allele1]
            for allele2 in locDict.keys():
                # needed for sequences with no overlaps in complete sequence
                if (allele1 != allele2):
                    seq2 = locDict[allele2]
                    seq_dist[allele2] = sum(1 for w, x in zip(seq1, seq2) if (w != '-') and (x != '-') and (w != x))/sum(1 for y,z in zip(seq1, seq2) if (x != '-') and (y != '-'))
                else:
                    seq_dist[allele2] = 0
            distance[allele1] = seq_dist
        disFrame = pd.DataFrame.from_dict(distance)
        return disFrame

    # Recursive algorithm for nearest 10 neighbor voting for sequence inference.    
    def seqvote(self, rVal, rKey, disFrame, locDict, start):
        rVal = rVal - 1
        votes = {}
        new_val = '-'

        for neighbor in list(disFrame.columns[start:start+10]):
            n_acid = locDict[neighbor][rVal]
            if (n_acid not in votes.keys() and n_acid != "-"):
                votes[n_acid] = abs(float(disFrame.loc[rKey][neighbor]))
            elif (n_acid == "-"):
                continue

            # for additional votes for a single position, add value
            # to increase the likelihood of selection
            else:
                votes[n_acid] += abs(float(disFrame.loc[rKey][neighbor]))

            try:
                new_val = max(votes, key = lambda x: votes[x])
            except ValueError:
                start += 10
                return self.seqvote(rVal, rKey, disFrame, locDict, start)
        return new_val

    # nearest 10 vote based on distance matrix
    def nearest10(self, loc, locDict, disFrame, missing_pos):
        #filter out null alleles
        disFrame = disFrame.filter(regex="[ABCD][PQR]?[AB]?\d?\*[\d\:]*(?!.*\w)", axis=1)
        with open(pathloc + "data/nearest/{}_topten.txt".format(loc), "w+") as handle:
            handle.write("NEAREST 10 NEIGHBORS FOR EACH IMPUTED ALLELE + \n")
            handle.write("Format for neighbors is $ALLELE (HamDist for NA Seq)\n")
            print("Nearest 10 imputation...")
            for inc_allele in missing_pos.keys():
                disFrame = disFrame.sort_values(by=inc_allele, axis=1, ignore_index=False)
                # pull top 10 closest alleles
                handle.write("Imputed allele:\t\t\t\t" + inc_allele + "\n")
                if str(disFrame.columns[0]) != inc_allele:
                    start = 0
                else:
                    start = 1
                i = 0
                for n in list(disFrame.columns[start:start+10]):
                    i += 1
                    handle.write("Neighbor # {}:\t\t\t\t{} ({})\n".format(str(i), n, str(disFrame.loc[inc_allele][n])))

                # infers sequence from nearest 10 neighbor vote
                for nemesis in missing_pos[inc_allele]:
                    new_val = self.seqvote(nemesis, inc_allele, disFrame, locDict, start)

                    locDict[inc_allele] = locDict[inc_allele][:nemesis] + \
                                new_val + locDict[inc_allele][nemesis+1:]
        return locDict

    def impute(self, loc, locDict, refseq, aaDict):
        seqIns = [m.start() for m in re.finditer('-', locDict[refseq])]
        aaIns = [m2.start() for m2 in re.finditer('-', str(aaDict[refseq].seq))]

        if os.path.isfile("{}data/locdist/{}.csv".format(pathloc,loc)):
            disFrame = pd.read_csv("{}data/locdist/{}.csv".format(pathloc,loc), index_col=0)
            if (set(disFrame.Index.tolist()) == set(locDict.keys())):
                next
            else:
                print("New alleles detected - recalculate distance matrix!")
                # TODO (gbiagini) - modify this to only add new alleles, not recalculate everything
                disFrame = self.distmatrix(locDict)
        else:
            disFrame = self.distmatrix(locDict)
            disFrame.to_csv(pathloc+"data/locdist/{}.csv".format(loc))

        missing_pos = {}
        missing_pos = {key: self.checkComplete(locDict[key], seqIns) for key in locDict.keys()}
        del(missing_pos)

        disFrame = disFrame.transpose()
        locDict = self.nearest10(loc, locDict, disFrame, missing_pos)
        locDict, incorrect = self.translate_nuc(locDict, seqIns)
        del(disFrame)
        translated = self.correction(incorrect, locDict, aaIns)
        return translated
        
    def main(self):
        aaDict = aa_mm.HLA_full_allele
        refseq = nuc_mm.refseq_full
        HLA_seq = nuc_mm.HLA_full_allele
        for loc in aa_mm.refseq:
            #if loc in ['A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4']:
            #    continue
            print("Processing locus " + loc + "...")
            newDict = {newKey: str(HLA_seq[newKey].seq) for newKey in HLA_seq.keys()}
            locDict = {locKey: newDict[locKey] for locKey in locDict.keys() if (
                locKey.split('*')[0] == loc)}
            del (newDict)
            imputed = self.impute(loc, locDict, refseq[loc], aaDict)
            # creates list from sequence strings for Pandas dataframe
            repDict = {repKey: list(imputed[repKey]) for repKey in imputed.keys()}
            del(imputed)
            repDict = self.finish_null(refseq[loc], repDict)
            impnew = {iKey: ''.join(repDict[iKey]) for iKey in repDict.keys()}
            my_records = []
            for pkey in impnew.keys():
                record = SeqRecord(
                    Seq(impnew[pkey]),
                    id=pkey,
                )
                my_records.append(record)
            SeqIO.write(my_records, "{}ifasta/{}_{}.fasta".format(pathloc,loc,nuc_mm.dbversion), "fasta")


Impute()
