#!/usr/bin/env python

###############################################################################
#   SCRIPT NAME:    nuc_matching_msf.py
#   DESCRIPTION:    Module for nucleotide matching functions
#   OUTPUT:
#   DATE:           September 01, 2020
#	LAST UPDATE:	May 1st, 2021
#	REASON:			Updates to incorporate into main library.
#   AUTHOR:         Giovanni Biagini (dbiagini@tulane.edu ; GitHub: gbiagini)
#   PI:             Loren Gragert, Ph.D.
#   ORGANIZATION:   Tulane University School of Medicine
#   NOTES:			
###############################################################################

import re
import os.path
from os import path
import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import requests

pathloc = str(os.getcwd()) + '/'

class NucMatch():

    def __init__(self, dbversion=3420, ungap=False):
        # if ungap = True, use complete reference sequence
        # if ungap = False, use gapped reference sequence
        self.dbversion = dbversion
        self.ungap = ungap
        hlaProteinOffset = {
            "A" : 0,
            "B" : 0,
            "C" : 0, 
            "DRA" : 0,
            "DRB1" : 0,
            "DRB3" : 0,
            "DRB4" : 0,
            "DRB5" : 0,
            "DQA1" : 0,
            "DQB1" : 0,
            "DPA1" : 0,
            "DPB1" : 0,
        }

        first_ten = {
            "A*01:01:01:01": "GGCTCCCACTCCATGAGGTATTTCTTCACA",
            "B*07:02:01:01": "GGCTCCCACTCCATGAGGTATTTCTACACC",
            "C*01:02:01:01": "TGCTCCCACTCCATGAAGTATTTCTTCACA",
            "DRB1*01:01:01:01": "GGGGACACCCGACCACGTTTCTTGTGGCAG",
            "DRB3*01:01:02:01": "GGGGACACCCGACCACGTTTCTTGGAGCTG",
            "DRB4*01:01:01:01": "GGGGACACCCAACCACGTTTCTTGGAGCAG",
            "DRB5*01:01:01:01": "GGGGACACCCGACCACGTTTCTTGCAGCAG",
            "DQA1*01:01:01:01": "GAAGACATTGTGGCTGACCACGTTGCCTCT",
            "DQB1*05:01:01:01": "AGAGACTCTCCCGAGGATTTCGTGTACCAG",
            "DPA1*01:03:01:01": "ATCAAGGCGGACCATGTGTCAACTTATGCC",
            "DPB1*01:01:01:01": "AGGGCCACTCCAGAGAATTACGTGTACCAG"
        }
        refseq_full = {
           "A" : "A*01:01:01:01",
            "B" : "B*07:02:01:01",
            "C" : "C*01:02:01:01",
            "DRB1" : "DRB1*01:01:01:01",
            "DRB3" : "DRB3*01:01:02:01",
            "DRB4" : "DRB4*01:01:01:01",
            "DRB5" : "DRB5*01:01:01:01",
            "DQA1" : "DQA1*01:01:01:01",
            "DQB1" : "DQB1*05:01:01:01",
            "DPA1" : "DPA1*01:03:01:01",
            "DPB1" : "DPB1*01:01:01:01",
        }

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

            # use only ARD positions
        ard_start_pos = {
            "A" : 1,
            "B" : 1,
            "C" : 1,
            "DRB1" : 1,
            "DRB3" : 1,
            "DRB4" : 1,
            "DRB5" : 1,	
            "DQA1" : 1,
            "DQB1" : 1,
            "DPA1" : 1,
            "DPB1" : 1,
        }
        ard_end_pos = {
            "A" : 182,
            "B" : 182,
            "C" : 182,
            "DRB1" : 94,
            "DRB3" : 94,
            "DRB4" : 94,
            "DRB5" : 94,	
            "DQA1" : 94,
            "DQB1" : 95, #increased by 1
            "DPA1" : 94,
            "DPB1" : 94,
        }

        # lots of incomplete ARD sequences in IMGT/HLA
        # this should be handled in reference alignment, not here
        ard_start_pos_incomplete = {
            "A" : 2, # A*02:50
            "B" : 2, # B*07:30
            "C" : 2, # C*01:10
            "DRB1" : 7, # DRB1*08:19
            "DRB3" : 2, # TBD
            "DRB4" : 2, # TBD
            "DRB5" : 2, # TBD
            "DQA1" : 6, # DQA1*01:06
            "DQB1" : 6, # DQB1*05:100
            "DPA1" : 11, # DPA1*01:03:02
            "DPB1" : 6, # DPB1*01:01:03
        }
        ard_end_pos_incomplete = {
            "A" : 182, # A*02:50
            "B" : 182, # B*07:30
            "C" : 182, # C*01:10
            "DRB1" : 92, # DRB1*08:19
            "DRB3" : 2, # TBD
            "DRB4" : 2, # TBD
            "DRB5" : 2, # TBD
            "DQA1" : 87, # DQA1*01:06
            "DQB1" : 94, # DQB1*05:100
            "DPA1" : 84, # DPA1*01:03:02
            "DPB1" : 92, # DPB1*01:01:03
        }

        self.refseq_full = refseq_full
        self.refseq = refseq
        self.ard_start_pos = ard_start_pos
        self.ard_end_pos = ard_end_pos
        self.ard_start_pos_incomplete = ard_start_pos_incomplete
        self.ard_end_pos_incomplete = ard_end_pos_incomplete
        self.hlaProteinOffset = hlaProteinOffset
        self.first_ten = first_ten
        self.main()

        if self.ungap == True:
            self.remove_gap()

    # generate a regular expression to find first ten residues of each locus
    # ignoring '-' characters.
    def regex_gen(self):
        regexes = {}
        for each in self.first_ten.keys():
            regex = ""
            for i in range(0, len(self.first_ten[each])):
                regex += str(self.first_ten[each][i])+"[^"+self.first_ten[each]+"]*?"
            regexes[each] = regex
        return regexes

    # Align the coordinate system as appropriate to the beginning of the mature
    # protein sequence.
    def coordinate(self, regex, sequence):
        o = re.search(regex, str(sequence))
        offset = o.start()
        return offset

    def getMatureProteinOffset(self, locus):
            return self.hlaProteinOffset.get(locus, "Invalid HLA Locus")
        
    def adjust_end(self, multipleseq, loc, ard_start_pos, ard_start_pos_incomplete, ard_end_pos,
                    ard_end_pos_incomplete, prev=0, prev_inc=0):
        loc_full_allele = self.refseq_full[loc]
        full_protein = multipleseq[id==loc_full_allele].seq
        mature_protein = full_protein[self.getMatureProteinOffset(loc):]
        start = ard_start_pos
        end = ard_end_pos
        start_inc = ard_start_pos_incomplete
        end_inc = ard_end_pos_incomplete
        count = mature_protein[start:end].count('-')
        count_inc = mature_protein[start_inc:end_inc].count('-')

        check = count - prev
        check_inc = count_inc - prev_inc

        if (check == 0): 
            new_end = end + (count-prev)
            new_end_inc = end_inc + (count_inc-prev_inc)
            newlist = [new_end, new_end_inc]
            return newlist
        else:
            new_end = end + (count - prev)
            prev = count
            if (check_inc != 0):
                new_end_inc = end_inc + (count_inc - prev_inc)
                prev_inc = count_inc
                return self.adjust_end(multipleseq, loc, start, start_inc, new_end,
                                    new_end_inc, prev, prev_inc)
            else:
                new_end_inc = end_inc
                return self.adjust_end(multipleseq, loc, start, start_inc, new_end,
                                    new_end_inc, prev, prev_inc=0)

    def remove_ins(self, loc_full_alseq):
        # need to remove the inserts from the reference sequence to print into
        # the IMGT/HLA .txt file
        gapframe = pd.DataFrame.from_dict(loc_full_alseq, orient="index")
        droplist = []
        for i, row in gapframe.iterrows():
            if i == self.refseq[loc]:
                for name, data in gapframe.iteritems():
                    if data[i] == '-':
                        droplist.append(name)
                    else:
                        continue
            else:
                continue
        ungapframe = gapframe.drop(droplist, axis=1)
        for j, jrow in ungapframe.iterrows():
            jrow = jrow.to_string(header=False, index=False)
            jrow = jrow.replace('\n', '')
            jrow = jrow.replace(' ', '')
            loc_full_alseq[j] = jrow
        
        return loc_full_alseq

    def reference(self, multipleseq, loc):
        loc_full_allele = self.refseq_full[loc]
        regexes = self.regex_gen()
        mslist = [x for x in multipleseq if x.id == loc_full_allele]
        loc_full_protein = mslist[0]
        #loc_full_protein = multipleseq[id==loc_full_allele]
        full_protein = loc_full_protein.seq
        
        # align coordinate system based on reference sequence.
        # uses two field allele if full allele not present.
        # does not overwrite full allele if present.
        # necessary for earlier database versions.
        offset = self.coordinate(regexes[self.refseq_full[loc]], full_protein)
        self.hlaProteinOffset[loc] = offset

        return

    def adjust(self, multipleseq, loc):
        # then adjust end position based on dashes.
        new_end, new_end_inc = self.adjust_end(multipleseq, loc, self.ard_start_pos[loc],
                                            self.ard_start_pos_incomplete[loc],
                                            self.ard_end_pos[loc],
                                            self.ard_end_pos_incomplete[loc],
                                            prev=0, prev_inc=0)
        self.ard_end_pos[loc] = new_end
        self.ard_end_pos_incomplete[loc] = new_end_inc
        return


    def generate_IMGT(self, HLA_full_alseq):
        outfile = open("./IMGT_HLA_Full_Protein_" + str(self.dbversion) + ".txt", "w+")
        outfile.write("Allele\tFull_Protein\n")
        for allele_loctype in HLA_full_alseq:
            outfile.write("HLA-" + allele_loctype + "\t" +
                        str(HLA_full_alseq[allele_loctype]) + "\n")
        outfile.close()
        return

    def remove_gap(self):
        loci =  ['A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
        refgaps2 = {}
        refgaps4 = {}
        for locus in loci:
            seq2 = str(self.HLA_seq[self.refseq[locus]].seq)
            seq4 = str(self.HLA_full_allele[self.refseq_full[locus]].seq)
            refgaps2[locus] = [i for i, ltr in enumerate(seq2) if ltr == '-']
            refgaps4[locus] = [j for j, lttr in enumerate(seq4) if lttr == '-']
        for x2 in self.HLA_seq.keys():
            locus = x2.split('*')[0]
            gapseq2 = str(self.HLA_seq[x2].seq)
            seq2x = ''.join([gapseq2[i] for i in range(len(gapseq2)) if not(gapseq2[i] == '-' and (i in refgaps2[locus]))])
            self.HLA_seq[x2].seq = Seq(seq2x)
        for x4 in self.HLA_full_allele.keys():
            locus = x4.split('*')[0]
            gapseq4 = str(self.HLA_full_allele[x4].seq)
            seq4x = ''.join([gapseq4[i] for i in range(len(gapseq4)) if not(gapseq4[i] == '-' and (i in refgaps4[locus]))])
            self.HLA_full_allele[x4].seq = Seq(seq4x)
        return

    def main(self):
        loci =  ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
        self.HLA_full_allele = {} # Full four-field allele names
        #self.HLA_full_alseq = {} # Only used for generating IMGTHLA files
        self.HLA_seq = {} # Two-field
        regex = '\w*\*\d*\:\d*'
        suffixes = ["L", "S", "C", "A", "Q", "N"]

        for locus in loci:
            loc_full_alseq = {}
            seq_filename = "{}nuc_msf/{}_nuc_{}.msf".format(pathloc,locus,self.dbversion)
            if not os.path.isfile(seq_filename):
                print("Downloading requested MSF nucleotide files for locus " + locus + "...")
                url = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/{}/msf/{}_nuc.msf".format(self.dbversion, locus)
                r = requests.get(url)
                with open(seq_filename, 'w') as f:
                    f.write(r.content)
            else: 
                print("MSF nucleotide files already downloaded")
            multipleseq = AlignIO.read(seq_filename, format="msf")
            
            if locus == 'DRB345':
                for loc in ['DRB3', 'DRB4', 'DRB5']:
                    self.reference(multipleseq, loc)
                    self.adjust(multipleseq, loc)
            else: 
                self.reference(multipleseq, locus)
                self.adjust(multipleseq, locus)

            for record in multipleseq:
                loc_full_allele = record.id
                (loc,full_allele) = loc_full_allele.split("*")
                # append the suffix - needed for null alleles
                # index starts at 1 since some loci characters are also suffixes
                separator = loc_full_allele.find("*")
                if any(x in loc_full_allele[separator:] for x in suffixes):
                    loc_two_field_allele = re.match(regex, loc_full_allele).group() + \
                                        loc_full_allele[-1]
                else:
                    loc_two_field_allele = re.match(regex, loc_full_allele).group()
                full_protein = record.seq
                nogap = full_protein

                loc_full_alseq[loc_full_allele] = SeqRecord(nogap)

                # skip missing sequences in IMGT/HLA file
                if (len(full_protein) <10): # #NAME?
                    print ("Missing Sequence:" + loc_full_allele)
                    continue

                mature_protein = full_protein[self.getMatureProteinOffset(loc):]
                mrecord = SeqRecord(mature_protein)


                # full allele name
                self.HLA_full_allele[loc_full_allele] = mrecord

                # don't overwrite two-field alleles with new sequences - more likely
                # to be incomplete
                if (loc_two_field_allele not in self.HLA_seq):
                    self.HLA_seq[loc_two_field_allele] = mrecord

                # print (HLA_seqrecord_dict[allele])
                
                # TODO - add feature annotation to SeqIO object
                # https://biopython.org/wiki/SeqRecord
                # e.g. - which allele contain Bw4/Bw6 epitopes, bind LILRB1, etc
                # https://www.biostars.org/p/57549/

                # print (HLA_seqrecord_dict[allele].seq)
            
            #!GB!# Commented out these two lines as well as the generate_IMGT()
            # function call, since the process is intensive and does not need to be
            #!GB!# run every time the code is used. Might should consider a switch to
            # splitting a HLA_full_alseq DataFrame by locus and doing all of the work in
            #!GB!# a single function definition.
            #loc_full_alseq = remove_ins(loc_full_alseq)
            #HLA_full_alseq.update(loc_full_alseq)

        #generate_IMGT(HLA_full_alseq)
