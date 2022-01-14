#!/usr/bin/env python

###############################################################################
#   SCRIPT NAME:    nuc_matching_msf.py
#   DESCRIPTION:    Module for nucleotide matching functions
#   OUTPUT:
#   DATE:           September 01, 2020
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
import random
import requests

# TODO (gbiagini) - wrap outer code in a main() function to avoid running
#  every time scripts are imported

hlaProteinOffset = {
    "A" : 23, # 365 versus 341 mature (decreased by 1))
    "B" : 24,
    "C" : 23, # decreased by 1
    "DRA" : 25,
    "DRB1" : 29,
    "DRB3" : 29,
    "DRB4" : 29,
    "DRB5" : 29,
    "DQA1" : 23,
    "DQB1" : 27, #(decreased by 5 to match RSNNS pat files)
    "DPA1" : 31,
    "DPB1" : 35, # increased by 6
}

def getMatureProteinOffset(locus):
	return hlaProteinOffset.get(locus, "Invalid HLA Locus")


def adjust_offset(loc, ard_start_pos, ard_start_pos_incomplete, ard_end_pos,
                  ard_end_pos_incomplete, prev=0, prev_inc=0):
	start = ard_start_pos
	end = ard_end_pos
	start_inc = ard_start_pos_incomplete
	end_inc = ard_end_pos_incomplete
	count = mature_protein[start:end].count('-')
	count_inc = mature_protein[start_inc:end_inc].count('-')

	check = count - prev
	check_inc = count_inc - prev_inc

	if (check == 0):
		new_end = end + (count - prev)
		new_end_inc = end_inc + (count_inc - prev_inc)
		newlist = [new_end, new_end_inc]
		return newlist
	else:
		new_end = end + (count - prev)
		prev = count
		if (check_inc != 0):
			new_end_inc = end_inc + (count_inc - prev_inc)
			prev_inc = count_inc
			return adjust_offset(loc, start, start_inc, new_end,
			                     new_end_inc, prev, prev_inc)
		else:
			new_end_inc = end_inc
			return adjust_offset(loc, start, start_inc, new_end,
			                     new_end_inc, prev, prev_inc=0)

# use only ARD positions
ard_start_pos = {
	"A": 1,
	"B": 1,
	"C": 1,
	"DRB1": 1,
	"DRB3": 1,
	"DRB4": 1,
	"DRB5": 1,
	"DQA1": 1,
	"DQB1": 1,
	"DPA1": 1,
	"DPB1": 1,
}
# FIXME (gbiagini) - Figure out why the shortened sequence error is only
#  happening with DQB1 and fix it.
ard_end_pos = {
	"A": 182,
	"B": 182,
	"C": 182,
	"DRB1": 94,
	"DRB3": 94,
	"DRB4": 94,
	"DRB5": 94,
	"DQA1": 94,
	"DQB1": 95, #increased by 1 because of weird issue stopping early only on
	# DQB1
	"DPA1": 94,
	"DPB1": 94,
}

# lots of incomplete ARD sequences in IMGT/HLA
# this should be handled in reference alignment, not here
ard_start_pos_incomplete = {
	"A": 2,  # A*02:50
	"B": 2,  # B*07:30
	"C": 2,  # C*01:10
	"DRB1": 7,  # DRB1*08:19
	"DRB3": 2,  # TBD
	"DRB4": 2,  # TBD
	"DRB5": 2,  # TBD
	"DQA1": 6,  # DQA1*01:06
	"DQB1": 6,  # DQB1*05:100
	"DPA1": 11,  # DPA1*01:03:02
	"DPB1": 6,  # DPB1*01:01:03
}
ard_end_pos_incomplete = {
	"A": 182,  # A*02:50
	"B": 182,  # B*07:30
	"C": 182,  # C*01:10
	"DRB1": 92,  # DRB1*08:19
	"DRB3": 2,  # TBD
	"DRB4": 2,  # TBD
	"DRB5": 2,  # TBD
	"DQA1": 87,  # DQA1*01:06
	"DQB1": 94,  # DQB1*05:100
	"DPA1": 84,  # DPA1*01:03:02
	"DPB1": 92,  # DPB1*01:01:03
}

# load protein sequence file to get full protein sequences into SeqRecord file
##seq_filename = "IMGT_HLA_Full_Protein_3330.txt"
##seqfile = open(seq_filename, "r")

loci = ['A', 'B', 'C', 'DRB1', 'DRB345', 'DQA1', 'DQB1', 'DPA1', 'DPB1']

refseq_full = {
	"A": "A*01:01:01:01",
	"B": "B*07:02:01:01",
	"C": "C*01:02:01:01",
	"DRB1": "DRB1*01:01:01:01",
	"DRB3": "DRB3*01:01:02:01",
	"DRB4": "DRB4*01:01:01:01",
	"DRB5": "DRB5*01:01:01:01",
	"DQA1": "DQA1*01:01:01:01",
	"DQB1": "DQB1*05:01:01:01",
	"DPA1": "DPA1*01:03:01:01",
	"DPB1": "DPB1*01:01:01:01",
}

refseq = {
	"A": "A*01:01",
	"B": "B*07:02",
	"C": "C*01:02",
	"DRB1": "DRB1*01:01",
	"DRB3": "DRB3*01:01",
	"DRB4": "DRB4*01:01",
	"DRB5": "DRB5*01:01",
	"DQA1": "DQA1*01:01",
	"DQB1": "DQB1*05:01",
	"DPA1": "DPA1*01:03",
	"DPB1": "DPB1*01:01",
}

HLA_full_allele = {}  # Full four-field allele names
HLA_full_alseq = {}
HLA_seq = {}  # Two-field
regex = '\w*\*\d*\:\d*'
suffixes = ["L", "S", "C", "A", "Q", "N"]

##placeholder to test removal of null alleles from imputation
#quest = input("Specify HLA/IMGT DB version? (y/n) ")
#if (quest == "y") or (quest == "Y"):
#	dbversion = input("HLA/IMGT DB version? (x.xx.x) ")
#	dbversion = dbversion.strip()
#	dbversion = dbversion.replace('.', '')
#else:
#	dbversion = "3400"

dbversion = "3400"

# The "$loc_nuc.msf" files are supposed to only contain coding sequences,
# so we can start with using those as a basis for translation,
# and translating the codons directly before making alterations.
for locus in loci:
	loc_full_alseq = {}
	seq_filename = "../aa-matching/nuc_msf/" + locus + "_nuc_" + str(
		dbversion) \
	               + ".msf"
	if path.exists(seq_filename) == False:
		print("Downloading requested MSF nucleotide files for locus " + locus +
		      "...")
		url = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/" + str(
			dbversion) + "/msf/" + locus + "_nuc.msf"
		r = requests.get(url)
		with open(seq_filename, 'wb+') as f:
			f.write(r.content)
	else:
		print("MSF nucleotide files already downloaded")
	multipleseq = AlignIO.read(seq_filename, format="msf")

	for record in multipleseq:
		loc_full_allele = record.id
		# append the suffix - needed for null alleles
		# index starts at 1 since some loci characters are also suffixes
		separator = loc_full_allele.find("*")
		if any(x in loc_full_allele[separator:] for x in suffixes):
			loc_two_field_allele = re.match(regex, loc_full_allele).group() +\
			                       loc_full_allele[-1]
		else:
			loc_two_field_allele = re.match(regex, loc_full_allele).group()
		full_protein = record.seq
		nogap = full_protein

		loc_full_alseq[loc_full_allele] = SeqRecord(nogap)

		# skip missing sequences in IMGT/HLA file
		if (len(full_protein) < 10):  # #NAME?
			print("Missing Sequence:" + loc_full_allele)
			continue

		(loc, full_allele) = loc_full_allele.split("*")

		# TODO (gbiagini) - Keep in mind that this will need to be
		#  reimplemented using the translated protein sequence
		# mature_protein = full_protein[getMatureProteinOffset(loc):]
		#
		mature_protein = full_protein
		if loc_full_allele == refseq_full[loc]:
		 	new_end, new_end_inc = adjust_offset(loc, ard_start_pos[loc],
		 	                                     ard_start_pos_incomplete[loc],
		 	                                     ard_end_pos[loc],
		 	                                     ard_end_pos_incomplete[loc],
		 	                                     prev=0, prev_inc=0)
		 	ard_end_pos[loc] = new_end
		 	ard_end_pos_incomplete[loc] = new_end_inc

		# print(mature_protein)

		mrecord = SeqRecord(mature_protein)

		# full allele name
		HLA_full_allele[loc_full_allele] = mrecord

		# don't overwrite two-field alleles with new sequences -  more likely
		# to be incomplete
		if (loc_two_field_allele not in HLA_seq):
			HLA_seq[loc_two_field_allele] = mrecord

		# print (HLA_seqrecord_dict[allele])

		# TODO - add feature annotation to SeqIO object
		# https://biopython.org/wiki/SeqRecord
		# e.g. - which allele contain Bw4/Bw6 epitopes, bind LILRB1, etc
		# https://www.biostars.org/p/57549/

		# print (HLA_seqrecord_dict[allele].seq)
