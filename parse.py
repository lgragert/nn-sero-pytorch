## Name: Giovanni Biagini
## PI: Loren Gragert, PhD
## Institution: Tulane University School of Medicine
## Department: Pathology and Laboratory Medicine
## Location: Louisiana Cancer Research Center
## Date Begun: 01/16/2020
## Purpose: To modernize SNNS

# figure out how to parse the training, testing, and validation files into vectors.
import torch
import numpy as np

# multi-purpose function to parse each type of file
# it might behoove us to attempt an OOP approach, or to incorporate regular expressions to simplify the parsing
def _parse(tng_file, tst_file, val_file):
    tng_dict = {}
    tst_dict = {}
    val_dict = {}
    AA_dict = {}

    # loop through each line of the training file
    for each in tng_file:
        # will only be executed once, to create a list of the polymorphic amino acids
        if each.find("No. of output") != -1:
            each = next(tng_file)
            AA = each.strip("# ")
            AA_list = AA.split()
            continue
        # main purpose, identifying the hashtags at the beginning of the significant lines
        if each.find('#') != -1:
            # gathering information on a specific allele
            if (each.find('# input') == -1) & (each.find('# output') == -1):
                temp = each.strip('# ')
                allele = temp.rstrip()
                continue
            # looking for the binary values corresponding to the amino acid polymorphisms
            elif each.find('# input') != -1:
                line = next(tng_file)
                # generating a list of the binary values
                bin_val = line.split()
                # the values for the specific amino acid are zipped into a dictionary with the polymorphic AAs as the keys
                AA_dict = dict(zip(AA_list, bin_val))
                # the dictionary of AA polymorphisms is added to the (nested) dictionary of all alleles, with its specific allele as its key
                tng_dict[allele] = AA_dict

    # basically the same, but to parse the testing file
    for each in tst_file:
        if each.find('#') != -1:
            if (each.find('# input') == -1) & (each.find('# output') == -1):
                temp = each.strip('# testing ')
                allele = temp.rstrip()
                continue
            elif each.find('# input') != -1:
                line = next(tst_file)
                bin_val = line.split()
                AA_dict = dict(zip(AA_list, bin_val))
                tst_dict[allele] = AA_dict

     # final loop to parse the validation file       
    for each in val_file:
        if each.find('#') != -1:
            if (each.find('# input') == -1) & (each.find('# output') == -1):
                temp = each.strip('# ')
                allele = temp.rstrip()
                continue
            elif each.find('# input') != -1:
                line = next(val_file)
                bin_val = line.split()
                AA_dict = dict(zip(AA_list, bin_val))
                val_dict[allele] = AA_dict
    return(tng_dict, tst_dict, val_dict)

def _file_handler()
    # opening files to send to the parser
    # there is almost definitely a much simpler way to code this
    A_tng_file = open("A.tng.pat", 'r')
    A_tst_file = open("A.tst.pat", 'r')
    A_val_file = open("A.val.pat", 'r')
    B_tng_file = open("B.tng.pat", 'r')
    B_tst_file = open("B.tst.pat", 'r')
    B_val_file = open("B.val.pat", 'r')
    C_tng_file = open("C.tng.pat", 'r')
    C_tst_file = open("C.tst.pat", 'r')
    C_val_file = open("C.val.pat", 'r')
    DPB1_tng_file = open("DPB1.tng.pat", 'r')
    DPB1_tst_file = open("DPB1.tst.pat", 'r')
    DPB1_val_file = open("DPB1.val.pat", 'r')
    DQB1_tng_file = open("DQB1.tng.pat", 'r')
    DQB1_tst_file = open("DQB1.tst.pat", 'r')
    DQB1_val_file = open("DQB1.val.pat", 'r')
    DRB1_tng_file = open("DRB1.tng.pat", 'r')
    DRB1_tst_file = open("DRB1.tst.pat", 'r')
    DRB1_val_file = open("DRB1.val.pat", 'r')

    # function call to the parser
    # this part is actually pretty succint, but it could be shorter
    A_tng_dict, A_tst_dict, A_val_dict = _parse(A_tng_file, A_tst_file, A_val_file)
    B_tng_dict, B_tst_dict, B_val_dict = _parse(B_tng_file, B_tst_file, B_val_file)
    C_tng_dict, C_tst_dict, C_val_dict = _parse(C_tng_file, C_tst_file, C_val_file)
    DPB1_tng_dict, DPB1_tst_dict, DPB1_val_dict = _parse(DPB1_tng_file, DPB1_tst_file, DPB1_val_file)
    DQB1_tng_dict, DQB1_tst_dict, DQB1_val_dict = _parse(DQB1_tng_file, DQB1_tst_file, DQB1_val_file)
    DRB1_tng_dict, DRB1_tst_dict, DRB1_val_dict = _parse(DRB1_tng_file, DRB1_tst_file, DRB1_val_file)
    
    # lists of dictionaries for easier return
    A = [A_tng_dict, A_tst_dict, A_val_dict]
    B = [B_tng_dict, B_tst_dict, B_val_dict]
    C = [C_tng_dict, C_tst_dict, C_val_dict]
    DPB1 = [DPB1_tng_dict, DPB1_tst_dict, DPB1_val_dict]
    DQB1 = [DQB1_tng_dict, DQB1_tst_dict, DQB1_val_dict]
    DRB1 = [DRB1_tng_dict, DRB1_tst_dict, DRB1_val_dict]

    # closing files after parsing
    # there is almost definitely a much simpler way to code this
    A_tng_file.close()
    A_tst_file.close()
    A_val_file.close()
    B_tng_file.close()
    B_tst_file.close()
    B_val_file.close()
    C_tng_file.close()
    C_tst_file.close()
    C_val_file.close()
    DPB1_tng_file.close()
    DPB1_tst_file.close()
    DPB1_val_file.close()
    DQB1_tng_file.close()
    DQB1_tst_file.close()
    DQB1_val_file.close()
    DRB1_tng_file.close()
    DRB1_tst_file.close()
    DRB1_val_file.close()

    # look into a potentially nested function that generates the locations for the input files and assigns to different variables in order
    # to cut down on the verbosity of the above code
    return(A, B, C, DPB1, DQB1, DRB1)
