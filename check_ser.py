# used to check for values in the training, testing, and
# validation sets

from parse import *
import pandas as pd

def _check(tng_dict, tst_dict, val_dict, tdict):
    before = len(tdict)
    wsero = {}
    tst_to_tng = {}
    
    for key in list(tdict):
        if (tdict[key] != '') & (tdict[key] != '?'):
            wsero[key] = tdict[key]
        else:
            continue

    for key in tng_dict['allele']:
        for allele in list(tdict):
            if key == allele:
                del tdict[allele]
        for wallele in list(wsero):
            if key == wallele:
                del wsero[wallele]
    
    for key in val_dict['allele']:
        for allele in list(tdict):
            if key == allele:
                del tdict[allele]
        for wallele in list(wsero):
            if key == wallele:
                del wsero[wallele]
    
    for key in tst_dict['allele']:
        for allele in list(tdict):
            if key == allele:
                del tdict[allele]                
        for wallele in list(wsero):
            if key == wallele:
                tst_to_tng[key] = wsero[key]


    after = len(tdict)
    same = before - after

    for key in list(tdict):
        if (tdict[key] == '') or (tdict[key] == '?'):
            del tdict[key]
    
    n_sero = len(tdict)

    print('# Different Alleles: ' + str(after))
    print('# Same Alleles: ' + str(same))
    print('# Different with Serologies: ' + str(n_sero))
    print('# Test Set with Serologies: ' + str(len(tst_to_tng)))

    #for value in tst_to_tng:
        #print(value)


    return(tdict)

ser_file = open("rel_dna_ser_3370.txt", 'r')
A_dict = {}
B_dict = {}
C_dict = {}
DPB1_dict = {}
DQB1_dict = {}
DRB1_dict = {}

for line in ser_file:
    if line.find('#') == -1:
        if (line.find('A*') != -1) & (line.find('MIC') == -1):
            letter = line.split(';')
            nums = letter[1].split(':')
            allele = nums[0] + ':' + nums[1]
            fixed = letter[0] + allele
            if fixed not in list(A_dict):
                A_dict[fixed] = letter[2]
        elif (line.find('B*') != -1) & (line.find('MIC') == -1):
            letter = line.split(';')
            nums = letter[1].split(':')
            allele = nums[0] + ':' + nums[1]
            fixed = letter[0] + allele
            if fixed not in list(B_dict):
                B_dict[fixed] = letter[2]
        elif line.find('C*') != -1:
            letter = line.split(';')
            nums = letter[1].split(':')
            allele = nums[0] + ':' + nums[1]
            fixed = letter[0] + allele
            if fixed not in list(C_dict):
                C_dict[fixed] = letter[2]
        elif line.find('DPB1') != -1:
            letter = line.split(';')
            nums = letter[1].split(':')
            allele = nums[0] + ':' + nums[1]
            fixed = letter[0] + allele
            if fixed not in list(DPB1_dict):
                DPB1_dict[fixed] = letter[2]
        elif line.find('DQB1') != -1:
            letter = line.split(';')
            nums = letter[1].split(':')
            allele = nums[0] + ':' + nums[1]
            fixed = letter[0] + allele
            if fixed not in list(DQB1_dict):
                DQB1_dict[fixed] = letter[2]
        elif line.find('DRB1') != -1:
            letter = line.split(';')
            nums = letter[1].split(':')
            allele = nums[0] + ':' + nums[1]
            fixed = letter[0] + allele
            if fixed not in list(DRB1_dict):
                DRB1_dict[fixed] = letter[2]


loci = ['A', 'B', 'C', 'DPB1', 'DQB1', 'DRB1']

for one in loci:
    tng_frame = pd.read_csv(one + "_train.csv")
    tst_frame = pd.read_csv(one + "_test.csv")
    val_frame = pd.read_csv(one + "_validation.csv")

    tng_dict = tng_frame.to_dict(orient='list')
    tst_dict = tst_frame.to_dict(orient='list')
    val_dict = val_frame.to_dict(orient='list')
    if one == 'A':
        print("Checking A...")
        A_dict = _check(tng_dict, tst_dict, val_dict, A_dict)
    elif one == 'B':
        print("Checking B...")
        B_dict = _check(tng_dict, tst_dict, val_dict, B_dict)
    elif one == 'C':
        print("Checking C...")
        C_dict = _check(tng_dict, tst_dict, val_dict, C_dict)
    elif one == 'DPB1':
        print("Checking DPB1...")
        DPB1_dict = _check(tng_dict, tst_dict, val_dict, DPB1_dict)
    elif one == 'DQB1':
        print("Checking DQB1...")
        DQB1_dict = _check(tng_dict, tst_dict, val_dict, DQB1_dict)
    elif one == 'DRB1':
        print("Checking DRB1...")
        DRB1_dict = _check(tng_dict, tst_dict, val_dict, DRB1_dict)
        print("Done.")
