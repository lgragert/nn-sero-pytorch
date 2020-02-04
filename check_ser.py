# used to check for values in the training, testing, and
# validation sets

from parse import _parse
from parse import _file_handler

def _check(tng_dict, tst_dict, val_dict, tdict):
    before = len(tdict)
    wsero = {}
    tst_to_tng = {}
    
    for key in list(tdict):
        if (tdict[key] != '') & (tdict[key] != '?'):
            wsero[key] = tdict[key]
        else:
            continue

    for key in tng_dict:
        for allele in list(tdict):
            if key == allele:
                del tdict[allele]
        for wallele in list(wsero):
            if key == wallele:
                del wsero[wallele]
    
    for key in tst_dict:
        for allele in list(tdict):
            if key == allele:
                del tdict[allele]
        for wallele in list(wsero):
            if key == wallele:
                tst_to_tng = wsero[key]

    for key in val_dict:
        for allele in list(tdict):
            if key == allele:
                del tdict[allele]

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

    for key in tst_to_tng:
        print(key)


    return(tdict)

A, B, C, DPB1, DQB1, DRB1 = _file_handler()
A_tng_dict, A_tst_dict, A_val_dict = A[0], A[1], A[2]
B_tng_dict, B_tst_dict, B_val_dict = B[0], B[1], B[2]
C_tng_dict, C_tst_dict, C_val_dict = C[0], C[1], C[2]
DPB1_tng_dict, DPB1_tst_dict, DPB1_val_dict = DPB1[0], DPB1[1], DPB1[2]
DQB1_tng_dict, DQB1_tst_dict, DQB1_val_dict = DQB1[0], DQB1[1], DQB1[2]
DRB1_tng_dict, DRB1_tst_dict, DRB1_val_dict = DRB1[0], DRB1[1], DRB1[2]

ser_file = open("rel_dna_ser_3370.txt", 'r')
A_dict = {}
B_dict = {}
C_dict = {}
DPB1_dict = {}
DQB1_dict = {}
DRB1_dict = {}

for line in ser_file:
    if line.find('#') == -1:
        if line.find('A*') != -1:
            allele = line.split(';')
            fixed = allele[1][:5]
            if fixed not in list(A_dict):
                A_dict[fixed] = allele[2]
        elif line.find('B*') != -1:
            allele = line.split(';')
            fixed = allele[1][:5]
            if fixed not in list(B_dict):
                B_dict[fixed] = allele[2]
        elif line.find('C*') != -1:
            allele = line.split(';')
            fixed = allele[1][:5]
            if fixed not in list(C_dict):
                C_dict[fixed] = allele[2]
        elif line.find('DPB1') != -1:
            allele = line.split(';')
            fixed = allele[1][:5]
            if fixed not in list(DPB1_dict):
                DPB1_dict[fixed] = allele[2]
        elif line.find('DQB1') != -1:
            allele = line.split(';')
            fixed = allele[1][:5]
            if fixed not in list(DQB1_dict):
                DQB1_dict[fixed] = allele[2]
        elif line.find('DRB1') != -1:
            allele = line.split(';')
            fixed = allele[1][:5]
            if fixed not in list(DRB1_dict):
                DRB1_dict[fixed] = allele[2]

print("Checking A...")
A_dict = _check(A_tng_dict, A_tst_dict, A_val_dict, A_dict)
print("Checking B...")
B_dict = _check(B_tng_dict, B_tst_dict, B_val_dict, B_dict)
print("Checking C...")
C_dict = _check(C_tng_dict, C_tst_dict, C_val_dict, C_dict)
print("Checking DPB1...")
DPB1_dict = _check(DPB1_tng_dict, DPB1_tst_dict, DPB1_val_dict, DPB1_dict)
print("Checking DQB1...")
DQB1_dict = _check(DQB1_tng_dict, DQB1_tst_dict, DQB1_val_dict, DQB1_dict)
print("Checking DRB1...")
DRB1_dict = _check(DRB1_tng_dict, DRB1_tst_dict, DRB1_val_dict, DRB1_dict)
print("Done.")