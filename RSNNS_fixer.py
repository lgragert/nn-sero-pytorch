import pandas as pd
import numpy as np

loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]
for locus in loci:
    old_trn = pd.read_csv("./old_sets/train/" + locus + "_train.csv", low_memory=False)
    old_val = pd.read_csv("./old_sets/train/" + locus + "_validation.csv", low_memory=False)
    old_tst = pd.read_csv("./old_sets/test/" + locus + "_test.csv", low_memory=False)
    
    new_trn = pd.read_csv("./training/" + locus + "_train.csv", low_memory=False)
    new_val = pd.read_csv("./training/" + locus + "_validation.csv", low_memory=False)
    new_tst = pd.read_csv("./testing/" + locus + "_test.csv", low_memory=False)

    trn_droplist = []
    for column in new_trn.columns:
        if column not in old_trn.columns:
            trn_droplist.append(column)

    print(trn_droplist)
    fixed_trn = new_trn.drop(trn_droplist, axis=1)
    fixed_trn.to_csv('./RSNNS_fixed/training/' + locus + '_train.csv')


    
    val_droplist = []
    for column in new_val.columns:
        if column not in old_val.columns:
            val_droplist.append(column)
    
    fixed_val = new_val.drop(val_droplist, axis=1)
    fixed_val.to_csv('./RSNNS_fixed/training/' + locus + '_validation.csv')

    tst_droplist = []
    for column in new_tst.columns:
        if column not in old_tst.columns:
            tst_droplist.append(column)
        else:
            next
    
    fixed_tst = new_tst.drop(tst_droplist, axis=1)
    fixed_tst.to_csv('./RSNNS_fixed/testing/' + locus + '_test.csv')
    