import pandas as pd
import numpy as np

loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]

for locus in loci:
    summary = {}
    trn_set = pd.read_csv('training/' + locus + '_train.csv')
    tst_set = pd.read_csv('testing/' + locus + '_test.csv')

    old_trn_set = pd.read_csv('old_sets/train/' + locus + '_train.csv')
    old_val_set = pd.read_csv('old_sets/train/' + locus + '_validation.csv')
    old_tst_set = pd.read_csv('old_sets/test/' + locus + '_test.csv')

    summary['# training/validation alleles'] = str(len(trn_set))
    summary['# testing alleles'] = str(len(tst_set))
    summary['# previous training alleles'] = str(len(old_trn_set))
    summary['# previous validation alleles'] = str(len(old_val_set))
    summary['# previous testing alleles'] = str(len(old_tst_set))

    sum_df = pd.DataFrame(summary.items(), columns=['data', 'value'])
    sum_df.to_csv('summary/' + locus + '_summary.csv', index=False)

