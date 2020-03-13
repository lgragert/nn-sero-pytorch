import pandas as pd
import numpy as np

loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]
summary = {}

for locus in loci:
    summary_data = {}
    trn_set = pd.read_csv('training/' + locus + '_train.csv')
    val_set = pd.read_csv('training/' + locus + '_validation.csv')
    tst_set = pd.read_csv('testing/' + locus + '_test.csv')

    old_trn_set = pd.read_csv('old_sets/train/' + locus + '_train.csv')
    old_val_set = pd.read_csv('old_sets/train/' + locus + '_validation.csv')
    old_tst_set = pd.read_csv('old_sets/test/' + locus + '_test.csv')

    trnlen = float(len(trn_set))
    vallen = float(len(val_set))
    tstlen = float(len(tst_set))
    polyAA = float(len(trn_set.iloc[0])) - 1
    oldtrnlen = float(len(old_trn_set))
    oldvallen = float(len(old_val_set))
    oldtstlen = float(len(old_tst_set))
    oldpolyAA = float(len(old_trn_set.iloc[0])) - 1

    summary_data['Number of Training Alleles'] = trnlen
    summary_data['R-SNNS Number of Training Alleles'] = oldtrnlen
    summary_data['Difference in Training Set'] = trnlen - oldtrnlen
    summary_data['Percent (%) Growth in Training Set'] = ((trnlen - oldtrnlen)/oldtrnlen) * 100
    summary_data['Number of Validation Alleles'] = vallen
    summary_data['R-SNNS Number of Validation Alleles'] = oldvallen
    summary_data['Difference in Validation Set'] = vallen - oldvallen
    summary_data['Percent (%) Growth in Validation Set'] = ((vallen - oldvallen)/oldvallen) * 100
    summary_data['Number of Testing Alleles'] = tstlen
    summary_data['R-SNNS Number of Testing Alleles'] = oldtstlen
    summary_data['Difference in Testing Set'] = tstlen - oldtstlen
    summary_data['Percent (%) Growth in Testing Set'] = ((tstlen - oldtstlen)/oldtstlen) * 100
    summary_data['Number of Polymorphisms'] = polyAA
    summary_data['R-SNNS Number of Polymorphisms'] = oldpolyAA
    summary_data['Difference in Polymorphisms'] = polyAA - oldpolyAA
    summary_data['Percent (%) Growth in Polymorphisms'] = ((polyAA - oldpolyAA)/oldpolyAA) * 100
    summary[locus] = summary_data
    
sum_df = pd.DataFrame(data=summary)

sum_df.to_csv('summary/summary.csv', index=True)

