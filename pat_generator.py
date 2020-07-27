## Generation of updated pattern files for use in the RSNNS implementation of NN-Serology

import pandas as pd
import numpy as np
from datetime import date

loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]
today = date.today()
date = today.strftime("%B %d, %Y")

# training files 
for locus in loci:

    loc_frame = pd.read_csv("./RSNNS_fixed/training/" + locus + "_train.csv")
    
    pat_file = open("./pat_new/" + locus + ".tng.pat", "w+")
    pat_file.write("SNNS pattern definition file V3.2\n")
    pat_file.write("generated at " + date + "\n\n\n")
    pat_file.write("No. of patterns: " + str(len(loc_frame)) + "\n")
    pat_file.write("No. of input units: " + str(len(loc_frame.columns) - 2) + "\n")
    uni_list = loc_frame.serology.unique()

    new_list = []
    for uni in uni_list:
        new = uni.split(';')
        [new_list.append(x) for x in new if x not in new_list]
    
    uni_list = new_list
    for i in range(0,len(uni_list)):
        uni_list[i] = uni_list[i].strip('a')

    pat_file.write("No. of output units: " + str(len(uni_list)) + "\n")

    AAs_list = loc_frame.columns.tolist()
    AAs_list.remove('allele')
    AAs_list.remove('serology')
    polyAAs = ' '.join(AAs_list)
    polyAAs = "# " + polyAAs
    pat_file.write(polyAAs + "\n")
    
    for allele in loc_frame.values.tolist():
        serologies = {}
        for j in uni_list:
            serologies[j] = 0.00
        allele.pop(0)
        name = allele.pop(0)
        name = name.split('*')[1]
        sero = allele.pop()
        sero = sero.split(';')
        for i in range(0,len(sero)):
            sero[i] = sero[i].strip('a')
        for value in serologies.keys():
            if value in sero:
                serologies[value] = 1.00
        serval = list(serologies.values())
        in_serkey = ' '.join(list(serologies.keys()))
        in_allele = ' '.join(map(str,allele))
        in_serval = ' '.join(map(str,serval))

        pat_file.write("# " + name + "\n")
        pat_file.write("# input" + "\n")
        pat_file.write(in_allele + "\n")
        pat_file.write('# output ' + in_serkey + "\n")
        pat_file.write(in_serval + "\n")

    val_frame = pd.read_csv("./RSNNS_fixed/training/" + locus + "_validation.csv")
    
    val_file = open("./pat_new/" + locus + ".val.pat", "w+")
    val_file.write("SNNS pattern definition file V3.2\n")
    val_file.write("generated at " + date + "\n\n\n")
    val_file.write("No. of patterns: " + str(len(val_frame)) + "\n")
    val_file.write("No. of input units: " + str(len(val_frame.columns) - 2) + "\n")
    val_file.write("No. of output units: " + str(len(uni_list)) + "\n")
    
    for val_allele in val_frame.values.tolist():
        serologies = {}
        for j in uni_list:
            serologies[j] = 0.00
        val_allele.pop(0)
        val_name = val_allele.pop(0)
        val_name = val_name.split('*')[1]
        val_sero = val_allele.pop()
        val_sero = val_sero.split(';')
        for i in range(0,len(val_sero)):
            val_sero[i] = val_sero[i].strip('a')
        for value in serologies.keys():
            if value in val_sero:
                serologies[value] = 1.00
        val_serval = list(serologies.values())
        val_in_serkey = ' '.join(list(serologies.keys()))
        val_in_allele = ' '.join(map(str,val_allele))
        val_in_serval = ' '.join(map(str,val_serval))

        val_file.write("# " + val_name + "\n")
        val_file.write("# input" + "\n")
        val_file.write(val_in_allele + "\n")
        val_file.write('# output ' + val_in_serkey + "\n")
        val_file.write(val_in_serval + "\n")


    tst_frame = pd.read_csv("./RSNNS_fixed/testing/" + locus + "_test.csv")
    
    tst_file = open("./pat_new/" + locus + ".tst.pat", "w+")
    tst_file.write("SNNS pattern definition file V3.2\n")
    tst_file.write("generated at " + date + "\n\n\n")
    tst_file.write("No. of patterns: " + str(len(tst_frame)) + "\n")
    tst_file.write("No. of input units: " + str(len(tst_frame.columns) - 2) + "\n")
    tst_file.write("No. of output units: " + str(len(uni_list)) + "\n")
    
    for tst_allele in tst_frame.values.tolist():
        tst_allele.pop(0)
        tst_name = tst_allele.pop(0)
        tst_name = tst_name.split('*')[1]
        tst_allele.pop()
        tst_in_allele = ' '.join(map(str,tst_allele))
        tst_serval = []
        for each in in_serkey:
            tst_serval.append(0.00)
        tst_serval = ' '.join(map(str,tst_serval))

        tst_file.write("# testing " + tst_name + "\n")
        tst_file.write("# input" + "\n")
        tst_file.write(tst_in_allele + "\n")
        tst_file.write('# output ' + in_serkey + "\n")
        tst_file.write(tst_serval + "\n")


