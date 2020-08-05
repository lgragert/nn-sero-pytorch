import pandas as pd
import numpy as np
from datetime import date
from collections import OrderedDict

def ser_parse():
    current = '3370'
    loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]
    ref_seq = ["A*01:01", "B*07:02", "C*01:02", "DPB1*01:01", "DQB1*05:01", "DRB1*01:01"]
    for locus in loci:
        indy = loci.index(locus)
        ser_handler = open('rel_dna_ser/rel_dna_ser_' + current + '.txt', 'r')
        loc_test = locus + '*'
        ser_dict = {}
        for line in ser_handler:
            serology = []
            if line.find('*') != -1:
                info = line.split(';')
                if info[0] == loc_test:
                    numsplit = info[1].split(':')
                    nums = numsplit[0] + ':' + numsplit[1]
                    allele = info[0] + nums
                    if info[2] == '?':
                        info[2] = np.nan
                    elif info[2] == '0':
                        info[2] = np.nan
                    elif info[2] == '':
                        info[2] = np.nan
                    else:
                        info[2] = info[2] + 'a'
                    serology.append(str(info[2]))
                    if locus == 'A':
                        if serology[0] in ['23a', '24a']:
                            serology.append('9a')
                        elif serology[0] in ['25a', '26a', '34a', '66a']:
                            serology.append('10a')
                        elif serology[0] in ['29a', '30a', '31a', '32a', '33a', '74a']:
                            serology.append('19a')
                        elif serology[0] in ['68a', '69a']:
                            serology.append('28a')
                    elif locus == 'B':
                        if serology[0] in ['51a', '52a']:
                            serology.append('5a')
                        elif serology[0] in ['44a', '45a']:
                            serology.append('12a')
                        elif serology[0] in ['64a', '65a']:
                            serology.append('14a')
                        elif serology[0] in ['62a', '63a', '75a', '76a', '77a']:
                            serology.append('15a')
                        elif serology[0] in ['38a', '39a']:
                            serology.append('16a')
                        elif serology[0] in ['57a', '58a']:
                            serology.append('17a')
                        elif serology[0] in ['49a', '50a']:
                            serology.append('21a')
                        elif serology[0] in ['54a', '55a', '56a']:
                            serology.append('22a')
                        elif serology[0] in ['60a', '61a']:
                            serology.append('40a')
                        elif serology[0] in ['71a', '72a']:
                            serology.append('70a')
                    elif locus == 'C':
                        if serology[0] in ['9a', '10a']:
                            serology.append('3a')
                    elif locus == 'DQB1':
                        if serology[0] in ['5a', '6a']:
                            serology.append('1a')
                        elif serology[0] in ['7a', '8a', '9a']:
                            serology.append('3a')
                    elif locus == 'DRB1':
                        if serology[0] in ['15a', '16a']:
                            serology.append('2a')
                        elif serology[0] in ['17a', '18a']:
                            serology.append('3a')
                        elif serology[0] in ['11a', '12a']:
                            serology.append('5a')
                        elif serology[0] in ['13a', '14a']:
                            serology.append('6a')
                    ser_dict[allele] = ';'.join(serology)

        serologies = pd.DataFrame(ser_dict.items(), columns=['allele', 'serology'])
        serologies.replace('', np.nan)
        serologies = serologies.set_index('allele')
        serologies.sort_values(by=['serology'], inplace=True, na_position='last')

        ser_handler.close()

        #open all of the current dataframes with old serology info
        oldtrn_frame = pd.read_csv('old_sets/train/' + locus + '_train.csv')
        oldtrn_frame.replace('', np.nan)
        oldtrn_frame = oldtrn_frame.set_index('allele')
        oldval_frame = pd.read_csv('old_sets/train/' + locus + '_validation.csv')
        oldval_frame.replace('', np.nan)
        oldval_frame = oldval_frame.set_index('allele')
        oldtst_frame = pd.read_csv('old_sets/test/' + locus + '_test.csv')
        oldtst_frame = oldtst_frame.replace('', np.nan)
        oldtst_frame = oldtst_frame.set_index('allele')


        #combine the current dataframes all together to make data wrangling work correctly
        train_frame = oldtrn_frame.append(oldval_frame)
        loc_frame = train_frame.append(oldtst_frame)
        loc_frame = loc_frame[~loc_frame.index.duplicated()]

        serologies = serologies[~serologies.index.duplicated()]
        serologies.update(loc_frame, overwrite=True)
        serologies.to_csv('ser/' + locus + '_ser.csv', index=True)
        df = pd.read_csv('aa-matching/output/' + locus + '_AA_poly.csv')
        df = df.set_index('allele')
        df['serology'] = np.nan
        df.update(serologies, overwrite=False)
        df = df.replace('', np.nan)
        df.sort_values(by=['serology'], inplace=True, na_position='last')
        df.to_csv('outframe/' + locus + '_full.csv')

        #combining new serology info with old frame to make new frame
        loc_frame.update(serologies, overwrite=True)
        loc_frame = loc_frame.replace('', np.nan)
        loc_frame.sort_values(by=['serology'], inplace=True, na_position='last')

        droplist = []
        for col in df.columns:
            if col.find('-') != -1:
                droplist.append(col)

        '''
        corfac = 0
        for vally in droplist:
            if df.loc[ref_seq[indy]][vally] == 0:
                droplist.remove(vally)
                corfac += 1
        '''

        fixed_df = df.drop(droplist, axis=1)

        loc_frame.to_csv('full_frames/' + locus + '_result.csv')


        #editing this section to create the new/old dataframes
        tst = fixed_df
        newtst1 = tst[~tst.index.isin(oldtrn_frame.index)]
        newtst = newtst1[~newtst1.index.isin(oldval_frame.index)]
        newtst = newtst[~newtst.index.duplicated()]
        newtst.sort_values(by=['allele'], inplace=True, na_position='last')
        newtst.to_csv('testing/' + locus + '_test.csv')

        trn = fixed_df[fixed_df['serology'] != 'nan']
        trn.sort_values(by=['allele'], inplace=True, na_position='last')
        trn = trn[~trn.index.duplicated()]
        trn = trn[~trn.serology.isnull()]
        newtrn = trn[trn.index.isin(oldtrn_frame.index)]
        newtrn.to_csv('training/' + locus + '_train.csv')
        newval = trn[trn.index.isin(oldval_frame.index)]
        newval.sort_values(by=['allele'], inplace=True, na_position='last')
        newval.to_csv('training/' + locus + '_validation.csv')
    return

def RSNNS_fixer():
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
        fixed_trn.to_csv('./RSNNS_fixed/training/' + locus + '_train.csv', index=False)


        
        val_droplist = []
        for column in new_val.columns:
            if column not in old_val.columns:
                val_droplist.append(column)
        
        fixed_val = new_val.drop(val_droplist, axis=1)
        fixed_val.to_csv('./RSNNS_fixed/training/' + locus + '_validation.csv', index=False)

        tst_droplist = []
        for column in new_tst.columns:
            if column not in old_tst.columns:
                tst_droplist.append(column)
            else:
                next
        
        fixed_tst = new_tst.drop(tst_droplist, axis=1)
        fixed_tst.to_csv('./RSNNS_fixed/testing/' + locus + '_test.csv', index=False)
    return

def pat_generator():
    loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]
    today = date.today()
    date = today.strftime("%B %d, %Y")

    # training files 
    for locus in loci:

        loc_frame = pd.read_csv("./RSNNS_fixed/training/" + locus + "_train.csv")
        
        pat_file = open("./pat_new/" + locus + ".tng.pat", "w+")
        pat_file.write("SNNS pattern definition file V3.2\n")
        pat_file.write("generated at " + date + "\n\n\n")
        pat_file.write("No. of patterns : " + str(len(loc_frame)) + "\n")
        pat_file.write("No. of input units : " + str(len(loc_frame.columns) - 2) + "\n")
        uni_list = loc_frame.serology.unique()

        new_list = []
        for uni in uni_list:
            new = uni.split(';')
            [new_list.append(x) for x in new if x not in new_list]
        
        uni_list = new_list
        for i in range(0,len(uni_list)):
            uni_list[i] = uni_list[i].strip('a')

        pat_file.write("No. of output units : " + str(len(uni_list)) + "\n")

        AAs_list = loc_frame.columns.tolist()
        AAs_list.remove('allele')
        AAs_list.remove('serology')
        polyAAs = ' '.join(AAs_list)
        polyAAs = "# " + polyAAs
        pat_file.write(polyAAs + "\n")
        
        for allele in loc_frame.values.tolist():
            serologies = {}
            for j in uni_list:
                serologies[j] = "0.00"
            #print(allele.pop(0))
            name = allele.pop(0)
            #print(name)
            name = name.split('*')[1]
            sero = allele.pop()
            sero = sero.split(';')
            for i in range(0,len(sero)):
                sero[i] = sero[i].strip('a')
            for value in serologies.keys():
                if value in sero:
                    serologies[value] = "1.00"

            serologies = OrderedDict(sorted(serologies.items(), key = lambda x: ((int(x[0])))))
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
        val_file.write("No. of patterns : " + str(len(val_frame)) + "\n")
        val_file.write("No. of input units : " + str(len(val_frame.columns) - 2) + "\n")
        val_file.write("No. of output units : " + str(len(uni_list)) + "\n")
        
        for val_allele in val_frame.values.tolist():
            serologies = {}
            for j in uni_list:
                serologies[j] = "0.00"
            #val_allele.pop(0)
            val_name = val_allele.pop(0)
            val_name = val_name.split('*')[1]
            val_sero = val_allele.pop()
            val_sero = val_sero.split(';')
            for i in range(0,len(val_sero)):
                val_sero[i] = val_sero[i].strip('a')
            for value in serologies.keys():
                if value in val_sero:
                    serologies[value] = "1.00"
            
            serologies = OrderedDict(sorted(serologies.items(), key = lambda x: ((int(x[0])))))
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
        tst_file.write("No. of patterns : " + str(len(tst_frame)) + "\n")
        tst_file.write("No. of input units : " + str(len(tst_frame.columns) - 2) + "\n")
        tst_file.write("No. of output units : " + str(len(uni_list)) + "\n")
        
        for tst_allele in tst_frame.values.tolist():
            #tst_allele.pop(0)
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
    return

ser_parse()
RSNNS_fixer()
pat_generator()