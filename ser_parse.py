import pandas as pd
import numpy as np

current = '3370'
loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]

for locus in loci:
    ser_handler = open('rel_dna_ser/rel_dna_ser_' + current + '.txt', 'r')
    loc_test = locus + '*'
    ser_dict = {}
    for line in ser_handler:
        if line.find('*') != -1:
            info = line.split(';')
            if info[0] == loc_test:
                numsplit = info[1].split(':')
                nums = numsplit[0] + ':' + numsplit[1]
                allele = info[0] + nums
                if info[2] == '?':
                    info[2] = ''
                elif info[2] == '0':
                    info[2] = ''
                elif info[2] == '':
                    info[2] == ''
                else:
                    info[2] = info[2] + 'a'
                serology = str(info[2])
                ser_dict[allele] = serology

    serologies = pd.DataFrame(ser_dict.items(), columns=['allele', 'serology'])
    serologies.replace('', np.nan)
    serologies.to_csv('ser/' + locus + '_ser.csv', index=False)

    ser_handler.close()

    #open all of the current dataframes with old serology info
    trn_frame = pd.read_csv('train/' + locus + '_train.csv')
    trn_frame.replace('', np.nan)
    val_frame = pd.read_csv('train/' + locus + '_validation.csv')
    val_frame.replace('', np.nan)
    tst_frame = pd.read_csv('test/' + locus + '_test.csv')
    tst_frame.replace('', np.nan)

    #combine the current dataframes all together to make data wrangling work correctly
    train_frame = trn_frame.append(val_frame)
    loc_frame = train_frame.append(tst_frame)

    #dropping the serology column from the old frame in order to combine with the new serology information
    #full_frame = full_frame.drop(columns=['serology'])

    #combining new serology info with old frame to make new frame
    loc_frame.update(serologies, overwrite=False)
    loc_frame.replace('', 'NaN')

    ordered = loc_frame.sort_values(by=['serology'])
    print(ordered.count())

    loc_frame.to_csv('full_frames/' + locus + '_result.csv', index=False)

