import pandas as pd
import numpy as np

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
    df = pd.read_csv('aa_matching/output/' + locus + '_AA_poly.csv')
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

    newtst = fixed_df[fixed_df['serology'] == 'nan']
    newtst = newtst.append(fixed_df[fixed_df.serology.isnull()])
    newtst.sort_values(by=['allele'], inplace=True, na_position='last')
    newtst.to_csv('testing/' + locus + '_test.csv')

    trn = fixed_df[fixed_df['serology'] != 'nan']
    trn.sort_values(by=['allele'], inplace=True, na_position='last')
    trn = trn[~trn.index.duplicated()]
    trn = trn[~trn.serology.isnull()]
    newtrn = trn[trn.index.isin(oldtrn_frame.index)]
    newtrn.to_csv('training/' + locus + '_train.csv')
    newval = trn[~trn.index.isin(oldtrn_frame.index)]
    newval.sort_values(by=['allele'], inplace=True, na_position='last')
    newval.to_csv('training/' + locus + '_validation.csv')
