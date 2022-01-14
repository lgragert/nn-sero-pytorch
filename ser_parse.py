import pandas as pd
import numpy as np


def ser_v1():
    current = '3430'
    loci = ["A", "B", "C", "DQB1", "DRB1"]
    for locus in loci:
        ser_handler = open('./rel_dna_ser/rel_dna_ser_' + current + '.txt', 'r')
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

                    #check if null allele - append if not already noted (as for two-field)
                    if (info[1][-1].isalpha() == True) & (allele[-1].isalpha() == False):
                        allele += info[1][-1]
                    if info[2] == '?':
                        info[2] = np.nan
                    elif info[2] == '0':
                        info[2] = np.nan
                    elif info[2] == '':
                        info[2] = np.nan
                    else:
                        info[2] = info[2]
                    serology.append(str(info[2]))
                    ser_dict[allele] = ';'.join(serology)

        serologies = pd.DataFrame(ser_dict.items(), columns=['allele', 'serology'])
        serologies.replace('', np.nan)
        serologies = serologies.set_index('allele')
        serologies.sort_values(by=['serology'], inplace=True, na_position='last')

        ser_handler.close()

        #open all of the current dataframes with old serology info
        oldtrn_frame = pd.read_csv('./old_sets/train/' + locus + '_train.csv')
        oldtrn_frame.replace('', np.nan)
        oldtrn_frame = oldtrn_frame.set_index('allele')
        oldval_frame = pd.read_csv('./old_sets/train/' + locus + '_validation.csv')
        oldval_frame.replace('', np.nan)
        oldval_frame = oldval_frame.set_index('allele')
        oldtst_frame = pd.read_csv('./old_sets/test/' + locus + '_test.csv')
        oldtst_frame = oldtst_frame.replace('', np.nan)
        oldtst_frame = oldtst_frame.set_index('allele')


        #combine the current dataframes all together to make data wrangling work correctly
        train_frame = oldtrn_frame.append(oldval_frame)
        loc_frame = train_frame.append(oldtst_frame)
        loc_frame = loc_frame[~loc_frame.index.duplicated()]

        serologies = serologies[~serologies.index.duplicated()]
        serologies.update(loc_frame, overwrite=True)
        serologies.to_csv('./ser/' + locus + '_ser.csv', index=True)
        df = pd.read_csv('./aa-matching/output/' + locus + '_AA_poly.csv')
        df = df.set_index('allele')
        df['serology'] = np.nan
        df.update(serologies, overwrite=False)
        df = df.replace('', np.nan)
        df.sort_values(by=['serology'], inplace=True, na_position='last')
        df.to_csv('./outframe/' + locus + '_full.csv')

        #combining new serology info with old frame to make new frame
        loc_frame.update(serologies, overwrite=True)
        loc_frame = loc_frame.replace('', np.nan)
        loc_frame.sort_values(by=['serology'], inplace=True, na_position='last')

        droplist = []
        for col in df.columns:
            if col.find('-') != -1:
                droplist.append(col)

        fixed_df = df.drop(droplist, axis=1)

        loc_frame.to_csv('./full_frames/' + locus + '_result.csv')

        newtst = fixed_df[fixed_df['serology'] == 'nan']
        newtst = newtst.append(fixed_df[fixed_df.serology.isnull()])
        newtst.sort_values(by=['allele'], inplace=True, na_position='last')
        newtst.to_csv('./testing/' + locus + '_test.csv')

        trn = fixed_df[fixed_df['serology'] != 'nan']
        trn.sort_values(by=['allele'], inplace=True, na_position='last')
        trn = trn[~trn.index.duplicated()]
        trn = trn[~trn.serology.isnull()]
        newtrn = trn[trn.index.isin(oldtrn_frame.index)]
        newtrn.to_csv('./training/' + locus + '_train.csv')
        newval = trn[~trn.index.isin(oldtrn_frame.index)]
        newval.sort_values(by=['allele'], inplace=True, na_position='last')
        newval.to_csv('./training/' + locus + '_validation.csv')
    return

def ser_v2():
    current = '3430'
    loci = ["A", "B", "C", "DQB1", "DRB1"]
    for locus in loci:
        ser_handler = open('./rel_dna_ser/rel_dna_ser_' + current + '.txt', 'r')
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

                    #check if null allele - append if not already noted (as for two-field)
                    if (info[1][-1].isalpha() == True) & (allele[-1].isalpha() == False):
                        allele += info[1][-1]
                    if info[2] == '?':
                        info[2] = np.nan
                    elif info[2] == '0':
                        info[2] = np.nan
                    elif info[2] == '':
                        info[2] = np.nan
                    else:
                        info[2] = info[2]
                    serology.append(str(info[2]))
                    ser_dict[allele] = ';'.join(serology)

        serologies = pd.DataFrame(ser_dict.items(), columns=['allele', 'serology'])
        serologies.replace('', np.nan)
        serologies = serologies.set_index('allele')
        serologies.sort_values(by=['serology'], inplace=True, na_position='last')

        ser_handler.close()

        serologies = serologies[~serologies.index.duplicated()]
        serologies.to_csv('./ser/' + locus + '_ser.csv', index=True)
        df = pd.read_csv('./aa-matching/output/' + locus + '_AA_poly.csv')
        df = df.set_index('allele')
        df['serology'] = np.nan
        df.update(serologies, overwrite=False)
        df = df.replace('', np.nan)
        df.sort_values(by=['serology'], inplace=True, na_position='last')
        df.to_csv('./outframe/' + locus + '_full.csv')

        '''
        #combining new serology info with old frame to make new frame
        loc_frame.update(serologies, overwrite=True)
        loc_frame = loc_frame.replace('', np.nan)
        loc_frame.sort_values(by=['serology'], inplace=True, na_position='last')

        droplist = []
        for col in df.columns:
            if col.find('-') != -1:
                droplist.append(col)

        fixed_df = df.drop(droplist, axis=1)

        loc_frame.to_csv('./full_frames/' + locus + '_result.csv')

        newtst = fixed_df[fixed_df['serology'] == 'nan']
        newtst = newtst.append(fixed_df[fixed_df.serology.isnull()])
        newtst.sort_values(by=['allele'], inplace=True, na_position='last')
        newtst.to_csv('./testing/' + locus + '_test.csv')

        trn = fixed_df[fixed_df['serology'] != 'nan']
        trn.sort_values(by=['allele'], inplace=True, na_position='last')
        trn = trn[~trn.index.duplicated()]
        trn = trn[~trn.serology.isnull()]
        newtrn = trn[trn.index.isin(oldtrn_frame.index)]
        newtrn.to_csv('./training/' + locus + '_train.csv')
        newval = trn[~trn.index.isin(oldtrn_frame.index)]
        newval.sort_values(by=['allele'], inplace=True, na_position='last')
        newval.to_csv('./training/' + locus + '_validation.csv')
        '''
    return

ser_v2()