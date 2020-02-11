import pandas as pd

# multi-purpose function to parse each type of file
def _parse(tng_file, tst_file, val_file, locus):
    tng_dict = {}
    tst_dict = {}
    val_dict = {}
    AA_dict = {}

    # loop through each line of the training file
    for each in tng_file:
        #values = []
        AAs = []
        trick = {}
        # will only be executed once, to create a list of the polymorphic amino acids
        if each.find("No. of output") != -1:
            each = next(tng_file)
            AA = each.strip("# ")
            AA_list = AA.split()
            continue
        # main purpose, identifying the hashtags at the beginning of the significant lines
        if each.find("#") != -1:
            # gathering information on a specific allele
            if (each.find("# input") == -1) & (each.find("# output") == -1):
                temp = each.strip("# ")
                allele = temp.rstrip()
                allele = locus + '*' + allele
                continue
            # looking for the binary values corresponding to the amino acid polymorphisms
            elif each.find("# input") != -1:
                line = next(tng_file)
                # generating a list of the binary values
                bin_val = line.split()
                bin_val = list(map(int, bin_val))
                line = next(tng_file)
                out_lines_idx = line
                specificities_n = out_lines_idx.strip("# output")
                specificities_f = specificities_n.rstrip()
                specificities = specificities_f.split()
                line = next(tng_file)
                out_vals = line.strip(' ').split()
                out_vals = list(map(float,out_vals))
                out_dict = dict(zip(specificities, out_vals))
                for spec in list(out_dict):
                  if out_dict[spec] == 0.00:
                    del(out_dict[spec])
                serology = list(out_dict)
                serology = list(map(str,serology))
                for val in range(len(serology)):
                  serology[val] += 'a'
                #values.append(out_dict)
                spacer = ';'
                serology = spacer.join(serology)
                AAs = ['allele'] + AA_list
                values = [allele] + bin_val
                trick = dict(zip(AAs, values))
                trick['serology'] = serology
                tng_dict[allele] = trick

    # basically the same, but to parse the testing file
    for each in tst_file:
        #values = []
        AAs = []
        trick = {}
        if each.find("#") != -1:
            if (each.find("# input") == -1) & (each.find("# output") == -1):
                temp = each.strip("# testing ")
                allele = temp.rstrip()
                allele = locus + '*' + allele
                continue
            elif each.find("# input") != -1:
                line = next(tst_file)
                bin_val = line.split()
                bin_val = list(map(int,bin_val))
                AAs = ['allele'] + AA_list + ['serology']
                values = [allele] + bin_val
                trick = dict(zip(AAs, values))
                trick['serology'] = None
                tst_dict[allele] = trick

     # final loop to parse the validation file       
    for each in val_file:
        #values = []
        AAs = []
        trick = {}
        if each.find("#") != -1:
            if (each.find("# input") == -1) & (each.find("# output") == -1):
                temp = each.strip("# ")
                allele = temp.rstrip()
                allele = locus + '*' + allele
                continue
            elif each.find("# input") != -1:
                line = next(val_file)
                bin_val = line.split()
                bin_val = list(map(int,bin_val))
                line = next(val_file)
                out_lines_idx = line
                specificities_n = out_lines_idx.strip("# output")
                specificities_f = specificities_n.rstrip()
                specificities = specificities_f.split()
                line = next(val_file)
                out_vals = line.strip(' ').split()
                out_vals = list(map(float,out_vals))
                out_dict = dict(zip(specificities, out_vals))
                for spec in list(out_dict):
                  if out_dict[spec] == 0.00:
                    del(out_dict[spec])
                serology = list(out_dict)
                serology = list(map(str,serology))
                for val in range(len(serology)):
                  serology[val] += 'a'
                #values.append(out_dict)
                spacer = ';'
                serology = spacer.join(serology)
                AAs = ['allele'] + AA_list
                values = [allele] + bin_val
                trick = dict(zip(AAs, values))
                trick['serology'] = serology

                val_dict[allele] = trick
    return(tng_dict, tst_dict, val_dict)

def _file_handler():

    loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]
    output_list = []
    
    for locus in loci:
        tng_AAs = []
        tst_AAs = []
        val_AAs = []
        AAs = []

        training_file = open(locus + ".tng.pat", 'r')
        testing_file = open(locus + ".tst.pat", 'r')
        validation_file = open(locus + ".val.pat", 'r')
        tng_dict, tst_dict, val_dict = _parse(training_file, testing_file, validation_file, locus)

        for key in tng_dict.keys():
            tng_AAs.append(tng_dict[key])

        for tst_key in tst_dict.keys():
            tst_AAs.append(tst_dict[tst_key])

        for val_key in val_dict.keys():
            val_AAs.append(val_dict[val_key])

        for one in tng_AAs[0].keys():
            if (one != 'serology') & (one != 'allele'):
                AAs.append(one)


        tng_frame = pd.DataFrame(data=tng_AAs)
        val_frame = pd.DataFrame(data=val_AAs)
        tst_frame = pd.DataFrame(data=tst_AAs)
        tng_frame.to_csv(locus + '_train.csv', index=False)
        val_frame.to_csv(locus + '_validation.csv', index=False)
        tst_frame.to_csv(locus + '_test.csv', index=False)
        training_file.close()
        testing_file.close()
        validation_file.close()

    return()

_file_handler()