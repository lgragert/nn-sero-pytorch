## Name: Giovanni Biagini
## PI: Loren Gragert, PhD
## Institution: Tulane University School of Medicine
## Department: Pathology and Laboratory Medicine
## Location: Louisiana Cancer Research Center
## Date Begun: 01/16/2020
## Purpose: To modernize SNNS

# figure out how to parse the training, testing, and validation files into vectors.


# multi-purpose function to parse each type of file
def _parse(tng_file, tst_file, val_file):
    tng_dict = {}
    tst_dict = {}
    val_dict = {}
    AA_dict = {}

    # loop through each line of the training file
    for each in tng_file:
        # will only be executed once, to create a list of the polymorphic amino acids
        if each.find("No. of output") != -1:
            each = next(tng_file)
            AA = each.strip("# ")
            AA_list = AA.split()
            continue
        # main purpose, identifying the hashtags at the beginning of the significant lines
        if each.find('#') != -1:
            # gathering information on a specific allele
            if (each.find('# input') == -1) & (each.find('# output') == -1):
                temp = each.strip('# ')
                allele = temp.rstrip()
                continue
            # looking for the binary values corresponding to the amino acid polymorphisms
            elif each.find('# input') != -1:
                line = next(tng_file)
                # generating a list of the binary values
                bin_val = line.split()
                # the values for the specific amino acid are zipped into a dictionary with the polymorphic AAs as the keys
                AA_dict = dict(zip(AA_list, bin_val))
                # the dictionary of AA polymorphisms is added to the (nested) dictionary of all alleles, with its specific allele as its key
                tng_dict[allele] = AA_dict

    # basically the same, but to parse the testing file
    for each in tst_file:
        if each.find('#') != -1:
            if (each.find('# input') == -1) & (each.find('# output') == -1):
                temp = each.strip('# testing ')
                allele = temp.rstrip()
                continue
            elif each.find('# input') != -1:
                line = next(tst_file)
                bin_val = line.split()
                AA_dict = dict(zip(AA_list, bin_val))
                tst_dict[allele] = AA_dict

     # final loop to parse the validation file       
    for each in val_file:
        if each.find('#') != -1:
            if (each.find('# input') == -1) & (each.find('# output') == -1):
                temp = each.strip('# ')
                allele = temp.rstrip()
                continue
            elif each.find('# input') != -1:
                line = next(val_file)
                bin_val = line.split()
                AA_dict = dict(zip(AA_list, bin_val))
                val_dict[allele] = AA_dict
    return(tng_dict, tst_dict, val_dict)

def _file_handler():
    # opening files to send to the parser
    # there is almost definitely a much simpler way to code this

    loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]
    output_list = []
    
    for each in loci:
        training_file = open(each + ".tng.pat", 'r')
        testing_file = open(each + ".tst.pat", 'r')
        validation_file = open(each + ".val.pat", 'r')
        tng_dict, tst_dict, val_dict = _parse(training_file, testing_file, validation_file)
        output_list.append(tng_dict)
        output_list.append(tst_dict)
        output_list.append(val_dict)
        training_file.close()
        testing_file.close()
        validation_file.close()

    A = output_list[0:3]
    B = output_list[3:6]
    C = output_list[6:9]
    DPB1 = output_list[9:12]
    DQB1 = output_list[12:15]
    DRB1 = output_list[15:18]

    return(A, B, C, DPB1, DQB1, DRB1)
