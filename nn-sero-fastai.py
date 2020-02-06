import pandas as pd
import numpy as np
from fastai import *
from fastai.tabular import *


# multi-purpose function to parse each type of file
def _parse(tng_file, tst_file, val_file):
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
                continue
            # looking for the binary values corresponding to the amino acid polymorphisms
            elif each.find("# input") != -1:
                line = next(tng_file)
                # generating a list of the binary values
                bin_val = line.split()
                bin_val = list(map(int, bin_val))
                # the values for the specific amino acid are zipped into a dictionary with the polymorphic AAs as the keys
                ##AA_dict = dict(zip(AA_list, bin_val))
                ##values.append(AA_dict)
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
                continue
            elif each.find("# input") != -1:
                line = next(tst_file)
                bin_val = line.split()
                bin_val = list(map(int,bin_val))
                #AA_dict = dict(zip(AA_list, bin_val))
                #values.append(AA_dict)
                '''
                line = next(tst_file)
                out_lines_idx = line
                specificities_n = out_lines_idx.strip("# output")
                specificities_f = specificities_n.rstrip()
                specificities = specificities_f.split()
                line = next(tst_file)
                out_vals = line.split()
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
                '''
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
                continue
            elif each.find("# input") != -1:
                line = next(val_file)
                bin_val = line.split()
                bin_val = list(map(int,bin_val))
                #AA_dict = dict(zip(AA_list, bin_val))
                #values.append(AA_dict)
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
                print(serology)
                AAs = ['allele'] + AA_list
                values = [allele] + bin_val
                trick = dict(zip(AAs, values))
                trick['serology'] = serology

                val_dict[allele] = trick
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




A, B, C, DPB1, DQB1, DRB1 = _file_handler()

A_tng = A[0]
A_tst = A[1]
A_val = A[2]
tng_alleles = []
tst_alleles = []
val_alleles = []
tng_AAs = []
tst_AAs = []
val_AAs = []
AAs = []
specificities = []

for key in A_tng.keys():
  tng_alleles.append(key)
  tng_AAs.append(A_tng[key])
  for each_dict in tng_AAs:
    for each in each_dict:
      if each_dict[each] == 0:
        each_dict[each] = False
      elif each_dict[each] == 1:
        each_dict[each] = True

for tst_key in A_tst.keys():
  tst_alleles.append(tst_key)
  tst_AAs.append(A_tst[tst_key])
  for each_dict in tst_AAs:
    for each in each_dict:
      if each_dict[each] == 0:
        each_dict[each] = False
      elif each_dict[each] == 1:
        each_dict[each] = True

for val_key in A_val.keys():
  val_alleles.append(val_key)
  val_AAs.append(A_val[val_key])
  for each_dict in val_AAs:
    for each in each_dict:
      if each_dict[each] == 0:
        each_dict[each] = False
      elif each_dict[each] == 1:
        each_dict[each] = True

for each in tng_AAs[0].keys():
  if (each != 'serology') & (each != 'allele'):
    AAs.append(each)


tng_df = pd.DataFrame(data = tng_AAs)
val_df = pd.DataFrame(data = val_AAs)
df = tng_df.append(val_df)
tng_idx = len(tng_df)
val_idx = len(val_df)
tst_df = pd.DataFrame(data = tst_AAs)

dep_var = 'serology'
cat_names = ['allele'] + AAs
procs = [FillMissing, Categorify]


cat_names = ['allele'] + AAs
test = TabularList.from_df(tst_df, path=Path(''), cat_names=cat_names)
data = (TabularList.from_df(df=df, path=Path(''), procs=procs, cat_names=cat_names)
                            .split_by_idx(list(range(tng_idx,tng_idx+val_idx)))
                            .label_from_df(cols=dep_var, label_delim=';')
                            .add_test(test)
                            .databunch(bs=37))
                          

data.show_batch(rows = 37)

acc_02 = partial(accuracy_thresh, thresh=0.2)
f_score = partial(fbeta, thresh=0.2)

learn = tabular_learner(data, layers=[30,30,30], metrics=[acc_02, f_score])
print(data.classes)

learn.lr_find()
learn.recorder.plot(suggestion=True)

learn.fit(75, lr=3.63e-2)

learn.model
learn.recorder.plot_losses()

learn.show_results(rows=10)
interp = ClassificationInterpretation.from_learner(learn)
interp.most_confused(min_val=5)

test_id = list(tst_df['allele'])
for j in range(0,len(test_id)):
  test_id[j] = 'A*' + test_id[j] 

classes = data.classes
predictions = []

for i in range(0,1692):
  category = str(learn.predict(tst_df.iloc[i], thresh=0.40)[0])
  sero = category.strip('MultiCategory ')
  sero = sero.replace(';',' ')
  predictions.append(sero.split())

# below code involved help from some website using fast.ai to demonstrate kaggle solutions
output_preds = pd.DataFrame({'allele': test_id, 'serology': predictions})
output_preds.to_csv('predictions.csv', index=False)
