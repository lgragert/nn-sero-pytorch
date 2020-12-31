import pandas as pd
import numpy as np
from fastai import *
from fastai.tabular import *
#from parse import _file_handler
from tqdm import tqdm



  # function to check if value can be an integer - to eliminate excess characters from serology labels
def checkInt(x):
    try:
        int(x)
        return True
    except ValueError:
        return False

# _file_handler()
base_dir = '/work/gbiagini/nn-sero-pytorch/'

loci = ['A', 'B', 'C', 'DPB1', 'DQB1', 'DRB1']
#loci = ["A"]
epoch = { "A":150, "B":22, "C":22, "DPB1":100, "DQB1":100, "DRB1":100 }
layer = { "A":[7500], "B":[2000, 1000], "C":[1000, 500], "DPB1":[700], "DQB1":[700], "DRB1":[700] }


for locus in loci:
  
  AAs = []
  tng_df = pd.read_csv(base_dir + 'RSNNS_fixed/training/' + locus + '_train.csv')
  tst_df = pd.read_csv(base_dir + 'RSNNS_fixed/testing/' + locus + '_test.csv')
  val_df = pd.read_csv(base_dir + 'RSNNS_fixed/training/' + locus + '_validation.csv')
  tng_idx = len(tng_df)
  val_len = len(val_df)
  val_idx = tng_idx + val_len
  tst_idx = len(tst_df)
  tbatch = tng_idx // 1
  if (tbatch <= 1):
    tbatch = tng_idx // 1
  vbatch = val_idx // 1
  if (vbatch <= 1):
    vbatch = val_idx // 1
  

  df = tng_df.append(val_df)

  for each in tng_df:
    if (each != 'allele') & (each != 'serology'):
      AAs.append(each)

  dep_var = 'serology'
  cat_names = ['allele'] + AAs
  procs = [FillMissing, Categorify]


  cat_names = ['allele'] + AAs
  test = TabularList.from_df(tst_df, path=Path(''), cat_names=cat_names)
  
  data = (TabularList.from_df(df=df, path=Path(''), procs=procs, cat_names=cat_names)
                              .split_by_idx(list(range(tng_idx,val_idx)))
                              .label_from_df(cols=dep_var, label_delim=';')
                              .add_test(test)
                              .databunch(bs=tbatch, val_bs=vbatch))
  '''
  data = (TabularList.from_df(df=df, path=Path(''), procs=procs, cat_names=cat_names)
                              .split_by_idx(list(range(tng_idx,val_idx)))
                              .label_from_df(cols=dep_var, label_delim=';')
                              .add_test(test)
                              .databunch(bs=tng_idx, val_bs=val_idx))
  '''
  acc_02 = partial(accuracy_thresh, thresh=0.99)
  f_score = partial(fbeta, thresh=0.55)

  pre_vote = {}
  avg_pred = {}
  all_models = []

  for n in range(1,2):
    #weights = torch.ones([data.c]).float().cuda()
    #loss = nn.BCEWithLogitsLoss(pos_weight=weights)
    #learn = tabular_learner(data, opt_func=optim.SGD, layers=layer[locus], metrics=[acc_02, f_score], loss_func=loss)
    learn = tabular_learner(data, opt_func=optim.SGD, layers=layer[locus], metrics=[acc_02, f_score])
    print(data.classes)

    lr = 0.5
    #learn.recorder.plot(suggestion=True)

    #learn.fit_one_cycle(epoch[locus], max_lr=slice(lr))
    learn.fit_one_cycle(epoch[locus], max_lr=slice(lr))
    #learn.fit_one_cycle(epoch[locus], 0.2)
    learn.model
    #learn.recorder.plot_losses()

    test_id = list(tst_df['allele'])

    classes = data.classes
    predictions = []
    print(classes)

   
    for i in tqdm(range(0,tst_idx)):
      category = str(learn.predict(tst_df.iloc[i], thresh=0.55)[0])
      sero = category.strip('MultiCategory ')
      sero = sero.replace(';',' ')
      sero = sero.replace('a','')
      predictions.append(sero.split())
    '''
    category = learn.get_preds()
    print(category)
    '''
    pre_vote = {test_id[j]: str(predictions[j]) for j in range(len(test_id)) }
    avg_pred[str(n)] = pre_vote
    

  avg_frame = pd.DataFrame.from_dict(avg_pred)
  mode = avg_frame.mode(axis=1)
  rmode = mode[0]
  rmode = pd.DataFrame(rmode)
  rmode = rmode.reset_index()
  rmode.columns=['allele', 'serology']

  #output_preds = pd.DataFrame({'allele': test_id, 'serology': predictions})
  rmode.to_csv(base_dir + 'predictions/' + locus + '_predictions.csv', index=False)


  newDict = {}
  simDict = {}
  diffDict = {}
  oldPredict = {}
  newPredict = {}
  oldPredFile = base_dir + "old-predictions/" + loc + ".chile"
  newPreds = pd.read_csv(base_dir + "predictions/" + loc + "_predictions.csv")
  newPreds = newPreds.set_index('allele')
  newPreds = newPreds.to_dict()
  newPredict = newPreds["serology"]
  for nKey in newPredict.keys():
      adjustMe = newPredict[nKey]
      adjustMe = adjustMe.replace('[','')
      adjustMe = adjustMe.replace(']','')
      adjustMe = adjustMe.replace(' ','')
      adjustMe = adjustMe.replace("'",'')
      adjustMe = adjustMe.split(',')
      newPredict[nKey] = [x.strip('a') for x in adjustMe if checkInt(x)]
  with open(oldPredFile, "r") as handle:
      for line in handle:
          if line.find('%') == -1:
            next
          else:
              line = line.split()
              if line == []:
                  next
              else:
                  line[:] = [x for x in line if (x != '[100.00%]')]
                  allele = loc + "*" + str(line[0][:-1])
                  oldPredict[allele] = line[1:]


  for each in oldPredict.keys():
      allDict = {}
      allDict["Allele"] = each
      allDict["Old Assignment"] = oldPredict[each]
      allDict["New Assignment"] = newPredict[each]
      if each not in newPredict.keys():
        next
      elif set(newPredict[each]) != set(oldPredict[each]):
          diffDict[each] = allDict
      elif set(newPredict[each]) == set(oldPredict[each]):
          simDict[each] = allDict
  diffFrame = pd.DataFrame.from_dict(diffDict)
  diffFrame = diffFrame.transpose()
  diffFrame.to_csv(base_dir + "comparison/" + loc + "_compfile.csv", index=False)
  simFrame = pd.DataFrame.from_dict(simDict)
  simFrame = simFrame.transpose()
  simFrame.to_csv(base_dir + "comparison/" + loc + "_similar.csv", index=False)
  

  for allele in newPredict.keys():
      allDict = {}
      allDict["Allele"] = allele
      allDict["Serologic Assignment"] = newPredict[allele]
      if allele not in oldPredict.keys():
          newDict[allele] = allDict
  newFrame = pd.DataFrame.from_dict(simDict)
  newFrame = newFrame.transpose()
  newFrame.to_csv(base_dir + "comparison/" + loc + "_newsies.csv", index=False)

  writefile = open(base_dir + "comparison/" + loc + "_concordance.txt", "w+")

  simLen = len(simFrame)
  diffLen = len(diffFrame)
  writefile.write("HLA-" +loc+ " Similar: " + str(simLen))
  writefile.write("HLA-" +loc+ " Different: " + str(diffLen))

  concordance = (simLen / (simLen + diffLen)) * 100
  writefile.write("HLA-" +loc+ " Concordance: " + str(concordance) + "%")
  writefile.close()