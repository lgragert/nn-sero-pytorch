#base_dir = '/home/gbiagini/dev/nn-sero-pytorch/'
base_dir = './'

import pandas as pd
import numpy as np
from fastai import *
from fastai.tabular.all import *
#from parse import _file_handler
from tqdm import tqdm



# _file_handler()

#loci = ['A', 'B', 'C', 'DPB1', 'DQB1', 'DRB1']
loci = ["A"]
epoch = { "A":22, "B":15, "C":13, "DPB1":13, "DQB1":20, "DRB1":13 }
layer = { "A":[1500], "B":[2500], "C":[1500], "DPB1":[1500], "DQB1":[2000], "DRB1":[1500] }


for locus in loci:
  
  AAs = []
  #tng_df = pd.read_csv(base_dir + 'RSNNS_fixed/training/' + locus + '_train.csv')
  #tst_df = pd.read_csv(base_dir + 'RSNNS_fixed/testing/' + locus + '_test.csv')
  #val_df = pd.read_csv(base_dir + 'RSNNS_fixed/training/' + locus + '_validation.csv')
  tng_df = pd.read_csv(base_dir + 'old_sets/train/' + locus + '_train.csv')
  tst_df = pd.read_csv(base_dir + 'old_sets/test/' + locus + '_test.csv')
  val_df = pd.read_csv(base_dir + 'old_sets/train/' + locus + '_validation.csv')
  tng_idx = len(tng_df)
  val_idx = len(val_df) + 1
  tst_idx = len(tst_df)
  tbatch = tng_idx // 2
  if (tbatch <= 1):
    tbatch = tng_idx // 2
  vbatch = val_idx // 2
  if (vbatch <= 1):
    vbatch = val_idx // 2
  
  df = tng_df.append(val_df)

  for each in tng_df:
    x = re.search("(^\w\d*$)", each)
    if x:
      AAs.append(each)
  
  df['serology'] = df['serology'].astype(str)

  cat_names = ['allele'] + AAs
  #dep_var = [x for x in tng_df.columns if (x not in cat_names)]
  dep_var = 'serology'
  cat_names = cat_names
  procs = [FillMissing, Categorify]
  splits = [list(range(0,tng_idx)), list(range(tng_idx,val_idx))]
  
  '''
  def get_x(r): return r['allele']
  def get_y(r): return r['serology'].split(' ')
  dblock = DataBlock(get_x = get_x, get_y = get_y)
  dblock = DataBlock(blocks=MultiCategoryBlock, get_x=get_x, get_y=get_y)
  dsets=dblock.datasets(df)
  print(dsets.train[0])
  '''

  cat_names = ['allele'] + AAs
  test = TabularPandas(tst_df, cat_names=cat_names)
  
  to = (TabularPandas(df=df, procs=procs, cat_names=cat_names, splits=splits, y_names=dep_var))

  dls = to.dataloaders(bs=tbatch, path=".", block_y=MultiCategoryBlock)
  print(dls)
  
  '''
  #acc_02 = partial(accuracy_thresh, thresh=0.99)
  #f_score = partial(fbeta, thresh=0.51)
  f_score = partial(FBetaMulti(0.51), thresh=0.51)
  acc_02 = partial(accuracy_multi, thresh=0.99)

  pre_vote = {}
  avg_pred = {}
  all_models = []

  for n in range(1,2):
    learn = tabular_learner(dls, opt_func = SGD, layers=layer[locus],
                              metrics=[acc_02, f_score],
                              loss_func=BCELossFlat)


    lr = 0.5
    #learn.recorder.plot(suggestion=True)

    #learn.fit(epoch[locus], lr)
    #learn.fit_one_cycle(epoch[locus], lr)
    learn.fit_one_cycle(epoch[locus])
    learn.model
    #learn.recorder.plot_losses()

    test_id = list(tst_df['allele'])

    classes = data.classes
    predictions = []
    print(classes)

   
    for i in tqdm(range(0,tst_idx)):
      category = str(learn.predict(tst_df.iloc[i], thresh=0.51)[0])
      sero = category.strip('MultiCategory ')
      sero = sero.replace(';',' ')
      sero = sero.replace('a','')
      predictions.append(sero.split())
    

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
  '''