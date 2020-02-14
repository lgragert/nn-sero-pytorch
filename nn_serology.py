import pandas as pd
import numpy as np
from fastai import *
from fastai.tabular import *
from parse import *

#_file_handler()

quest = 'y'
while (quest == 'y') or (quest == 'Y'):
  locus = input("Locus? (A/B/C/DPB1/DQB1/DRB1) ")

  tng_df = pd.read_csv('train/' + locus + '_train.csv')
  tst_df = pd.read_csv('test/' + locus + '_test.csv')
  val_df = pd.read_csv('train/' + locus + '_validation.csv')
  tng_idx = len(tng_df)
  val_idx = len(val_df) + 1
  tst_idx = len(tst_df)
  df = tng_df.append(val_df)

  AAs = []
  for each in tng_df:
    if (each != 'allele') & (each != 'serology'):
      AAs.append(each)

  dep_var = 'serology'
  cat_names = ['allele'] + AAs
  procs = [FillMissing, Categorify]

  test = TabularList.from_df(tst_df, path=Path(''), cat_names=cat_names, cont_names=None)
  data = (TabularList.from_df(df=df, path=Path(''), procs=procs, cat_names=cat_names, cont_names=None)
                              .split_by_idx(list(range(tng_idx,(tng_idx+val_idx-1))))
                              .label_from_df(cols=dep_var, label_delim=';')
                              .add_test(test)
                              .databunch(bs=37, num_workers=0))                          

  acc_02 = partial(accuracy_thresh, thresh=0.501)
  f_score = partial(fbeta, thresh=0.501)

  trainio = 'y'
  while (trainio == 'y') or (trainio == 'Y'): 
    layers = input('Size of layer? (Int) ')
    learn = tabular_learner(data, layers=[int(layers)], loss_func=MSELossFlat(), metrics=[acc_02,f_score])

    epochs = int(input('# of Epochs? (Int) '))
    lr = float(input('Learning Rate? (Int or Float) '))

    learn.fit(epochs, lr=lr)
    
    trainio = input('Reset and train again? (y/n) ')

  print('predicting on ' + locus + '...')
  test_id = list(tst_df['allele'])
  for j in range(0,len(test_id)):
    test_id[j] = test_id[j] 

  classes = data.classes
  predictions = []

  for i in range(0, tst_idx):
    category = str(learn.predict(tst_df.iloc[i], thresh=0.35)[0])
    sero = category.strip('MultiCategory ')
    sero = sero.replace(';',' ')
    predictions.append(sero.split())

  # below code involved help from some website using fast.ai to demonstrate kaggle solutions
  output_preds = pd.DataFrame({'allele': test_id, 'serology': predictions})
  output_preds.to_csv('snns_predictions/' + locus + '_predictions.csv', index=False)

  quest = input('Predict for another locus? (y/n) ')