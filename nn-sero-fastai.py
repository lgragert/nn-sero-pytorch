import pandas as pd
import numpy as np
from fastai import *
from fastai.tabular import *
from parse import _file_handler
from tqdm import tqdm

#below function directly from Andrew Chang in fast.ai forums (https://forums.fast.ai/t/automated-learning-rate-suggester/44199)
def find_appropriate_lr(model:Learner, lr_diff:int = 15, loss_threshold:float = .05, adjust_value:float = 1, plot:bool = False) -> float:
    #Run the Learning Rate Finder
    model.lr_find()
    
    #Get loss values and their corresponding gradients, and get lr values
    losses = np.array(model.recorder.losses)
    assert(lr_diff < len(losses))
    loss_grad = np.gradient(losses)
    lrs = model.recorder.lrs
    
    #Search for index in gradients where loss is lowest before the loss spike
    #Initialize right and left idx using the lr_diff as a spacing unit
    #Set the local min lr as -1 to signify if threshold is too low
    r_idx = -1
    l_idx = r_idx - lr_diff
    while (l_idx >= -len(losses)) and (abs(loss_grad[r_idx] - loss_grad[l_idx]) > loss_threshold):
        local_min_lr = lrs[l_idx]
        r_idx -= 1
        l_idx -= 1

    lr_to_use = local_min_lr * adjust_value
    
    if plot:
        # plots the gradients of the losses in respect to the learning rate change
        plt.plot(loss_grad)
        plt.plot(len(losses)+l_idx, loss_grad[l_idx],markersize=10,marker='o',color='red')
        plt.ylabel("Loss")
        plt.xlabel("Index of LRs")
        plt.show()

        plt.plot(np.log10(lrs), losses)
        plt.ylabel("Loss")
        plt.xlabel("Log 10 Transform of Learning Rate")
        loss_coord = np.interp(np.log10(lr_to_use), np.log10(lrs), losses)
        plt.plot(np.log10(lr_to_use), loss_coord, markersize=10,marker='o',color='red')
        plt.show()
        
    return lr_to_use

# _file_handler()

loci = ['A', 'B', 'C', 'DPB1', 'DQB1', 'DRB1']
for locus in loci:
  AAs = []
  tng_df = pd.read_csv('RSNNS_fixed/training/' + locus + '_train.csv')
  tst_df = pd.read_csv('RSNNS_fixed/testing/' + locus + '_test.csv')
  val_df = pd.read_csv('RSNNS_fixed/training/' + locus + '_validation.csv')
  tng_idx = len(tng_df)
  val_idx = len(val_df) + 1
  tst_idx = len(tst_df)
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
                              .databunch(bs=tng_idx, val_bs=(val_idx-1)))
                            


  acc_02 = partial(accuracy_thresh, thresh=0.51)
  f_score = partial(fbeta, thresh=0.2)

  learn = tabular_learner(data, layers=[200,50], metrics=[acc_02, f_score])
  print(data.classes)

  lr = find_appropriate_lr(model=learn)
  learn.recorder.plot(suggestion=True)

  learn.fit(15, lr=lr)

  learn.model
  learn.recorder.plot_losses()

  test_id = list(tst_df['allele'])

  classes = data.classes
  predictions = []
  print(classes)

  for i in tqdm(range(0,tst_idx)):
    category = str(learn.predict(tst_df.iloc[i], thresh=0.75)[0])
    sero = category.strip('MultiCategory ')
    sero = sero.replace(';',' ')
    predictions.append(sero.split())


  # below code involved help from some website using fast.ai to demonstrate kaggle solutions
  output_preds = pd.DataFrame({'allele': test_id, 'serology': predictions})
  output_preds.to_csv('predictions/' + locus + '_predictions.csv', index=False)
