import pandas as pd
import numpy as np
from fastai import *
from fastai.tabular import *
from _parse.py import *

_file_handler()

tng_df = pd.read_csv('A_train.csv')
tst_df = pd.read_csv('A_test.csv')

dep_var = 'serology'
cat_names = ['allele'] + AAs
procs = [FillMissing, Categorify]


cat_names = ['allele'] + AAs
test = TabularList.from_df(tst_df, path=Path(''), cat_names=cat_names)
data = (TabularList.from_df(df=df, path=Path(''), procs=procs, cat_names=cat_names)
                            .split_by_idx(list(range(65,100)))
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
