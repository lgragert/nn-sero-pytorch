#base_dir = '/home/gbiagini/dev/nn-sero-pytorch/'
base_dir = './'

import pandas as pd
import numpy as np
import torch
import random
from fastai import *
from fastai.basics import *
from fastai.tabular import *

def random_seed(seed_value, use_cuda):
    np.random.seed(seed_value) # cpu vars
    torch.manual_seed(seed_value) # cpu  vars
    random.seed(seed_value) # Python
    if use_cuda: 
        torch.cuda.manual_seed(seed_value)
        torch.cuda.manual_seed_all(seed_value) # gpu vars
        torch.backends.cudnn.deterministic = True  #needed
        torch.backends.cudnn.benchmark = False
    
    return

loci = ['A', 'B', 'C', 'DQB1', 'DRB1']
#loci = ['A']

# function to check if value can be an integer - to eliminate excess characters from serology labels
def metrics(print_all='no'):
    loci = ['A', 'B', 'C', 'DQB1', 'DRB1']
    #loci = ['A']

    # function to check if value can be an integer - to eliminate excess characters from serology labels
    def checkInt(x):
        try:
            int(x)
            return True
        except ValueError:
            return False

    concordances = {}

    for loc in loci:
        newDict = {}
        simDict = {}
        diffDict = {}
        oldPredict = {}
        newPredict = {}
        oldPredFile = base_dir + "/old-predictions/" + loc + ".chile"
        newPreds = pd.read_csv(base_dir + "predictions/" + loc + "_predictions.csv")
        newPreds = newPreds.set_index('allele')
        newPreds = newPreds.to_dict()
        newPredict = newPreds["serology"]
        for nKey in newPredict.keys():
            adjustMe = str(newPredict[nKey])
            adjustMe = adjustMe.replace('[','')
            adjustMe = adjustMe.replace(']','')
            adjustMe = adjustMe.replace("'",'')
            adjustMe = adjustMe.split(' ')
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
            if each not in newPredict.keys():
                next
            else:
                allDict["New Assignment"] = newPredict[each]
                if set(newPredict[each]) != set(oldPredict[each]):
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

        simLen = len(simFrame)
        diffLen = len(diffFrame)
        with open(base_dir + "comparison/" + loc + "_concordance.txt", "w+") as fhandle:
            fhandle.write("HLA-" +loc+ " Similar: " + str(simLen)+'\n')
            fhandle.write("HLA-" +loc+ " Different: " + str(diffLen)+'\n')
            concordance = (simLen / (simLen + diffLen)) * 100
            concordances[loc] = concordance
            fhandle.write("HLA-" +loc+ " Concordance: " + str(concordance) + "%"+'\n')
            if print_all == "yes":
                print("HLA-" +loc+ " Similar: " + str(simLen))
                print("HLA-" +loc+ " Different: " + str(diffLen))
                print("HLA-" +loc+ " Concordance: " + str(concordance) + "%")
    return concordances

# below function directly from Andrew Chang in fast.ai forums (https://forums.fast.ai/t/automated-learning-rate-suggester/44199)
def find_appropriate_lr(model: Learner, lr_diff: int = 50, loss_threshold: float = .05, adjust_value: float = 1, plot: bool = False) -> float:
    # Run the Learning Rate Finder
    model.lr_find()

    # Get loss values and their corresponding gradients, and get lr values
    losses = np.array(model.recorder.losses)
    assert (lr_diff < len(losses))
    loss_grad = np.gradient(losses)
    lrs = model.recorder.lrs

    # Search for index in gradients where loss is lowest before the loss spike
    # Initialize right and left idx using the lr_diff as a spacing unit
    # Set the local min lr as -1 to signify if threshold is too low
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
        plt.plot(len(losses) + l_idx, loss_grad[l_idx], markersize=10, marker='o', color='red')
        plt.ylabel("Loss")
        plt.xlabel("Index of LRs")
        plt.show()

        plt.plot(np.log10(lrs), losses)
        plt.ylabel("Loss")
        plt.xlabel("Log 10 Transform of Learning Rate")
        loss_coord = np.interp(np.log10(lr_to_use), np.log10(lrs), losses)
        plt.plot(np.log10(lr_to_use), loss_coord, markersize=10, marker='o', color='red')
        plt.show()

    return lr_to_use

pre_concord = metrics()

#loci = ['A', 'B', 'C', 'DQB1', 'DRB1']
loci = ['B']
epoch = { "A":60, "B":120, "C":50, "DPB1":100, "DQB1":100, "DRB1":100 }
layer = { "A":[2000, 1500], "B":[5000, 3000, 1000], "C":[150], "DPB1":[1000], "DQB1":[1000], "DRB1":[1000] }


for locus in loci:

    random_seed(50,use_cuda=True)

    AAs = []
    tng_df = pd.read_csv(base_dir + 'old_sets/train/' + locus + '_train.csv')
    tst_df = pd.read_csv(base_dir + 'old_sets/test/' + locus + '_test.csv')
    val_df = pd.read_csv(base_dir + 'old_sets/train/' + locus + '_validation.csv')
    tng_idx = len(tng_df)
    val_len = len(val_df)
    val_idx = tng_idx + val_len
    tst_idx = len(tst_df)
    tbatch = int(tng_idx // 1.5)
    if (tbatch <= 1):
        tbatch = tng_idx // 1
    vbatch = int(val_idx // 1.5)
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
                                .label_from_df(cols=dep_var, label_delim=' ')
                                .add_test(test)
                                .databunch(bs=tbatch, val_bs=vbatch))

    acc_02 = partial(accuracy_thresh, thresh=0.99)
    #acc_m = partial()
    f_score = partial(fbeta, thresh=0.52)

    pre_vote = {}
    avg_pred = {}
    all_models = []
    #weights = torch.ones([data.c]).float().cuda()
    #loss = nn.BCEWithLogitsLoss(pos_weight=weights)
    #learn = tabular_learner(data, opt_func=optim.SGD, layers=layer[locus], metrics=[acc_02, f_score], loss_func=loss)
    learn = tabular_learner(data, opt_func=optim.SGD, layers=layer[locus], metrics=[acc_02, f_score])
    #learn = tabular_learner(data, layers=layer[locus], metrics=[acc_02, f_score])

    lr = 0.5
    
    #lr = find_appropriate_lr(model=learn)

    #learn.recorder.plot(suggestion=True)
    #learn.fit_one_cycle(epoch[locus], lr)
    learn.fit_one_cycle(epoch[locus], max_lr=slice(lr))
    learn.model
    #learn.recorder.plot_losses()

    test_id = list(tst_df['allele'])

    classes = data.classes
    predictions = {}


    category = learn.get_preds(DatasetType.Test)
    preds = category[0].tolist()
    for count in range(0,len(preds)):
        predictions[test_id[count]] = ' '.join([classes[i] for i in range(len(preds[count])) if (preds[count][i] >= 0.52)])
    

    frame = pd.Series(predictions, name='serology')
    frame.index.name = 'allele'
    frame.to_csv(base_dir + 'predictions/' + locus + '_predictions.csv', index=True)


post_concord = metrics()

for loc in loci:
	print(loc + " Concordance:\t\t\t\t" + str(post_concord[loc])[:5] + "%")
	change = post_concord[loc] - pre_concord[loc]
	print("% Change:\t\t\t\t" + str(change)[:5] + "%")