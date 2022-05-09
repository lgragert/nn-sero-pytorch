# %% [markdown]
# # Random Forest for Serology Prediction

# %% [markdown]
# ## Setup Steps

from pathlib import Path
root_dir = './'
base_dir = root_dir + 'randomforest/'
path = Path(base_dir)
NN_dir = './'

# %% [markdown]
# ## Random Forest Modeling

# %% [markdown]
# ### Prior Concordance Assessment

# %%
import pandas as pd
import numpy as np
import sys
import math
import lime
import lime.lime_tabular
from tqdm import tqdm
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from collections import defaultdict
from scipy.stats import spearmanr
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from sklearn.inspection import permutation_importance
#from sklearn.model_selection import train_test_split


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
        oldPredFile = Path(NN_dir + "old-predictions/" + loc + ".chile")
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

        if loc == 'C':
            skipc = ['C*01', 'C*02', 'C*03', 'C*04', 'C*05', 'C*06', 'C*07', 'C*08']
            oldPredict = {k:v for k,v in oldPredict.items() if k[:4] in skipc}
            newPredict = {k:v for k,v in newPredict.items() if k[:4] in skipc}


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
            fhandle.write("HLA-" +loc+ " Similar: " + str(simLen))
            fhandle.write("HLA-" +loc+ " Different: " + str(diffLen))
            concordance = (simLen / (simLen + diffLen)) * 100
            concordances[loc] = concordance
            fhandle.write("HLA-" +loc+ " Concordance: " + str(concordance) + "%")
            if print_all == "yes":
                print("HLA-" +loc+ " Similar: " + str(simLen))
                print("HLA-" +loc+ " Different: " + str(diffLen))
                print("HLA-" +loc+ " Concordance: " + str(concordance) + "%")
    return concordances

#main(print_all="yes")

# %% [markdown]
# ### Additional Data

# %%
# All data here from Sigma Aldrich
# https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html

mol_wghts = {
    'A': 89.10,
    'R': 174.20,
    'N': 132.12,
    'D': 133.11,
    'C': 121.16,
    'E': 147.13,
    'Q': 146.15,
    'G': 75.07,
    'H': 155.16,
    #'O': 131.13,
    'I': 131.18,
    'L': 131.18,
    'K': 146.19,
    'M': 149.21,
    'F': 165.19,
    'P': 115.13,
    #'U': 139.11,
    'S': 105.09,
    'T': 119.12,
    'W': 204.23,
    'Y': 181.19,
    'V': 117.15,
    '-': 0,
    'X': 0,
}

pKa = {
    'A': 2.34,
    'R': 2.17,
    'N': 2.02,
    'D': 1.88,
    'C': 1.96,
    'E': 2.19,
    'Q': 2.17,
    'G': 2.34,
    'H': 1.82,
    #'O': 1.82,
    'I': 2.36,
    'L': 2.36,
    'K': 2.18,
    'M': 2.28,
    'F': 1.83,
    'P': 1.99,
    #'U': 0,
    'S': 2.21,
    'T': 2.09,
    'W': 2.83,
    'Y': 2.20,
    'V': 2.32,
    '-': 0,
    'X': 0,
}

pKb = {
    'A': 9.69,
    'R': 9.04,
    'N': 8.80,
    'D': 9.60,
    'C': 10.28,
    'E': 9.67,
    'Q': 9.13,
    'G': 9.60,
    'H': 9.17,
    #'O': 9.65,
    'I': 9.60,
    'L': 9.60,
    'K': 8.95,
    'M': 9.21,
    'F': 9.13,
    'P': 10.60,
    #'U': 0,
    'S': 9.15,
    'T': 9.10,
    'W': 9.39,
    'Y': 9.11,
    'V': 9.62,
    '-': 0,
    'X': 0,
}

pKx = {
    'A': 0,
    'R': 12.48,
    'N': 0,
    'D': 3.65,
    'C': 8.18,
    'E': 4.25,
    'Q': 0,
    'G': 0,
    'H': 6.00,
    #'O': 0,
    'I': 0,
    'L': 0,
    'K': 10.53,
    'M': 0,
    'F': 0,
    'P': 0,
    #'U': 0,
    'S': 0,
    'T': 0,
    'W': 0,
    'Y': 10.07,
    'V': 0,
    '-': 0,
    'X': 0,
}

pI = {
    'A': 6.00,
    'R': 10.76,
    'N': 5.41,
    'D': 2.77,
    'C': 5.07,
    'E': 3.22,
    'Q': 5.65,
    'G': 5.97,
    'H': 7.59,
    #'O': 0,
    'I': 6.02,
    'L': 5.98,
    'K': 9.74,
    'M': 5.74,
    'F': 5.48,
    'P': 6.30,
    #'U': 5.68,
    'S': 5.68,
    'T': 5.60,
    'W': 5.89,
    'Y': 5.66,
    'V': 5.96,
    '-': 0,
    'X': 0,
}

# Hydrophobicity Index at pH 2
HI2 = {
    'A': 47,
    'R': -26,
    'N': -18,
    'D': -18,
    'C': 52,
    'E': 8,
    'Q': -18,
    'G': 0,
    'H': -42,
    #'O': 0,
    'I': 100,
    'L': 100,
    'K': -37,
    'M': 74,
    'F': 92,
    'P': -46,
    #'U': 0,
    'S': -7,
    'T': 13,
    'W': 84,
    'Y': 49,
    'V': 79,
    '-': 0,
    'X': 0,
}

# Hydrophobicity index at pH 7
HI7 = {
    'A': 41,
    'R': -14,
    'N': -28,
    'D': -55,
    'C': 49,
    'E': -31,
    'Q': -10,
    'G': 0,
    'H': 8,
    #'O': 0,
    'I': 99,
    'L': 97,
    'K': -23,
    'M': 74,
    'F': 100,
    'P': -46, # SA used pH 2
    #'U': 0,
    'S': -5,
    'T': 13,
    'W': 97,
    'Y': 63,
    'V': 76,
    '-': 0,
    'X': 0,
}


# %% [markdown]
# ### Data Preprocessing

# %%
np.set_printoptions(threshold=sys.maxsize)

def one_hot_decode(df):
	df['serology']=''

	for col in df.columns:
		df.loc[df[col]==1,'serology'] = df['serology']+col+';'

	return df

def fix_data(uniques, data, loc, iset, ident):
    sero = {}
    for row in data.itertuples(name='Pandas'):
        sero[row.allele] = str(row.serology)
        #sero[row[1]] = str(row[-1])

    data = data.drop('serology', axis=1)

    for key in sero.keys():
        '''
    	# not applicable for old_sets train/test
        if (sero[key].find(';') != -1):
            sero[key] = sero[key].replace('a','')
            sero[key] = sero[key].split(';')
        else:
            sero[key] = sero[key].replace('a','')
            sero[key] = [sero[key]]
        '''

        #for old_sets train/test
        sero[key] = sero[key].split(' ')

        for x in sero[key]:
            if (x not in uniques):
                uniques.append(x)
            else:
                continue

    uniques = list(map(int, uniques))
    uniques.sort()
    uniques = list(map(str, uniques))

    for y in uniques:
        data[y] = 0

    one_sero = {}
    for key in sero.keys():
        one_sero[key] = { some_key : ("1" if (some_key in sero[key]) else "0")
  	  	                    for some_key in uniques }
    one_df = pd.DataFrame.from_dict(one_sero)
    one_df = one_df.transpose()
    one_df.index.name = "allele"
    data = data.set_index('allele')
    data.update(one_df, overwrite=True)
    data.to_csv(base_dir + 'randfor/'+iset+'/'+loc+'_'+ident+'.csv', index=True)
    return data, uniques


def add_data(df):
    OHcols = list(df.columns)
    OHcols.remove('serology')
    OHcols.remove('allele')

    for col in OHcols:
        MWname = 'MW'+str(col)
        df[MWname] = df.apply(lambda row: mol_wghts[row[col]], axis=1)
        pKaname = 'pKa'+str(col)
        df[pKaname] = df.apply(lambda row: pKa[row[col]], axis=1)
        pKbname = 'pKb'+str(col)
        df[pKbname] = df.apply(lambda row: pKb[row[col]], axis=1)
        pKxname = 'pKx'+str(col)
        df[pKxname] = df.apply(lambda row: pKx[row[col]], axis=1)
        pIname = 'pI'+str(col)
        df[pIname] = df.apply(lambda row: pI[row[col]], axis=1)
        HI2name = 'HI2_'+str(col)
        df[HI2name] = df.apply(lambda row: HI2[row[col]], axis=1)
        HI7name = 'HI7_'+str(col)
        df[HI7name] = df.apply(lambda row: HI7[row[col]], axis=1)



    df = df.drop(columns=OHcols, axis=1)
    df = df.drop('serology', axis=1)

    return df


# %% [markdown]
# ### Current Best Parameters for Random Forest

# %%

RSEED = 0

pre_concord = metrics()

loci = ["A", "B", "C", "DQB1", "DRB1"]
#loci = ['DQB1']
nest = {
    'A': 40, #10
    'B': 25, #199
    'C': 5, #14
    'DQB1': 64, #10
    'DRB1': 16 #250
}

modfeat = {
    'A':  'auto', #auto
    'B': 1, #auto
    'C': 'log2', #log2
    'DQB1': 1, #log2
    'DRB1': 'log2' #log2
}


numfeat = {
    'A':  100, #False
    'B': 200, #False
    'C': 200, #False
    'DQB1': 100, #False
    'DRB1': False #log2
}

boot = {
    'A': False,
    'B': False,
    'C': False,
    'DQB1': False,
    'DRB1': False
}

#boot = {
#    'A': True,
#    'B': True,
#    'C': True,
#    'DQB1': True,
#    'DRB1': True
#}

criteria = {
    'A': 'gini',
    'B': 'gini',
    'C': 'gini',
    'DQB1': 'gini',
    'DRB1': 'gini'
}

oob = {
    'A': False,
    'B': False,
    'C': False,
    'DQB1': False,
    'DRB1': False
}

#oob = {
#    'A': True,
#    'B': True,
#    'C': True,
#    'DQB1': True,
#    'DRB1': True
#}

# %% [markdown]
# ### Random Forest Classifiers

# %%
import warnings
warnings.filterwarnings('ignore')
RSEED = 0

pre_concord = metrics()

print("Predicting...")
for loc in tqdm(loci):
    uniques = []
    print(loc+'...', end='')
    features = pd.read_csv(base_dir + "OHtraining/" + loc + "_train.csv")
    vfeatures = pd.read_csv(base_dir + "OHtraining/" + loc + "_validation.csv")
    test = pd.read_csv(base_dir + "OHtesting/" + loc + "_test.csv")
    features, sers = fix_data(uniques, features,loc,iset='training',ident='train')
    vfeatures, vsers = fix_data(uniques, vfeatures,loc,iset='training',ident='validation')

    test = test.drop('serology', axis=1)
    test.to_csv(base_dir + 'randfor/testing/'+loc+'_test.csv', index=True)


    features = features.append(vfeatures)
    # had to change following two lines from sers to vsers to account for additional validation data
    labels = np.array(features[vsers])
    features = features.drop(vsers, axis=1)
    features = features.reset_index()
    indices = features["allele"]
    indices = list(indices)
    features = features.drop('allele', axis=1)
    feature_list = list(features.columns)
    n_features = len(feature_list)
    #maxfeat = int(math.sqrt(n_features))

    all_features = features
    features = np.array(features)
    labels[labels!=labels]='0'
    features[features!=features]='0'
    features = features.astype(int)
    labels = labels.astype(int)

    test_idcs = test['allele']
    test = test.drop('allele', axis=1)
    #print(test.head(100))
    all_test = test
    test_list = list(test.columns)
    test = np.array(test)
    test[test!=test]='0'
    test = test.astype(int)
    ind_labels = [str(x) for x in sers]


    all_predictions = []
    for idx in range(0,len(ind_labels)):
        ilabels = labels[:,idx]
        forest = RandomForestClassifier(n_estimators=nest[loc],
                                        criterion=criteria[loc],
                                        bootstrap=boot[loc],
                                        oob_score=oob[loc],
                                        max_features=modfeat[loc],
                                        random_state=RSEED,
                                        n_jobs=-1)
        forest.fit(features,ilabels)
        #predictions = forest.predict(test)
        threshold = 0.42

        #predictions = forest.predict_proba(test)
        #predictions[:,0] = (predictions[:,0] < threshold).astype('int')
        #predictions = (predictions[:,1] >= threshold).astype('int')
        #all_predictions.append(predictions)


        # Feature/Permutation Importance
        #feat_indices = np.argsort(forest.feature_importances_)[::-1]
        if numfeat[loc] != False:
            tree_importance_sorted_idx = np.argsort(forest.feature_importances_)
            tree_importance_sorted_idx = tree_importance_sorted_idx[-numfeat[loc]:]
            less_features = np.array(feature_list)[tree_importance_sorted_idx]
            new_features = all_features[less_features]
            new_features = np.array(new_features)
            new_features[new_features!=new_features]='0'
            new_features = new_features.astype(int)
            less_test = np.array(test_list)[tree_importance_sorted_idx]
            new_test = all_test[less_test]
            new_test = np.array(new_test)
            new_test[new_test!=new_test]='0'
            new_test = new_test.astype(int)
            new_forest = RandomForestClassifier(n_estimators=nest[loc],
                                            criterion=criteria[loc],
                                            bootstrap=boot[loc],
                                            oob_score=oob[loc],
                                            max_features=modfeat[loc],
                                            random_state=RSEED,
                                            n_jobs=-1)
            new_forest.fit(new_features,ilabels)
            #threshold = 0.42

            predictions = new_forest.predict_proba(new_test)
            predictions[:,0] = (predictions[:,0] < threshold).astype('int')
            predictions = (predictions[:,1] >= threshold).astype('int')
            all_predictions.append(predictions)
        else:
            predictions = forest.predict_proba(test)
            predictions[:,0] = (predictions[:,0] < threshold).astype('int')
            predictions = (predictions[:,1] >= threshold).astype('int')
            all_predictions.append(predictions)


    all_predictions = np.asarray(all_predictions)
    all_predictions = np.transpose(all_predictions)

    explainer = lime.lime_tabular.LimeTabularExplainer(features,feature_names=feature_list,class_names=ind_labels,kernel_width=5)
    for rowexp in range(0,1):
      exp = explainer.explain_instance(test[rowexp], forest.predict_proba, num_features=n_features)
      exp.save_to_file('{}_{}.lime.html'.format(loc, str(rowexp)), show_table=True)

    #preds_output = pd.DataFrame(predictions, index=test_idcs, columns=ind_labels)
    preds_output = pd.DataFrame(all_predictions, index=test_idcs, columns=ind_labels)
    preds_output = one_hot_decode(preds_output)
    preds_output = preds_output.drop(ind_labels, axis=1)
    preds_output.index.name = 'allele'
    preds_output = preds_output.apply(lambda x: str((x['serology'].split(';'))[:-1]), result_type='broadcast', axis=1)
    preds_output.to_csv(base_dir + 'predictions/'+loc+'_predictions.csv', index=True)

print("Done.")

# %% [markdown]
# ### Concordance Checker

# %%
post_concord = metrics()

for loc in loci:
	print(loc + " Concordance:\t\t\t\t" + str(post_concord[loc])[:5] + "%")
	change = post_concord[loc] - pre_concord[loc]
	print("% Change:\t\t\t\t" + str(change)[:5] + "%")

# %% [markdown]
# ### Accuracy Checker

# %%
import pandas as pd
import numpy as np

loci = ['A', 'B', 'C', 'DQB1', 'DRB1']
summary = {}

# dict of dicts to store splits of broad specificities
broad_split = {
    "A" : {
        "9" : ["23", "24"],
        "10" : ["25", "26", "34", "66"],
        "19" : ["29", "30", "31", "32", "33", "74"],
        "28" : ["68", "69"]
    },
    "B" : {
        "5" : ["51", "52"],
        "12" : ["44", "45"],
        "14" : ["64", "65"],
        "15" : ["62", "63", "75", "76", "77"],
        "16" : ["38", "39"],
        "17" : ["57", "58"],
        "21" : ["49", "50"],
        "22" : ["54", "55", "56"],
        "40" : ["60", "61"],
        "70" : ["71", "72"]
    },
    "C" : {
        "3" : ["9", "10"]
    },
    "DQB1" : {
        "1" : ["5", "6"],
        "3" : ["7", "8", "9"]
    },
    "DRB1" : {
        "2" : ["15", "16"],
        "3" : ["17", "18"],
        "5" : ["11", "12"],
        "6" : ["13", "14"]
    }
}

# dict of dict to store broad specificity for each split
split_broad = {}

for alphakey in broad_split.keys():
    sb = {}

    for betakey in broad_split[alphakey].keys():
        for value in broad_split[alphakey][betakey]:
            sb[value] = betakey

    split_broad[alphakey] = sb

# function to check if value can be an integer - to eliminate excess characters from serology labels
def checkInt(x):
    try:
        int(x)
        return True
    except ValueError:
        return False

# function to eliminate any serological assignments with under a 95% likelihood
def chance(x, line):
    if (line[x].find("%") != -1):
        x = float(line[x][:-1])
        if 51 <= x:
            test = True
        else:
            test = False
    else:
        test = False
    return test

# function to generate dataframes to contain SNNS predictions
def SNNS_preds(loci=loci):
    for loc in loci:
        oldPredict = {}
        oldPredFile = NN_dir + "old-predictions/" + loc + ".chile"
        with open(oldPredFile, "r") as handle:
            for line in handle:
                if line.find('%') != -1:
                    line = line.split()
                    if line != []:
                        line[:] = [x for x in line if x != '[100.00%]']
                        allele = loc + "*" + str(line[0][:-1])
                        oldPredict[allele] = ' '.join(line[1:])
                else:
                    next

        opseries = pd.Series(oldPredict, name="serology")
        opseries.index.name = "allele"
        opseries.reset_index()
        opseries.to_csv(NN_dir+"old-predictions/"+loc+"_predictions.csv", line_terminator='\n')
    return

# function to measure concordance between old SNNS and new ML models
def concordance(loci=loci):
    for loc in loci:
        comparison = open(NN_dir+"comparison/" + loc + "_compfile.txt", "w+")
        newsies = open(NN_dir+"comparison/" + loc + "_newsies.txt", "w+")
        similarities = open(NN_dir+"comparison/" + loc + "_similar.txt", "w+")
        oldPredict = {}
        newPredict = {}
        oldPredFile = NN_dir+"old-predictions/" + loc + ".chile"
        newPreds = pd.read_csv(NN_dir+"predictions/" + loc + "_predictions.csv")
        newPreds = newPreds.set_index('allele')
        newPreds = newPreds.to_dict()
        newPredict = newPreds["serology"]
        for nKey in newPredict.keys():
            adjustMe = str(newPredict[nKey])
            adjustMe = adjustMe.replace('[','')
            adjustMe = adjustMe.replace(']','')
            adjustMe = adjustMe.replace('a','')
            adjustMe = adjustMe.replace("'",'')
            adjustMe = adjustMe.split(' ')
            newPredict[nKey] = [x.strip('a') for x in adjustMe if checkInt(x)]
        with open(oldPredFile, "r") as handle:
            for line in handle:
                if line.find('%') == -1:
                    next
                else:
                    line = handle.readline()
                    line = line.split()
                    if line == []:
                        next
                    else:
                        line[:] = [x for x in line if x != '[100.00%]']
                        allele = loc + "*" + str(line[0][:-1])
                        oldPredict[allele] = line[1:]

        for each in oldPredict.keys():
            if each not in newPredict.keys():
                next
            elif set(newPredict[each]) != set(oldPredict[each]):
                comparison.write("Different: " + str(each) + "\n")
                comparison.write("Old Serologic Assignment: " + str(oldPredict[each]) + "\n")
                comparison.write("New Serologic Assignment: " + str(newPredict[each]) + "\n")
            elif set(newPredict[each]) == set(oldPredict[each]):
                similarities.write("Same: " + str(each) + "\n")
                similarities.write("Old Serologic Assignment: " + str(oldPredict[each]) + "\n")
                similarities.write("New Serologic Assignment: " + str(newPredict[each]) + "\n")
        comparison.close()
        similarities.close()

        for allele in newPredict.keys():
            if allele not in oldPredict.keys():
                newsies.write("NEW: " + str(allele) + "\n")
                newsies.write("Serologic Assignment: " + str(newPredict[allele]) + "\n")
        newsies.close()

    return

def summary_table(loci=loci, summary=summary):
    for locus in loci:
        summary_data = {}
        trn_set = pd.read_csv('training/' + locus + '_train.csv')
        val_set = pd.read_csv('training/' + locus + '_validation.csv')
        tst_set = pd.read_csv('testing/' + locus + '_test.csv')

        old_trn_set = pd.read_csv('old_sets/train/' + locus + '_train.csv')
        old_val_set = pd.read_csv('old_sets/train/' + locus + '_validation.csv')
        old_tst_set = pd.read_csv('old_sets/test/' + locus + '_test.csv')

        trnlen = float(len(trn_set))
        vallen = float(len(val_set))
        tstlen = float(len(tst_set))
        polyAA = float(len(trn_set.iloc[0])) - 1
        oldtrnlen = float(len(old_trn_set))
        oldvallen = float(len(old_val_set))
        oldtstlen = float(len(old_tst_set))
        oldpolyAA = float(len(old_trn_set.iloc[0])) - 1

        summary_data['Number of Training Alleles'] = trnlen
        summary_data['R-SNNS Number of Training Alleles'] = oldtrnlen
        summary_data['Difference in Training Set'] = trnlen - oldtrnlen
        summary_data['Percent (%) Growth in Training Set'] = ((trnlen - oldtrnlen)/oldtrnlen) * 100
        summary_data['Number of Validation Alleles'] = vallen
        summary_data['R-SNNS Number of Validation Alleles'] = oldvallen
        summary_data['Difference in Validation Set'] = vallen - oldvallen
        summary_data['Percent (%) Growth in Validation Set'] = ((vallen - oldvallen)/oldvallen) * 100
        summary_data['Number of Testing Alleles'] = tstlen
        summary_data['R-SNNS Number of Testing Alleles'] = oldtstlen
        summary_data['Difference in Testing Set'] = tstlen - oldtstlen
        summary_data['Percent (%) Growth in Testing Set'] = ((tstlen - oldtstlen)/oldtstlen) * 100
        summary_data['Number of Polymorphisms'] = polyAA
        summary_data['R-SNNS Number of Polymorphisms'] = oldpolyAA
        summary_data['Difference in Polymorphisms'] = polyAA - oldpolyAA
        summary_data['Percent (%) Growth in Polymorphisms'] = ((polyAA - oldpolyAA)/oldpolyAA) * 100
        summary[locus] = summary_data

    sum_df = pd.DataFrame(data=summary)

    sum_df.to_csv(NN_dir+'comparison/summary.csv', index=True)

    return

def evaluate(loc, p_allele, relser, right, wrong, partial, bad, broad_split=broad_split, split_broad=split_broad):
    p_ser = str(p_allele.serology).replace("'",'').replace('[','').replace(']','').replace('"','').replace(',','')
    p_ser = set(p_ser.split(' '))

    if p_allele.allele in relser.index:
        newser = set(relser.loc[p_allele.allele].serology.split(' '))
        if p_ser == newser:
            right.append(p_allele.allele)
        elif p_ser != newser:
            wrong.append(p_allele.allele)
            if any([w in newser for w in p_ser]):
                partial.append(p_allele.allele)
            else:
                bad.append(p_allele.allele)

        return right, wrong, partial, bad

    else:

        return right, wrong, partial, bad

def met_pct(datalist, right, wrong, partial, bad):
    n_alleles = len(datalist)
    #print(n_alleles)
    n_r = len(right)
    n_w = len(wrong)
    n_p = len(partial)
    n_b = len(bad)

    p_r = (n_r / n_alleles) * 100
    p_w = (n_w / n_alleles) * 100
    p_p = (n_p / n_alleles) * 100
    p_b = (n_b / n_alleles) * 100

    p_dict = {
        "All Calls Correct" : p_r,
        "Incorrect" : p_w,
        "At Least One Correct Call" : p_p,
        "All Calls Incorrect" : p_b,
    }

    return p_dict

def accuracy(loc, dataframe, relser):
    right = []
    wrong = []
    partial = []
    bad = []

    #print(loc)

    for all in dataframe.iloc:
        # FIXME - A*23:19Q does not appear in rel_dna_ser (A*23:19N instead)
        # FIXME - B*07:44 does not appear in rel_dna_ser (B*07:44N instead)
        # FIXME - B*08:06 does not appear in rel_dna_ser at all
        # FIXME - B*49:15 does not appear in rel_dna_ser at all
        # FIXME - C*03:23 does not appear in rel_dna_ser (C*03:23N instead)
        # FIXME - C*03:99 does not appear in rel_dna_ser at all
        # FIXME - C*05:02 does not appear in rel_dna_ser at all
        # FIXME - C*07:226 does not appear in rel_dna_ser (C*07:226Q instead)
        if all.allele in ["A*23:19Q", "B*07:44", "B*08:06", "B*49:15", "C*03:23", "C*03:99", "C*05:02", "C*07:226"]:
            continue

        right, wrong, partial, bad = evaluate(loc, all, relser, right, wrong, partial, bad)

    df = relser[relser.index.isin(dataframe.allele)]

    met_dict = met_pct(df, right, wrong, partial, bad)

    return met_dict


def check_acc_all(loci=loci):
    mets = {
        "Old NN" : {},
        "New NN" : {},
        "Random Forest" : {},
    }

    for loc in loci:
        old_nn_preds = pd.read_csv(NN_dir+"old-predictions/"+loc+"_predictions.csv", dtype=str)
        new_nn_preds = pd.read_csv(NN_dir+"predictions/"+loc+"_predictions.csv", dtype=str)
        new_nn_preds = new_nn_preds[new_nn_preds.allele.isin(old_nn_preds.allele)]
        rf_preds = pd.read_csv(NN_dir+"randomforest/predictions/"+loc+"_predictions.csv", dtype=str)
        rf_preds = rf_preds[rf_preds.allele.isin(old_nn_preds.allele)]
        # added this filter to make sure all comparisons are identical
        old_nn_preds = old_nn_preds[old_nn_preds.allele.isin(rf_preds.allele)]
        #print('{}:\t\t{} alleles'.format(loc, str(len(old_nn_preds))))
        relser = pd.read_csv(NN_dir+"ser/"+loc+"_ser.csv", dtype=str)
        relser = relser.set_index('allele')
        relser = relser.dropna()

        mets["Old NN"][loc] = accuracy(loc, old_nn_preds, relser)
        mets["New NN"][loc] = accuracy(loc, new_nn_preds, relser)
        mets["Random Forest"][loc] = accuracy(loc, rf_preds, relser)

    return mets

def compare_acc(mets, opt1, opt2, loci=loci):
    c_dict = {}
    for loc in loci:
        l_dict = {}
        r1 = mets[opt1][loc]['All Calls Correct']
        w1 = mets[opt1][loc]['Incorrect']
        p1 = mets[opt1][loc]['At Least One Correct Call']
        b1 = mets[opt1][loc]['All Calls Incorrect']
        r2 = mets[opt2][loc]['All Calls Correct']
        w2 = mets[opt2][loc]['Incorrect']
        p2 = mets[opt2][loc]['At Least One Correct Call']
        b2 = mets[opt2][loc]['All Calls Incorrect']

        l_dict['All Calls Correct'] = r2-r1
        l_dict['Incorrect'] = w2-w1
        l_dict['At Least One Correct Call'] = p2-p1
        l_dict['All Calls Incorrect'] = b2-b1
        c_dict[loc] = l_dict

    return c_dict

def compare_acc_all(mets):
    cond1 = "New_vs_Old_NN"
    cond2 = "Random_Forest_vs_Old_NN"
    cond3 = "Random_Forest_vs_New_NN"

    opt1 = "Old NN"
    opt2 = "New NN"
    opt3 = "Random Forest"

    comp_dict = {
        cond1 : {},
        cond2 : {},
        cond3 : {},
    }

    comp_dict[cond1] = compare_acc(mets, opt1, opt2)
    cframe1 = pd.DataFrame.from_dict(comp_dict[cond1])
    cframe1.to_csv(NN_dir+'comparison/'+cond1+'.csv', index=True)
    comp_dict[cond2] = compare_acc(mets, opt1, opt3)
    cframe2 = pd.DataFrame.from_dict(comp_dict[cond2])
    cframe2.to_csv(NN_dir+'comparison/'+cond2+'.csv', index=True)
    comp_dict[cond3] = compare_acc(mets, opt2, opt3)
    cframe3 = pd.DataFrame.from_dict(comp_dict[cond3])
    cframe3.to_csv(NN_dir+'comparison/'+cond3+'.csv', index=True)

    return comp_dict

concordance()
mets = check_acc_all()
mframe1 = pd.DataFrame.from_dict(mets['Old NN'])
mframe1.to_csv(NN_dir+'comparison/OldNN_mets.csv', index=True)
mframe2 = pd.DataFrame.from_dict(mets['New NN'])
mframe2.to_csv(NN_dir+'comparison/NewNN_mets.csv', index=True)
mframe3 = pd.DataFrame.from_dict(mets['Random Forest'])
mframe3.to_csv(NN_dir+'comparison/RF_mets.csv', index=True)
print("Random Forest:")
print(mframe3.to_string())
print('\n')
print("RSNNS:")
print(mframe1.to_string())
comp_dict = compare_acc_all(mets)

