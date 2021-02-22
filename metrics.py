import pandas as pd
import numpy as np

loci = ['A', 'B', 'C', 'DRB1', 'DQB1']
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
    print(x)
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
        oldPredFile = "./old-predictions/" + loc + ".chile"
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
        opseries.to_csv("./old-predictions/"+loc+"_predictions.csv", line_terminator='\n')
    return

# function to measure concordance between old SNNS and new ML models
def concordance(loci=loci):
    for loc in loci:
        comparison = open("./comparison/" + loc + "_compfile.txt", "w+")
        newsies = open("./comparison/" + loc + "_newsies.txt", "w+")
        similarities = open("./comparison/" + loc + "_similar.txt", "w+")
        oldPredict = {}
        newPredict = {}
        oldPredFile = "./old-predictions/" + loc + ".chile"
        newPreds = pd.read_csv("./predictions/" + loc + "_predictions.csv")
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
            print(newPredict[nKey])
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
                print(each)
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

    sum_df.to_csv('./comparison/summary.csv', index=True)
    
    return

def evaluate(loc, p_allele, relser, right, wrong, partial, close, bad, broad_split=broad_split, split_broad=split_broad):
    p_ser = p_allele.serology.replace("'",'').replace('[','').replace(']','').replace('"','').replace(',','')
    p_ser = set(p_ser.split(' '))
    
    if p_allele.allele in relser.index:
        newser = set(relser.loc[p_allele.allele].serology.split(' '))
        if p_ser == newser:
            right.append(p_allele.allele)
        elif p_ser != newser:
            wrong.append(p_allele.allele)
            if any(w in newser for w in p_ser):
                partial.append(p_allele.allele)
            else:
                switch1 = "no"
                for oldval in p_ser:
                    if (oldval in list(broad_split[loc].keys())) or (oldval in list(split_broad[loc].keys())):
                        if oldval in list(broad_split[loc].keys()):
                            if any(x in newser for x in broad_split[loc][oldval]):
                                close.append(p_allele.allele)
                                switch1 = "yes"
                        elif oldval in list(split_broad[loc].keys()):
                            if any(y in newser for y in split_broad[loc][oldval]):
                                close.append(p_allele.allele)
                                switch1 = "yes"
                if switch1 == "no":
                    bad.append(p_allele.allele)

        return right, wrong, partial, close, bad
    
    else:
        
        return right, wrong, partial, close, bad

def met_pct(datalist, right, wrong, partial, close, bad):
    print(len(datalist))
    n_alleles = len(datalist)
    n_r = len(right)
    n_w = len(wrong)
    n_p = len(partial)
    n_c = len(close)
    n_b = len(bad)

    p_r = (n_r / n_alleles) * 100
    p_w = (n_w / n_alleles) * 100
    p_p = (n_p / n_alleles) * 100
    p_c = (n_c / n_alleles) * 100
    p_b = (n_b / n_alleles) * 100

    p_dict = {
        "Right" : p_r,
        "Wrong" : p_w,
        "Partial" : p_p,
        "Close" : p_c,
        "Bad" : p_b,
    }

    return p_dict

def accuracy(loc, dataframe, relser):
    right = []
    wrong = []
    partial = []
    close = []
    bad = []

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
        
        right, wrong, partial, close, bad = evaluate(loc, all, relser, right, wrong, partial, close, bad)

    df = relser[relser.index.isin(dataframe.allele)]

    met_dict = met_pct(df, right, wrong, partial, close, bad)

    return met_dict


def check_acc_all(loci=loci):
    mets = {
        "Old NN" : {},
        "New NN" : {},
        "Random Forest" : {},
    }

    for loc in loci:    
        old_nn_preds = pd.read_csv("./old-predictions/"+loc+"_predictions.csv", dtype=str)
        new_nn_preds = pd.read_csv("./predictions/"+loc+"_predictions.csv", dtype=str)
        new_nn_preds = new_nn_preds[new_nn_preds.allele.isin(old_nn_preds.allele)]
        rf_preds = pd.read_csv("./randomforest/predictions/"+loc+"_predictions.csv", dtype=str)
        rf_preds = rf_preds[rf_preds.allele.isin(rf_preds.allele)]
        relser = pd.read_csv("./ser/"+loc+"_ser.csv", dtype=str)
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
        r1 = mets[opt1][loc]['Right']
        w1 = mets[opt1][loc]['Wrong']
        p1 = mets[opt1][loc]['Partial']
        c1 = mets[opt1][loc]['Close']
        b1 = mets[opt1][loc]['Bad']
        r2 = mets[opt2][loc]['Right']
        w2 = mets[opt2][loc]['Wrong']
        p2 = mets[opt2][loc]['Partial']
        c2 = mets[opt2][loc]['Close']
        b2 = mets[opt2][loc]['Bad']

        l_dict['Right'] = r1-r2
        l_dict['Wrong'] = w1-w2
        l_dict['Partial'] = p1-p2
        l_dict['Close'] = c1-c2
        l_dict['Bad'] = b1-b2
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
    cframe1.to_csv('./comparison/'+cond1+'.csv', index=True)
    comp_dict[cond2] = compare_acc(mets, opt1, opt3)
    cframe2 = pd.DataFrame.from_dict(comp_dict[cond2])
    cframe2.to_csv('./comparison/'+cond2+'.csv', index=True)
    comp_dict[cond3] = compare_acc(mets, opt2, opt3)
    cframe3 = pd.DataFrame.from_dict(comp_dict[cond3])
    cframe3.to_csv('./comparison/'+cond3+'.csv', index=True)

    return comp_dict

mets = check_acc_all()
mframe1 = pd.DataFrame.from_dict(mets['Old NN'])
mframe1.to_csv('./comparison/OldNN_mets.csv', index=True)
mframe2 = pd.DataFrame.from_dict(mets['New NN'])
mframe2.to_csv('./comparison/NewNN_mets.csv', index=True)
mframe3 = pd.DataFrame.from_dict(mets['Random Forest'])
mframe3.to_csv('./comparison/RF_mets.csv', index=True)
comp_dict = compare_acc_all(mets)