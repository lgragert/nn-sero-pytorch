import pandas as pd
import numpy as np

def main(print_all='no'):
    loci = ['A', 'B', 'C', 'DPB1', 'DQB1', 'DRB1']
    #loci = ['A']

    base_dir = './'
    RF_dir = './randomforest/'

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
        oldPredFile = base_dir + "old-predictions/" + loc + ".chile"
        newPreds = pd.read_csv(RF_dir + "predictions/" + loc + "_predictions.csv")
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
        diffFrame.to_csv(RF_dir + "comparison/" + loc + "_compfile.csv", index=False)
        simFrame = pd.DataFrame.from_dict(simDict)
        simFrame = simFrame.transpose()
        simFrame.to_csv(RF_dir + "comparison/" + loc + "_similar.csv", index=False)
        

        for allele in newPredict.keys():
            allDict = {}
            allDict["Allele"] = allele
            allDict["Serologic Assignment"] = newPredict[allele]
            if allele not in oldPredict.keys():
                newDict[allele] = allDict
        newFrame = pd.DataFrame.from_dict(simDict)
        newFrame = newFrame.transpose()
        newFrame.to_csv(RF_dir + "comparison/" + loc + "_newsies.csv", index=False)

        simLen = len(simFrame)
        diffLen = len(diffFrame)
        with open(RF_dir + "comparison/" + loc + "_concordance.txt", "w+") as fhandle:
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