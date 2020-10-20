import pandas as pd
import numpy as np

loci = ['A', 'B', 'C', 'DPB1', 'DRB1', 'DQB1']



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