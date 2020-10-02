from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

def fix_data(data, loc, iset, ident):
	sero = {}
	for row in data.itertuples():
		print(row)
		sero[row.allele] = row.serology
	data = data.drop('serology', axis=1)
	uniques = []
	for key in sero.keys():
		if (sero[key].find(';') != -1):
			sero[key] = sero[key].replace('a','')
			sero[key] = sero[key].split(';')
		else:
			sero[key] = sero[key].replace('a','')
			sero[key] = list(sero[key])

		for x in sero[key]:
			if (x not in uniques):
				uniques.append(x)
			else:
				next

	uniques = list(map(int, uniques))
	uniques.sort()
	uniques = list(map(str, uniques))
	for y in uniques:
		data[y] = 0

	one_sero = {}
	for key in sero.keys():
		one_sero[key] = { some_key : ("1" if (some_key in sero[key]) else "0")
		                  for some_key in uniques }

	print(one_sero)
	one_df = pd.DataFrame.from_dict(one_sero)
	print(one_df)
	one_df = one_df.transpose()
	one_df.index.name = "allele"
	print(one_df)
	data = data.set_index('allele')
	data.update(one_df, overwrite=True)
	data.to_csv('./randfor/' + iset + '/' + loc + '_' + ident + '.csv' )
	return data

RSEED = 50

loci = ["A", "B", "C", "DQB1", "DRB1", "DPB1"]
for loc in loci:
	features = pd.read_csv("./training/" + loc + "_train.csv")
	features = fix_data(features,loc,iset='training',ident='train')
	vfeatures = pd.read_csv("./training/" + loc + "_validation.csv")
	vfeatures = fix_data(vfeatures,loc,iset='training',ident='validation')
	other = pd.read_csv("./testing/" + loc + "_test.csv")
	other = other.drop('serology', axis=1)
	other.to_csv('./randfor/testing/' + loc + '_test.csv')
	# features = features.append(vfeatures)
	# labels = np.array(features['serology'])
	# features = features.drop('serology', axis=1)
	# indices = features["allele"]
	# features = features.drop('allele', axis=1)
	# feature_list = list(features.columns)
	# features = np.array(features)
	#
	# train_features, test_features, train_labels, test_labels = \
	# 	train_test_split(features, labels, test_size=0.25,
	# 	                 random_state=RSEED)
	#
	# newlabels = []
	# drops = []
	# for label in train_labels:
	# 	if (label.find(';') != -1):
	# 		alter = label.split(';')
	# 		drops.append(label)
	# 		newlabels.append(alter)
	# for each in drops:
	# 	train_labels = train_labels[train_labels != each]
	# np.append(train_labels, newlabels)
	#
	# rf = RandomForestClassifier(n_estimators=1000, random_state=RSEED,
	#                             max_features='sqrt', n_jobs=-1, verbose=1)
	# rf.fit(train_features, train_labels)
	# predictions = rf.predict(test_features)
	# print(classification_report(test_labels,predictions))
