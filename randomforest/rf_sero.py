import pandas as pd
import numpy as np
import sys
import metrics_RF as metrics
from tqdm import tqdm
from sklearn.ensemble import RandomForestClassifier
from sklearn.multioutput import MultiOutputClassifier

np.set_printoptions(threshold=sys.maxsize)

def one_hot_decode(df):
	df['serology']=''

	for col in df.columns:
		df.loc[df[col]==1,'serology'] = df['serology']+col+';'

	return df

def fix_data(uniques, data, loc, iset, ident):
	sero = {}
	for row in data.itertuples():
		sero[row.allele] = str(row.serology)
	data = data.drop('serology', axis=1)
	for key in sero.keys():
		# not applicable for old_sets train/test
		'''
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
	data.to_csv('./randomforest/randfor/'+iset+'/'+loc+'_'+ident+'.csv', index=True)
	return data, uniques

RSEED = 5

pre_concord = metrics.main()

loci = ["A", "B", "C", "DQB1", "DRB1", "DPB1"]
print("Predicting...")
for loc in tqdm(loci):
	uniques = []
	#print(loc)
	features = pd.read_csv("./randomforest/training/" + loc + "_train.csv")
	features, sers = fix_data(uniques, features,loc,iset='training',ident='train')
	vfeatures = pd.read_csv("./randomforest/training/" + loc + "_validation.csv")
	vfeatures, vsers = fix_data(uniques, vfeatures,loc,iset='training',ident='validation')
	test = pd.read_csv("./randomforest/testing/" + loc + "_test.csv")
	test = test.drop('serology', axis=1)
	test.to_csv('./randomforest/randfor/testing/'+loc+'_test.csv', index=True)

	features = features.append(vfeatures)
	labels = np.array(features[sers])
	features = features.drop(sers, axis=1)
	features = features.reset_index()
	indices = features["allele"]
	indices = list(indices)
	features = features.drop('allele', axis=1)
	feature_list = list(features.columns)
	n_features = len(feature_list)
	maxfeat = int(n_features)
	features = np.array(features)
	labels[labels!=labels]='0'
	features[features!=features]='0'
	features = features.astype(np.int)
	labels = labels.astype(np.int)

	test_idcs = test['allele']
	test = test.drop('allele', axis=1)
	test_list = list(test.columns)
	test = np.array(test)
	test[test!=test]='0'
	test = test.astype(np.int)

	

	forest = RandomForestClassifier(n_estimators=5000, random_state=RSEED, max_features=maxfeat, n_jobs=-1)
	multi_target_forest = MultiOutputClassifier(forest, n_jobs=-1)
	multi_target_forest.fit(features,labels)
	predictions = multi_target_forest.predict(test)

	ind_labels = [str(x) for x in sers]
	preds_output = pd.DataFrame(predictions, index=test_idcs, columns=ind_labels)
	preds_output = one_hot_decode(preds_output)
	preds_output = preds_output.drop(ind_labels, axis=1)
	preds_output.index.name = 'allele'
	preds_output = preds_output.apply(lambda x: str((x['serology'].split(';'))[:-1]), result_type='broadcast', axis=1)
	#print(preds_output)
	preds_output.to_csv('./randomforest/predictions/'+loc+'_predictions.csv', index=True)


post_concord = metrics.main()

print("Done.")

for loc in loci:
	print(loc + " Concordance:\t\t\t\t" + str(post_concord[loc])[:5] + "%")
	change = post_concord[loc] - pre_concord[loc]
	print("% Change:\t\t\t\t" + str(change)[:5] + "%")
