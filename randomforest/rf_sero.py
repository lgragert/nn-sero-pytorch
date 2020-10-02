from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np
from sklearn.multioutput import MultiOutputClassifier
from sklearn.model_selection import train_test_split


loci = ["A", "B", "C", "DQB1", "DRB1", "DPB1"]
for loc in loci:
	features = pd.read_csv("./training/" + loc + "_train.csv")
	vfeatures = pd.read_csv("./training/" + loc + "_validation.csv")
	features = features.append(vfeatures)
	labels = np.array(features['serology'])
	features = features.drop('serology', axis=1)
	feature_list = list(features.columns)
	features = np.array(features)

	train_features, test_features, train_labels, test_labels = \
		train_test_split(features, labels, test_size = 0.25, random_state = 42)

	rf = MultiOutputClassifier(RandomForestClassifier(n_estimators=1000,
	                                                  random_state=42))
	rf.fit(train_features, train_labels)
	predictions = rf.predict(test_features)
	errors = abs(predictions - test_labels)
	print("Mean Absolute Error: ", round(np.mean(errors), 2))
