# -*- coding: utf-8 -*-
"""
@author: PRABHATH
"""

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler


data_location = "Outlier handeling stage/binned_data/"
data_set_name = "binned_points - Copy.csv"
number_of_features = 2
length_consider = 1000


def get_dataset_name():
    return data_location + data_set_name


# Define the encoder function.
def one_hot_encode(labels):
    n_labels = len(labels)
    n_unique_labels = len(np.unique(labels))
    one_hot_encod = np.zeros((n_labels, n_unique_labels))
    one_hot_encod[np.arange(n_labels), labels] = 1
    return one_hot_encod


# Importing the dataset
dataset = pd.read_csv(get_dataset_name())
df = dataset[['bin', 'ofdeg', 'gc', 'r', 'q']]

X = df[df.columns[1:number_of_features + 1]].values
y = df[df.columns[0]]

# Feature Scaling
sc = StandardScaler()
X_train = X
y_train = y

# Fitting classifier to the Training set
classifier = KNeighborsClassifier(n_neighbors=10, metric='minkowski', p=2)
classifier.fit(X_train, y_train)

# file contain details about all the contigs
unbinned_contigs_dataframe = pd.read_csv('sample_data/simBG/simBG_output.unbinned_contiges', sep='\t')
unbinned_contig_count = len(unbinned_contigs_dataframe.index)
drop_index_names_1 = unbinned_contigs_dataframe[unbinned_contigs_dataframe['length'] < length_consider].index
unbinned_contigs_dataframe.drop(drop_index_names_1, inplace=True)
unbinned_contigs_dataframe.reset_index(inplace=True, drop=True)
predict_data_set = unbinned_contigs_dataframe[['ofdeg', 'gc', 'r', 'q']]

contig_bins_dataframe = pd.read_csv('sample_data/simBG/simBG_output.L2_BINS', sep='\t', header=None)

contig_ans = pd.read_csv('sample_data/simBG/sim.contig.ans', sep='\t', header=None)

merged_unbinned_dataframe_with_ans = pd.merge(unbinned_contigs_dataframe, contig_ans, left_on='id',
                                              right_on=0, how='left').drop(0, axis=1)
# merged_unbinned_dataframe_with_ans = merged_unbinned_dataframe_with_ans.dropna()
predict_X = df[df.columns[0:number_of_features]].values

binned_contigs_with_ = pd.merge(dataset, contig_ans, left_on='id', right_on=0, how='left').drop(0, axis=1)
binned_contigs_combinations = binned_contigs_with_[['bin', 1]]
# binned_contigs_combinations = binned_contigs_combinations.groupby(['bin', 1])
binned_contigs_combinations = binned_contigs_combinations.drop_duplicates(subset='bin', keep="last")

merged_unbinned_dataframe_with_ans = pd.merge(merged_unbinned_dataframe_with_ans, binned_contigs_combinations,
                                              on=1, how='left')
merged_unbinned_dataframe_with_ans = merged_unbinned_dataframe_with_ans.dropna()
print(merged_unbinned_dataframe_with_ans)

unbinned_df = merged_unbinned_dataframe_with_ans[['bin', 'ofdeg', 'gc', 'r', 'q']]
X_unbinned = unbinned_df[unbinned_df.columns[1:number_of_features + 1]].values
y_unbinned = unbinned_df[unbinned_df.columns[0]]

y_pred = classifier.predict(X_unbinned)
confusion_matrix_ = confusion_matrix(y_unbinned, y_pred)
accuracy = accuracy_score(y_unbinned, y_pred)


