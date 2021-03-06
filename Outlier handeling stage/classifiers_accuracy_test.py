# -*- coding: utf-8 -*-
"""
@author: PRABHATH
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split


data_location = "binned_data/"
data_set_name = "binned_points - Copy.csv"
number_of_features = 2


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

# Splitting the dataset into the Training set and Test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=0)

# Feature Scaling
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)


classifiers = []
confusion_matrices = []
accuracy_scores = []

# Fitting classifier to the Training set
classifier = GaussianNB()
classifier.fit(X_train, y_train)
classifiers.append(classifier)

classifier = DecisionTreeClassifier(criterion='entropy', random_state=0)
classifier.fit(X_train, y_train)
classifiers.append(classifier)

classifier = LogisticRegression(random_state=0)
classifier.fit(X_train, y_train)
classifiers.append(classifier)

classifier = KNeighborsClassifier(n_neighbors=5, metric='minkowski', p=2)
classifier.fit(X_train, y_train)
classifiers.append(classifier)

classifier = SVC(kernel='linear', random_state=0)
classifier.fit(X_train, y_train)
classifiers.append(classifier)

classifier = RandomForestClassifier(n_estimators=10, criterion='entropy', random_state=0)
classifier.fit(X_train, y_train)
classifiers.append(classifier)

accuracies = []


for classifier in classifiers:
    # Predicting the Test set results
    y_pred = classifier.predict(X_test)

    # Making the Confusion Matrix
    confusion_matrix_ = confusion_matrix(y_test, y_pred)
    confusion_matrices.append(confusion_matrix_)
    accuracy = accuracy_score(y_test, y_pred)
    print("Test Accuracy", accuracy)
    accuracies.append([classifier, accuracy])

    y_pred = classifier.predict(X_train)

    accuracy = accuracy_score(y_train, y_pred)
    print("Train Accuracy", accuracy)

for n in range(1, 11):
    classifier = KNeighborsClassifier(n_neighbors=n, metric='minkowski', p=2)
    classifier.fit(X_train, y_train)
    classifiers.append(classifier)

    y_pred = classifier.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(n, accuracy)
