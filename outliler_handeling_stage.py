# -*- coding: utf-8 -*-
"""
Created on Sat Jan 05 11:44:20 2019
@author: PRABHA
"""

import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split

number_of_features = 3


# Reading the dataset
def read_dataset():
    raw_df = pd.read_csv("binned_points.csv")
    df = raw_df[['bin', 'ofdeg', 'gc', 'length']]

    X = df[df.columns[1:number_of_features+1]].values
    y = df[df.columns[0]]

    # Encode the dependent variable
    encoder = LabelEncoder()
    encoder.fit(y)
    y = encoder.transform(y)
    Y = one_hot_encode(y)

    return X, Y


# Define the encoder function.
def one_hot_encode(labels):
    n_labels = len(labels)
    n_unique_labels = len(np.unique(labels))
    one_hot_encod = np.zeros((n_labels, n_unique_labels))
    one_hot_encod[np.arange(n_labels), labels] = 1
    return one_hot_encod


# Read the dataset
X, Y = read_dataset()

# Shuffle the dataset
X, Y = shuffle(X, Y, random_state=1)

# Convert the dataset train and test part
train_x, test_x, train_y, test_y = train_test_split(X, Y, test_size=0.20, random_state=415)

# Inspect the shape of the testing and training
print('train_x shape', train_x.shape)
print('train_y shape', train_y.shape)
print('test_x shape', test_x.shape)

# Define the important parameters and variable to work with the tensors
learning_rate = 0.3
training_epochs = 700
cost_history = np.empty(shape=[1], dtype=float)
n_dim = X.shape[1]
print('n_dim', n_dim)
n_class = 7
model_path = 'Model/ANN_model'

# Define the hidden layers

n_hidden_1 = 100

x = tf.placeholder(tf.float32, [None, n_dim])
W = tf.Variable(tf.zeros([n_dim, n_class]))
b = tf.Variable(tf.zeros([n_class]))
y_ = tf.placeholder(tf.float32, [None, n_class])


# Define the model
def multilayer_perceptron(x, weights, biases):

    layer_1 = tf.add(tf.matmul(x, weights['h1']), biases['b1'])
    layer_1 = tf.nn.sigmoid(layer_1)

    out_layer = tf.matmul(layer_1, weights['out']) + biases['out']
    return out_layer


# Define the weights and the biases for each layer

weights = {
    'h1': tf.Variable(tf.truncated_normal([n_dim, n_hidden_1])),
    'out': tf.Variable(tf.truncated_normal([n_hidden_1, n_class]))
}
biases = {
    'b1': tf.Variable(tf.truncated_normal([n_hidden_1])),
    'out': tf.Variable(tf.truncated_normal([n_class]))
}

# Initialize all the variables
init = tf.global_variables_initializer()
saver = tf.train.Saver()

# Call modle
y = multilayer_perceptron(x, weights, biases)

# Define the cost function and optimizer
cost_function = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=y, labels=y_))
training_step = tf.train.GradientDescentOptimizer(learning_rate).minimize(cost_function)

sess = tf.Session()
sess.run(init)

# Calculate the cost and thr accuracy for each epoch
mse_history = []
accuracy_history = []

for epoch in range(training_epochs):
    sess.run(training_step, feed_dict={x: train_x, y_: train_y})
    cost = sess.run(cost_function, feed_dict={x: train_x, y_: train_y})
    cost_history = np.append(cost_history, cost)
    correct_prediction = tf.equal(tf.argmax(y, 1), tf.argmax(y_, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

    pred_y = sess.run(y, feed_dict={x: test_x})
    mse = tf.reduce_mean(tf.square(pred_y - test_y))
    mse_ = sess.run(mse)
    mse_history.append(mse_)
    accuracy = (sess.run(accuracy, feed_dict={x: train_x, y_: train_y}))
    accuracy_history.append(accuracy)

    print('epoch : ', epoch, '-', 'cost: ', cost, ' - MSE:', mse_, '- Train Accuracy: ', accuracy)

save_path = saver.save(sess, model_path)
print("Model saved in file: %s" % save_path)


plt.plot(mse_history, 'r')
plt.show()
plt.plot(accuracy_history)
plt.show()

# Print the final accuracy
correct_prediction = tf.equal(tf.argmax(y, 1), tf.argmax(y_, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
print("Test Accuracy: ", (sess.run(accuracy, feed_dict={x: test_x, y_:test_y})))

# Print the final mean square error
pred_y = sess.run(y, feed_dict={x: test_x})
mse = tf.reduce_mean(tf.square(pred_y - test_y))
print("MSE: %.4f" % sess.run(mse))
