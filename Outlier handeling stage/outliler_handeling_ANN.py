# -*- coding: utf-8 -*-
"""
@author: PRABHATH
"""

import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import pandas as pd
from PIL import Image
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split

number_of_features = 4
data_location = "binned_data/"
data_set_name = "binned_points - Copy.csv"


def get_dataset_name():
    return data_location + data_set_name


def get_color(index):
    switcher = {
        0: 'r',
        1: 'g',
        2: 'b',
        3: 'k',
        4: 'c',
        5: 'crimson',
        6: 'dodgerblue',
        7: 'deeppink',
        8: 'fuchsia',
        9: 'olive',
        10: 'coral',
        11: 'orange',
        12: 'navy',
        13: 'darkviolet',
        14: 'maroon',
        15: 'darkorange',
    }
    return switcher.get(index, 'g')


def get_color_for_bins(bin_set, bins):
    colors = []

    for bin in bins.tolist():
        colors.append(get_color(bin_set.index(bin)))

    return colors


# Visualize data
def visualize():
    data_frame = pd.read_csv(get_dataset_name())
    n = data_frame.shape[0]

    bin_values = data_frame['bin'].values
    bin_list = bin_values.tolist()
    bin_set = list(set(bin_list))
    colors_ = get_color_for_bins(bin_set, data_frame['bin'])

    plt.subplot(2, 2, 1)
    plt.scatter(data_frame['gc'], data_frame['ofdeg'], s = 5, c=colors_, alpha = 0.5, edgecolors='none')

    plt.subplot(2, 2, 2)
    plt.scatter(data_frame[['r']], data_frame[['q']], s = 5, c=colors_, alpha = 0.5, edgecolors='none')

    plt.subplot(2, 2, 3)
    plt.scatter(data_frame[['r']], data_frame[['gc']], s = 10, c=colors_, alpha = 0.5, edgecolors='none')

    plt.subplot(2, 2, 4)
    plt.scatter(data_frame[['gc']], data_frame[['q']], s = 10, c=colors_, alpha = 0.5, edgecolors='none')
    plt.show()


visualize()


# Reading the dataset
def read_dataset():
    raw_df = pd.read_csv(get_dataset_name())
    df = raw_df[['bin', 'ofdeg', 'gc', 'r', 'q', 'NA', 'NA.1', 'ns', 'rms']]
    # df = raw_df[['bin', 'ofdeg', 'gc', 'r', 'q']]

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
train_x, test_x, train_y, test_y = train_test_split(X, Y, test_size=0.01, random_state=415)

# Inspect the shape of the testing and training
print('train_x shape', train_x.shape)
print('train_y shape', train_y.shape)
print('test_x shape', test_x.shape)

# Define the important parameters and variable to work with the tensors
learning_rate = 0.3
training_epochs = 300
cost_history = np.empty(shape=[1], dtype=float)
n_dim = X.shape[1]
print('n_dim', n_dim)
n_class = Y.shape[1]
print('n_class', n_class)
model_path = 'Model/ANN_model'

# Define the hidden layers

n_hidden_1 = 5000
n_hidden_2 = 100

x = tf.placeholder(tf.float32, [None, n_dim])
W = tf.Variable(tf.zeros([n_dim, n_class]))
b = tf.Variable(tf.zeros([n_class]))
y_ = tf.placeholder(tf.float32, [None, n_class])


# Define the model
def multilayer_perceptron(x, weights, biases):

    layer_1 = tf.add(tf.matmul(x, weights['h1']), biases['b1'])
    layer_1 = tf.nn.sigmoid(layer_1)

    layer_2 = tf.add(tf.matmul(layer_1, weights['h2']), biases['b2'])
    layer_2 = tf.nn.sigmoid(layer_2)

    out_layer = tf.add(tf.matmul(layer_2, weights['out']), biases['out'])
    return out_layer


# Define the weights and the biases for each layer

weights = {
    'h1': tf.Variable(tf.truncated_normal([n_dim, n_hidden_1])),
    'h2': tf.Variable(tf.truncated_normal([n_hidden_1, n_hidden_2])),
    'out': tf.Variable(tf.truncated_normal([n_hidden_2, n_class]))
}
biases = {
    'b1': tf.Variable(tf.truncated_normal([n_hidden_1])),
    'b2': tf.Variable(tf.truncated_normal([n_hidden_2])),
    'out': tf.Variable(tf.truncated_normal([n_class]))
}


def tensor_board_visualization(tf_sess):
    File_Writer = tf.summary.FileWriter('graph', tf_sess.graph)


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
print("Test size: ", np.sum(test_y))
print(np.sum(test_y, 0))
print("Train size: ", np.sum(train_y))
print(np.sum(train_y, 0))
tensor_board_visualization(sess)


for epoch in range(training_epochs):
    sess.run(training_step, feed_dict={x: train_x, y_: train_y})
    cost = sess.run(cost_function, feed_dict={x: train_x, y_: train_y})
    cost_history = np.append(cost_history, cost)
    print(y, y_)
    correct_prediction = tf.equal(tf.argmax(y, 1), tf.argmax(y_, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

    pred_y = sess.run(y, feed_dict={x: test_x})
    mse = tf.reduce_mean(tf.square(pred_y - test_y))
    mse_ = sess.run(mse)
    mse_history.append(mse_)
    accuracy = (sess.run(accuracy, feed_dict={x: test_x, y_: test_y}))
    accuracy_history.append(accuracy)

    print('epoch : ', epoch, '-', 'cost: ', cost, ' - MSE:', mse_, '- Train Accuracy: ', accuracy * 100)

save_path = saver.save(sess, model_path)
print("Model saved in file: %s" % save_path)

# print(sess.run(prediction, feed_dict={x:test_x}))
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
