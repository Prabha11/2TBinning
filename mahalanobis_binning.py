import sys
import numpy as np
import pandas as pd
import scipy.linalg as sp
import seaborn as sns
from pylab import *
import matplotlib.pyplot as plt
from scipy.stats import chi2
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.gridspec import GridSpec
from scipy.interpolate import UnivariateSpline

def calculate_binwise_dist (label) :
    arr = df_2Tbinned[df_2Tbinned['bin'] == label]
    return arr

def calculate_mahalanobis_dist (x=None, data=None, cov=None) :
    """Compute the Mahalanobis Distance between each row of x and the data  
    x    : vector or matrix of data with, say, p columns.
    data : ndarray of the distribution from which Mahalanobis distance of each observation of x is to be computed.
    cov  : covariance matrix (p x p) of the distribution. If None, will be computed from data.
    D^2 = (x−m)^T * C^(−1)* (x−m)
    """
    x_minus_mu = x - np.mean(data)
    if not cov:
        cov = np.cov(data.values.T)
    inv_covmat = sp.inv(cov)
    left_term = np.dot(x_minus_mu, inv_covmat)
    mahal = np.dot(left_term, x_minus_mu.T)
    return mahal.diagonal()

critical_value = chi2.ppf((1-0.01), df=2)

# UNBINNED DATA FEATURE EXTRACTION
# critical_value = 4
length_consider = 1000

contig_nameread = ['id', 'Actual taxon','option','len']
file1 = pd.read_csv('sample_data/simBG/simBG_unbinned_contigs.OFDEG', delimiter = ',') # file contain details about all the contigs
unbinned_contig_count = len(file1.index)
dropindexNames1 = file1[ file1['length'] < length_consider ].index
file1.drop(dropindexNames1 , inplace=True)
file1.reset_index(inplace = True, drop = True)
# print(file1)

file2 = pd.read_csv('sample_data/simBG/simBG_unbinned_contigs.n4', delimiter = ',') # file contain details about all the contigs
dropindexNames2 = file2[ file2['length'] < length_consider ].index
file2.drop(dropindexNames2 , inplace=True)
file2.reset_index(inplace = True, drop = True)
# print(file2)
df_tnf = file2.drop(['length', 'ns'], axis=1)
# print(df_tnf)

df_2Tunbinned = file1[['id', 'ofdeg', 'gc', 'length']]

# tnf pca for unbinned data
ids = df_tnf.loc[:,['id']] # Separating out the ids
tnt_frquencies = df_tnf.drop(['id'], axis=1).values # Separating out the tnf freq
tranformed_tnf = StandardScaler().fit_transform(tnt_frquencies) # Standardizing the features
pca = PCA(n_components=3)
principalComponents = pca.fit_transform(tranformed_tnf)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2','principal component 3'])
principalDf['id'] = ids
df_2Tunbinned = df_2Tunbinned.merge(principalDf, on="id", how = 'inner') # add principle component cols
# df_2Tunbinned = file1[['id', 'ofdeg', 'gc', 'length']];
# print(df_2Tunbinned.head(10))

# BINNED DATA FEATURE EXTRACTION

binnedFile = pd.read_csv("binned_points.csv") # binned file
df_2Tbinned = binnedFile[['id', 'ofdeg', 'gc', 'length','bin']]
file4 = pd.read_csv('sample_data/simBG/simBG_view2.n4', delimiter = ',') # file contain details about all the contigs
df_tnf_binned = file4.drop(['length', 'ns'], axis=1)

# tnf pca for binned data
ids = df_tnf_binned.loc[:,['id']] # Separating out the ids
tnt_frquencies = df_tnf_binned.drop(['id'], axis=1).values # Separating out the tnf freq
tranformed_tnf = StandardScaler().fit_transform(tnt_frquencies) # Standardizing the features
pca = PCA(n_components=3)
principalComponents = pca.fit_transform(tranformed_tnf)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2','principal component 3'])
principalDf['id'] = ids
df_2Tbinned = df_2Tbinned.merge(principalDf, on="id", how = 'inner') # add principle component cols
# print(df_2Tbinned.head(10))

file3 =  pd.read_csv('sample_data/simBG/sim.contig.ans', delimiter = '\t', names=contig_nameread, header=None) # file contain names
df_names = file3.iloc[: , [0, 1]]

binned_contig_count = len(binnedFile.index)

# COMPUTATION

df= df_2Tunbinned # for debugging   
bins_array = df_2Tbinned['bin'].unique() # get the bins
i = 0
newdf = pd.DataFrame(columns=['id', 'ofdeg', 'gc', 'length','bin'])
mahala_dist = []

while i < len(df.index): 
    row = df.loc[i,]
    lowest_variance = sys.float_info.max
    assigned_bin = None
    goodForBin = 0
    distance = None
    minD = sys.float_info.max; # this variable for store each unbinned contigs mahalanobis dist
    
    for bin in bins_array :
        label_bin = calculate_binwise_dist(bin)
        df1 = label_bin[['id', 'ofdeg', 'gc', 'length','principal component 1', 'principal component 2','principal component 3']].append(row, ignore_index=True) # add each unbinned contigs to bin and calculate distance
        df2 = df1[['ofdeg', 'gc','principal component 1', 'principal component 2','principal component 3']]
        df1['mahala'] = calculate_mahalanobis_dist(x=df2, data=df2[['ofdeg', 'gc','principal component 1', 'principal component 2','principal component 3']]) #dataframe with distance column
        dist = df1.loc[df1.index[-1], "mahala"]
        d = dist
        
        if d < minD : # get the minimum mahalanobis dist
            minD = d
        
        if dist < critical_value : # this will check whether contig has smaller distance or not
            goodForBin = 1
            x = dist; # contig mahalanobis distances
            m = df1['mahala'].mean();  # mean
            n = len(df1['mahala'].index); #sample number 
            variance = (x-m)*(x-m) / n 
            
            if variance < lowest_variance :
                lowest_variance = variance
                assigned_bin = bin
                distance = dist
    i += 1
    
    if goodForBin == 1 : # contigs bin in our model
        row['bin'] = assigned_bin
        row['mahala'] = distance
        newdf = newdf.append(row, ignore_index=True)
    
    mahala_dist.append(minD)

newdf = newdf.merge(df_names, on="id", how = 'inner')
# print(newdf[['id', 'bin','mahala', 'Actual taxon']]); 
pd.options.mode.chained_assignment = None  # default='warn'
df['mahala'] = mahala_dist # assign distance for to column 'mahala' in df dataframe
# print(df)
# print(newdf)


# check results for dataset

summery = binnedFile.drop_duplicates(subset='bin', keep="last")
summery_of_bins = summery[['bin','taxon']]
summery_of_bins.reset_index(inplace = True, drop = True)
# print(summery_of_bins)
results = newdf.merge(summery_of_bins, on="bin", how = 'inner')
# print(results[['id','bin','taxon','taxon']])
results = results[['id','length','bin','mahala','taxon','Actual taxon']]
results['Is prediction correct']= (results['taxon']==results['Actual taxon'])
# print(results[['id','bin','taxon','taxon','Is prediction correct']])
 
true_prediction = results['Is prediction correct'].values.sum() # true count
false_prediction = (~results['Is prediction correct']).values.sum() # false count
binned_count_our_model = false_prediction + true_prediction
accuracy = 100 * true_prediction / binned_count_our_model
print("Accuracy is:% 5.2f" % accuracy, "%")


# A N A L Y S E

plt.figure(1, figsize=(13,10))
the_grid = GridSpec(2, 2)

# chart 1
labels = ['Unbinned contigs', 'Binned contigs']
sizes = [unbinned_contig_count, binned_contig_count]
colors = ['mediumseagreen','yellowgreen']
plt.subplot(the_grid[0, 0], aspect=1)
plt.pie(sizes, colors=colors,autopct='%1.1f%%', labels= sizes, counterclock=False)
plt.title('2T Binning method results', fontsize=16)
plt.legend(labels,loc="best")

# # chart 2
# labels = ['2T binning method', 'Our method']
# sizes = [binned_contig_count, binned_count_our_model]
# colors = ['burlywood', 'peru']
# plt.subplot(the_grid[0, 1], aspect=1)
# plt.pie(sizes, colors=colors,autopct='%1.1f%%', labels= sizes, counterclock=False)
# plt.title('2T Binning vs Our model', fontsize=16)
# plt.legend(labels,loc="best")

# chart 3
labels = ['Unbinned contigs', 'Proposed method']
sizes = [unbinned_contig_count - binned_count_our_model, binned_contig_count+ binned_count_our_model]
colors = ['steelblue', 'lightskyblue']
plt.subplot(the_grid[0, 1], aspect=1)
plt.pie(sizes, colors=colors,autopct='%1.1f%%', labels= sizes, counterclock=False)
plt.title('Proposed method results (critical value = 9 & length >1000)', fontsize=16)
plt.legend(labels,loc="best")


plt.suptitle('Analyze number of contigs binned', fontsize=20)
plt.show()

#  H I S T O G R A M

# fig, axes = plt.subplots(1, 2, figsize=(10,2.5), dpi=100, sharex=True, sharey=True)
x1 = df[ 'mahala']
plt.rcParams["figure.figsize"] = [15,7]
plt.suptitle('Probability Density function and Histogram of Mahalanobis distance', fontsize=20)

subplot(1,3,1)
ax = x1.plot.kde()
title('All the unbinneed contigs')
xlabel('Distance')

x2 = newdf['mahala']
subplot(1,3,2)
# plt.subplot(the_grid[0, 1], aspect=1)
ax = x2.plot.kde()
title('Binned contigs in our model')
xlabel('Distance')

subplot(1,3,3)
plt.hist(x1, color = 'blue', edgecolor = 'black', bins = 50)
# seaborn histogram
sns.distplot(x1, hist=True, kde=False, bins=50, color = 'blue', hist_kws={'edgecolor':'black'})
plt.title('Histogram of Distances')
plt.xlabel('Distances')
plt.ylabel('Frequency')


# # A N A L Y S E   W I T H   C H A N G I N G   L E N G T H
# x = length_consider
# y= accuracy
# z = binned_contig_count+ binned_count_our_model;
# data = [[x, y, z]] 
  
# # Create the pandas DataFrame 
# df = pd.DataFrame(data, columns = ['length', 'Accuracy', 'Number of contigs binned'])  # first time with headers afetr that without
# df.to_csv('sample_data/simBG/change_length.csv', header=False, mode = 'a')


# A N A L Y S E   W I T H   T H R E S H O L D   M A H A L A N O B I S    D I S T A N C E
# x = critical_value
# y= accuracy
# z = binned_contig_count+ binned_count_our_model;
# data = [[x, y, z]] 
  
# # Create the pandas DataFrame 
# df = pd.DataFrame(data, columns = ['critical_value', 'Accuracy', 'Number of contigs binned'])  # first time with headers afetr that without
# df.to_csv('sample_data/simBG/change_critical_value.csv', header=False, mode = 'a')


# # A N A L Y S E   W I T H   C H A N G I N G   L E N G T H

df_len = pd.read_csv('sample_data/simBG/change_length.csv')
length = df_len.iloc[:, 1].to_numpy()
accuracy = df_len.iloc[:, 2].to_numpy()
lenN = df_len.iloc[:, 3]
lenN = (100 * lenN) / 40000
number = lenN.to_numpy()
# print(len)


x_min = df_len.iloc[:, 1].min() -100
x_max = df_len.iloc[:, 1].max() +100
y1_min =  df_len.iloc[:, 2].min() - 5
y1_max = df_len.iloc[:, 2].max() + 5
y2_min =  lenN - 2
y2_max = lenN.max() + 2


# # scale 'number' with 'accuracy'
# for index, val in np.ndenumerate(number):
#     number[index]= y1_min + ((val - y2_min) / (y2_max - y2_min)) * (y1_max - y1_min)

    
# calculate polynomial 1
z1 = np.polyfit(length, accuracy, 2)
f1 = np.poly1d(z1) 

# calculate polynomial 2
z2 = np.polyfit(length, number, 2)
f2 = np.poly1d(z2) 


# calculate new x's and y's
x_new = np.linspace(x_min, x_max, 100)
y1_new = f1(x_new)
y2_new = f2(x_new)

# plt.title('Accuracy & Number of contigs binned with Different lengths',fontsize=18)
# plt.xlabel('Lengths',fontsize=16)
# plt.ylabel('Percentages',fontsize=16)
# plt.plot(length, accuracy, 'o',markerfacecolor='red',markersize=10)
# plt.plot( x_new, y1_new, '-',color='red',linewidth=3)
# plt.plot(length, number, 'o',markerfacecolor='blue',markersize=10)
# plt.plot( x_new, y2_new, '-',color='blue',linewidth=3)
# plt.legend(['data points of Accuarcy', 'Interpolated curve of Accuracy', '% of Binned contig count', 'Interpolated curve of % Binned contig coun'], loc = 'best')
# plt.show()

fig, ax1 = plt.subplots()
plt.title('Accuracy & Number of contigs binned with Different lengths (Critical value = 9)',fontsize=18)
color = 'tab:red'
ax1.set_xlabel('Length of contigs',fontsize='large', fontweight='bold')
ax1.set_ylabel('Accuracy (%)', color=color,fontsize='large', fontweight='bold')
plt.plot(length, accuracy, 'o',markerfacecolor='red',markersize=10)
plt.plot( x_new, y1_new, '-',color='red',linewidth=3)
ax1.tick_params(axis='y', labelcolor=color)

ax1.legend(['Accuarcy', 'Interpolated curve of Accuracy'], loc = 1)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Contigs binned %', color=color,fontsize='large', fontweight='bold')  # we already handled the x-label with ax1
plt.plot(length, number, 'o',markerfacecolor='blue',markersize=10)
plt.plot( x_new, y2_new, '-',color='blue',linewidth=3)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.legend(['Contigs %', 'Interpolated curve of contigs binned %',], loc = 2)

plt.show()

# A N A L Y S E   W I T H   T H R E S H O L D   M A H A L A N O B I S    D I S T A N C E

df_threshold = pd.read_csv('sample_data/simBG/change_critical_value.csv')
critical_value = df_threshold.iloc[:, 1].to_numpy()
accuracy = df_threshold.iloc[:, 2].to_numpy()
lenN = df_threshold.iloc[:, 3]
lenN = (100 * lenN) / 40000
number = lenN.to_numpy()

x_min = df_threshold.iloc[:, 1].min()-1
x_max = df_threshold.iloc[:, 1].max() +1
y1_min =  df_threshold.iloc[:, 2].min() - 5
y1_max = df_threshold.iloc[:, 2].max() + 5
y2_min =  lenN - 2
y2_max = lenN.max() + 2


# # scale 'number' with 'accuracy'
# for index, val in np.ndenumerate(number):
#     number[index]= y1_min + ((val - y2_min) / (y2_max - y2_min)) * (y1_max - y1_min)

    
# calculate polynomial 1
z1 = np.polyfit(critical_value, accuracy, 3)
f1 = np.poly1d(z1) 

# calculate polynomial 2
z2 = np.polyfit(critical_value, number, 3)
f2 = np.poly1d(z2) 


# calculate new x's and y's
x_new = np.linspace(x_min, x_max, 100)
y1_new = f1(x_new)
y2_new = f2(x_new)


# plt.title('Accuracy & Number of contigs binned with Different critical values',fontsize=18)
# plt.xlabel('Critical Values',fontsize=16)
# plt.ylabel('Percentages',fontsize=16)
# plt.plot(critical_value, accuracy, 'o',markerfacecolor='yellow',markersize=10)
# plt.plot( x_new, y1_new, '-',color='green',linewidth=3)
# plt.plot(critical_value, number, 'o',markerfacecolor='blue',markersize=10)
# plt.plot( x_new, y2_new, '-',color='blue',linewidth=3)
# plt.legend(['data points of Accuarcy', 'Interpolated curve of Accuracy', '% of Binned contig count', 'Interpolated curve of % Binned contig coun'], loc = 'best')
# plt.show()

fig, ax1 = plt.subplots()
plt.title('Accuracy & Number of contigs binned with Different critical values (Length > 1000)',fontsize=18)
color = 'tab:green'
ax1.set_xlabel('Critical Values',fontsize='large', fontweight='bold')
ax1.set_ylabel('Accuracy (%)', color=color,fontsize='large', fontweight='bold')
plt.plot(critical_value, accuracy, 'o',markerfacecolor=color,markersize=10)
plt.plot( x_new, y1_new, '-',color=color,linewidth=3)
ax1.tick_params(axis='y', labelcolor=color)

ax1.legend(['Accuarcy ', 'Interpolated curve of Accuracy'], loc = 2)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Contigs binned %', color=color,fontsize='large', fontweight='bold')  # we already handled the x-label with ax1
plt.plot(critical_value, number, 'o',markerfacecolor=color,markersize=10)
plt.plot( x_new, y2_new, '-',color=color,linewidth=3)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.legend(['Contigs %', 'Interpolated curve of contigs binned %',], loc = 1)

plt.show()
