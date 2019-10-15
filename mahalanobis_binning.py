import numpy as np
import pandas as pd
import scipy.linalg as sp
from scipy.stats import chi2
import sys

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

critical_value = chi2.ppf((1-0.01), df=1)

contig_nameread = ['id', 'Actual taxon','option','len']
file1 = pd.read_csv('sample_data/simBG/simBG_view1.OFDEG', delimiter = ',') # file contain details about all the contigs
file2 = pd.read_csv('sample_data/simBG/simBG_view2.n4', delimiter = ',') # file contain details about all the contigs
file3 =  pd.read_csv('sample_data/simBG/sim.contig.ans', delimiter = '\t', names=contig_nameread, header=None) # file contain names

binnedFile = pd.read_csv("binned_points.csv")
df_2Tbinned = binnedFile[['id', 'ofdeg', 'gc', 'length','bin']]
df_2Tunbinned = file1[['id', 'ofdeg', 'gc', 'length']]
df_names = file3.iloc[: , [0, 1]]
# print(df_names);

df= df_2Tunbinned.head(50) # for debugging   
bins_array = df_2Tbinned['bin'].unique() # get the bins
i = 0
newdf = pd.DataFrame(columns=['id', 'ofdeg', 'gc', 'length','bin'])

while i < len(df.index): 
    row = df.loc[i,]
    lowest_variance = sys.float_info.max
    assigned_bin = None
    goodForBin = 0
    distance = None

    for bin in bins_array :
        label_bin = calculate_binwise_dist(bin)
        df1 = label_bin[['id', 'ofdeg', 'gc', 'length']].append(row, ignore_index=True) # add each unbinned contigs to bin and calculate distance
        df2 = df1[['ofdeg', 'gc']]
        df1['mahala'] = calculate_mahalanobis_dist(x=df2, data=df2[['ofdeg', 'gc']]) #dataframe with distance column
        dist = df1.loc[df1.index[-1], "mahala"]
        
        if dist < critical_value : # this will check whether contig has smaller distance or not
            goodForBin = 1
            x = dist # contig mahalanobis distances
            m = df1['mahala'].mean()  # mean
            n = len(df1['mahala'].index) #sample number 
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
        

newdf = newdf.merge(df_names, on="id", how = 'inner')
# print(newdf[['id', 'bin','mahala', 'Actual taxon']]) 

# check results for dataset
summery = binnedFile.drop_duplicates(subset='bin', keep="last")
summery_of_bins = summery[['bin','taxon']]
summery_of_bins.reset_index(inplace = True, drop = True)
# print(summery_of_bins)
results = newdf.merge(summery_of_bins, on="bin", how = 'inner')
# print(results[['id','bin','taxon','Actual taxon']])
results = results[['id','length','bin','mahala','taxon','Actual taxon']]
results['Is prediction correct']= (results['taxon']==results['Actual taxon'])
print(results[['id','length','bin','mahala','Is prediction correct']])

true_prediction = results['Is prediction correct'].values.sum() # true count
false_prediction = (~results['Is prediction correct']).values.sum() # false count
accuracy = 100 * true_prediction / (false_prediction + true_prediction)
print("Accuracy is: ", accuracy)