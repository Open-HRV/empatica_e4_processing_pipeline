import os.path
import shutil
from os import listdir

import heartpy as hp
import neurokit2 as nk
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


import mat4py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
from matplotlib import pyplot
from scipy import signal
import statistics as stat

def filter_high_values(df):
    result = pd.DataFrame([], columns=["rest_e4_ppg_HRV_RMSSD", "rest_e4_ppg_HRV_SDNN", "rest_e4_ppg_HRV_LFn", "rest_e4_ppg_HRV_HFn", "rest_e4_ppg_HRV_LnHF", "rest_e4_ppg_HRV_LFHF"])
    for row, value in df.iterrows():
        good = True
        for i, val in enumerate(value):
            t1 = means[i] + 1.5 * stds[i]
            t2 = means[i] - 1.5 * stds[i]
            if val > t1 or val < t2:
                print("wrong")
                good = False
                break
        if good:
            df2 = pd.DataFrame(value).T.reset_index(drop=True)
            result = pd.concat([df2, result], ignore_index=True)
            good = True
    return result
    #for value in row


data = pd.read_csv(r"biopac_vs_e4_full.csv", sep=",")
data4 = data[["rest_e4_ppg_HRV_RMSSD", "rest_e4_ppg_HRV_SDNN", "rest_e4_ppg_HRV_LFn", "rest_e4_ppg_HRV_HFn", "rest_e4_ppg_HRV_LnHF", "rest_e4_ppg_HRV_LFHF"]]

means = data4.mean(axis=0)
stds = data4.std(axis=0)

data5 = filter_high_values(data4)
print(data5)
hist = data5.hist(bins = 20)
plt.show()

data5 = StandardScaler().fit_transform(data5)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(data5)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['pc1', 'pc2'])
print(principalDf)

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)

ax.scatter(principalDf['pc1'], principalDf['pc2'])
plt.show()
#print(data5)


