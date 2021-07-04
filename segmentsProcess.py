import bioread
import os.path
import shutil
from os import listdir

import heartpy as hp
import neurokit2 as nk
from hrv.filters import moving_median
from hrv.filters import quotient
from hrv.filters import threshold_filter
from hrv.filters import moving_average

import mat4py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
from matplotlib import pyplot
from scipy import signal
import statistics as stat
from functions import *


path = "C:\\Users\\Przemek\\Documents\\magisterka\\Segments dla Przemka\\"

segments = [name for name in listdir(path)]

for segment in segments:
    ppg_data = pd.read_csv(path + segment, sep=",")
    ppg = ppg_data['x']
   # plt.plot(ppg)
    filtered_e4 = hp.filter_signal(ppg, cutoff=[0.7,2],
                                     sample_rate=sampling_rate_use_acc, order=5, filtertype='bandpass', return_top=False)
    #filtered_e4_2 = hp.filter_signal(filtered_e4, cutoff=2, sample_rate=64, filtertype='lowpass', order=3)

    ppg_cleaned = nk.ppg_clean(filtered_e4, sampling_rate=64)
    #ppg_cleaned = hp.enhance_peaks(ppg_cleaned, iterations=1)

    signals, e4_peaks = nk.ppg_process(ppg_cleaned, sampling_rate=sampling_rate_use_acc)
    plt.plot(signals)
    plt.show()
    new_peaks, clean_factor = filter_peaks(e4_peaks)
    print(clean_factor)
    try:
        hrv = nk.hrv(new_peaks, sampling_rate=sampling_rate_use_acc, show=False)#show=True)
    except ValueError:
        pass

    hrv["clean_factor"] = clean_factor
    print(hrv)
   # saveToFile(segment, "e4", hrv)
# plt.plot(ppg_cleaned)
    #plt.show()