# -*- coding: UTF-8 -*-
# Used packages and libraries
import bioread
import os.path
import shutil
import sys
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

sampling_rate=64

def preprocess_single(csv_data):
    filtered_signal = hp.filter_signal(csv_data, cutoff=[0.5, 1.5], sample_rate=sampling_rate, order=5, filtertype='bandpass', return_top=False)
    signals, e4_peaks = nk.ppg_process(filtered_signal, sampling_rate=sampling_rate)
    new_peaks_e4, clean_factor, threshold = filter_peaks(e4_peaks["PPG_Peaks"])
    hrv = nk.hrv(new_peaks_e4, sampling_rate=sampling_rate, show=False)
    print(hrv)

# convert RR's millisecond values to location of peaks
def rr_to_peaks(rri, sampl):
    res = [0]
    for n, rr in enumerate(rri):
        if n > 0:
            res.append((res[n - 1] + ((rr / 1000) * sampl)).astype(int))
    return res

# Function which is filterring the values of peaks and replace it with calculated ones
def filter_peaks(peaks):
    rri = np.diff(peaks) / sampling_rate * 1000
    # Variable to keep iteration of the filtering process
    thres = 1
    try:
        # Using threshold filter from HRV package. When strong (150ms) threshold is to low to interpolate,
        # it is increased
        filtered_rri = threshold_filter(rri, threshold='strong',
                                       local_median_size=3)
    except ValueError:
        try:
            # Using threshold filter from HRV package. When medium (250ms) threshold is to low to interpolate,
            # it is increased
            print("changing threshold to medium (250ms)")
            filtered_rri = threshold_filter(rri, threshold='medium',
                                            local_median_size=3)
            thres = 2
        except ValueError:
            try:
                # Using threshold filter from HRV package. When low (350ms) threshold is to low to interpolate,
                # it is increased
                print("changing threshold to low (350ms)")
                filtered_rri = threshold_filter(rri, threshold='low',
                                               local_median_size=3)
                thres = 3
            except ValueError:
                try:
                    # Using threshold filter from HRV package. When very low (450ms) threshold is to low to interpolate,
                    # the filter is changed to the moving median error
                    print("changing threshold to very low (450ms)")
                    filtered_rri = threshold_filter(rri, threshold='very low',
                                                   local_median_size=3)
                    thres = 4
                except ValueError:
                    try:
                        # Using moving median filter from HRV package, when the interpolation is failed,
                        # the whole pipeline process for this subject is failed and its rejected from further analysis
                        filtered_rri = moving_median(rri, order=3)
                        print("changing to moving avarege filter")
                        thres = 5
                    except ValueError:
                        raise ValueError("Subjects RRi's cannot be filtered break")
    # Calculating clean factor for each filtering process
    clean_factor = feature_filter_rri(rri, filtered_rri)
    # Ge new peaks from filtered RRs
    new_peaks = rr_to_peaks(filtered_rri, sampling_rate)
    new_peaks = new_peaks + peaks[:][0]
    # return new filtered peaks location, clean factor of filtering, iteration of filtering which passed
    return new_peaks, clean_factor, thres

# Calculating clean factor as a sum of differences between original and cleaned peaks
def feature_filter_rri(rri, filtered_rri):
    res = 0
    for rr, filtr_rr in zip(rri, filtered_rri):
        res += abs(rr - filtr_rr)
    return res / len(rri)

if __name__ == "__main__":
    path = sys.argv[1]
    data = pd.read_csv(path, sep=",")
    preprocess_single(data['x'])
    
