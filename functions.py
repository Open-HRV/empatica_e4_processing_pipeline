#functions.py

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

run_rest = True; #tbd
run_ert = False; #tbd
plot_hrv = False;
sampling_rate_use_acc=64

def downsample(base_signal, signal_under, f_rate):
    indexes = (signal_under.index * f_rate).astype(int)
    indexes = indexes[indexes < base_signal.shape[0]]
    return base_signal.iloc[indexes].reset_index(drop=True)


def printHead(signal):
    plt.plot(signal.head(1000))
    plt.show()

def rrToPeaks(rri, sampl):
    res = [0]
    for n, rr in enumerate(rri):
        if n > 0:
            res.append((res[n - 1] + ((rr / 1000) * sampl)).astype(int))
    return res

def filter_peaks(peaks):
    rri = np.diff(peaks / sampling_rate_use_acc * 1000)

    try:
        filtered_rri = threshold_filter(rri, threshold='strong',
                                       local_median_size=10)
    except ValueError:
        try:
            print("changing threshold to medium (250ms)")
            filtered_rri = threshold_filter(rri, threshold='medium',
                                            local_median_size=10)
        except ValueError:
            try:
                print("changing threshold to low (350ms)")
                filtered_rri = threshold_filter(rri, threshold='low',
                                               local_median_size=10)
            except ValueError:
                try:
                    print("changing threshold to very low (450ms)")
                    filtered_rri = threshold_filter(rri, threshold='very low',
                                                   local_median_size=10)
                except ValueError:
                    try:
                        filtered_rri = moving_median(rri, order=3)
                        print("changing to moving avarege filter")
                    except ValueError:
                        raise ValueError("Subjects RRi's cannot be filtered break")
    clean_factor = feature_filter_rri(rri, filtered_rri)
    print(clean_factor)
    new_peaks = rrToPeaks(filtered_rri, sampling_rate_use_acc)
    new_peaks = new_peaks + peaks[:][0]

    return new_peaks, clean_factor


def feature_filter_rri(rri, filtered_rri):
    res = 0;
    for rr, filtr_rr in zip(rri, filtered_rri):
        res += abs(rr - filtr_rr)
    return res / len(rri)

def filter_signals(bp_ecg, bp_ppg, e4_ppg, use_acc=False):
    ecg_down = downsample(bp_ecg, e4_ppg, 2000 / 64)
    ppg_down = downsample(bp_ppg, e4_ppg, 2000 / 64)

    sampling_rate_use_acc = 64
    # Acceleromer is sampled by 32 hz so we do downsampling again
    if use_acc:
        sampling_rate_use_acc = 32  # 32 while using acc,  64 without
        e4_ppg = downsample(e4_ppg, acc_e4_df, 64 / sampling_rate_use_acc)
        ppg_down = downsample(ppg_down, acc_e4_df, 64 / sampling_rate_use_acc)
        acc_e4_df = pd.DataFrame(data=avg, index=range(1, len(avg) + 1), columns=['acc'])

    # Filtering by heartpy algorithm
    filtered_e4_3 = hp.filter_signal(e4_ppg[0].values, cutoff=[0.5, 1.5],
                                     sample_rate=sampling_rate_use_acc, order=3, filtertype='bandpass', return_top=False)

    filtered_biopac_ppg = hp.filter_signal(ppg_down[0].values, cutoff=[0.5, 1.5],
                                           sample_rate=sampling_rate_use_acc, order=3, filtertype='bandpass', return_top=False)

    signals, e4_peaks = nk.ppg_process(filtered_e4_3, sampling_rate=sampling_rate_use_acc)
    new_peaks_e4, clean_factor = filter_peaks(e4_peaks["PPG_Peaks"])

    ### Using HRV methotds to filter rri's. Using threshlod filter to replace undetected peaks with median


   # print(feature_filter_rri(rri_e4,filtered_rri_e4))

    signals_2, biopac_ppg_peaks = nk.ppg_process(filtered_biopac_ppg, sampling_rate=sampling_rate_use_acc)
    new_peaks_bp_ppg, clean_factor = filter_peaks(biopac_ppg_peaks["PPG_Peaks"])

    info, biopac_ecg_peaks = nk.ecg_peaks(ecg_down[0].to_numpy(), sampling_rate=sampling_rate_use_acc)
    new_peaks_bp_ecg, clean_factor = filter_peaks(biopac_ecg_peaks["ECG_R_Peaks"])


    return pd.concat([pd.DataFrame({"E4_PPG": new_peaks_e4}), pd.DataFrame({"BP_PPG": new_peaks_bp_ppg}), pd.DataFrame({"BP_ECG": new_peaks_bp_ecg})], axis=1)



def saveToFile(subject, channel, hrv):
    hrv.set_axis([subject], inplace=True)
    if subject == 'A001' or subject == 'EP02_1593189235.csv':
        hrv.to_csv("results\hrv_" + channel + ".csv", mode='w', header=True, index=True)
    else:
        hrv.to_csv("results\hrv_" + channel + ".csv", mode='a', header=False, index=True)
