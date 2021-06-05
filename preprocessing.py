# -*- coding: UTF-8 -*-
import bioread
import os.path
import shutil
from os import listdir

import heartpy as hp
import neurokit2 as nk
from hrv.filters import moving_median
from hrv.filters import moving_median
from hrv.filters import threshold_filter

import mat4py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
from matplotlib import pyplot
from scipy import signal

## Run settings, set first
run_rest = True; #tbd
run_ert = False; #tbd

path_of_subjects = 'C:\\Users\\Przemek\\Documents\\PILTOR_bysubject\\'

distored_signal_subjects = ['A018','A022','A027','A042','A054','A056','A066','A069','A072','A077','A083','A091','A093','A097','A105','A108','A121','A122']

subjects = [name for name in listdir(path_of_subjects) if name not in distored_signal_subjects and name[0] == "A"]
# subjects =
print(subjects)
def downsample(base_signal, signal_under, f_rate):
    indexes = (signal_under.index * f_rate).astype(int)
    indexes = indexes[indexes < base_signal.shape[0]]
    return base_signal.iloc[indexes].reset_index(drop=True)


def printHead(signal):
    plt.plot(signal.head(1000))

def rrToPeaks(rri, sampl):
    res = [0]
    for n, rr in enumerate(rri):
        if n > 0:
            res.append((res[n - 1] + ((rr / 1000) * sampl)).astype(int))
    return res

for subject in subjects:
    print('subject ' + subject)

    ppg = pd.read_csv(
        r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ppg\\" + subject + "_rest_ppg.csv",
        sep=" ", header=None)
    empatica_ppg = pd.read_csv(
        r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_rest_e4_ppg.csv",
        sep=" ", header=None)
    ecg = pd.read_csv(
        r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ecg\\" + subject + "_rest_ecg.csv",
        sep=" ", header=None)


    # accelometer
    acc_e4 = pd.read_csv(
        r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_rest_e4_acc.csv",
        sep=",", header=None)

    acc_x = acc_e4.loc[2:, 0].reset_index(drop=True);
    acc_y = acc_e4.loc[2:, 1].reset_index(drop=True);
    acc_z = acc_e4.loc[2:, 2].reset_index(drop=True);

    acc_x_tmp = 0
    acc_y_tmp = 0
    acc_z_tmp = 0
    avg_sec = 0

    avg = np.zeros(len(acc_x))

    # Empatica algorithm for calculating avarrege of acceleration from X, Y, Z factors
    for n in range(len(acc_x)):
        if (n % 32 == 0):
            avg[n] = avg[n - 1] * 0.9 + (avg_sec / 32) * 0.1
            avg_sec = 0
        elif (n > 1):
            avg[n] = avg[n - 1]
        avg_sec += max((acc_x[n] - acc_x_tmp), (acc_y[n] - acc_y_tmp), (acc_z[n] - acc_z_tmp), key=abs)
        acc_x_tmp = acc_x[n]
        acc_y_tmp = acc_y[n]
        acc_z_tmp = acc_z[n]



    # Sampling rate to 2000 Hz for biopac and 64 Hz for E4. Downsampling...
    #ecg_down = downsample(ecg, empatica_ppg, 2000 / 64)
    ppg_down = downsample(ppg, empatica_ppg, 2000 / 64)
    #print(ecg_down.shape)


# Acceleromer is sampled by 32 hz so we do downsampling again

    sampling_rate_use_acc = 64 # 32 while using acc,  64 without

# correlation seems to be a little bit absurd feature diffrence to feature score
    acc_e4_df = pd.DataFrame(data=avg, index=range(1, len(avg) + 1), columns=['acc'])
    e4_down_acc = downsample(empatica_ppg, acc_e4_df, 64/sampling_rate_use_acc)
    ppg_down_acc_biopac = downsample(ppg_down, acc_e4_df, 64/sampling_rate_use_acc)

    # Filtering by heartpy algorithm
    filtered_e4_3 = hp.filter_signal(e4_down_acc[0].values.tolist(), cutoff=[0.5, 1.5], sample_rate = sampling_rate_use_acc,  order=3,
                                           filtertype='bandpass', return_top = False)

    filtered_biopac_ppg = hp.filter_signal(ppg_down_acc_biopac[0].values.tolist(), cutoff=[0.5, 1.5], sample_rate = sampling_rate_use_acc, order=3,
                                           filtertype='bandpass', return_top = False)

    signals, e4_peaks = nk.ppg_process(filtered_e4_3, sampling_rate=sampling_rate_use_acc)

    ### Using HRV methotds to filter rri's. Using threshlod filter to replace undetected peaks with median

    rri_e4 = np.diff(e4_peaks['PPG_Peaks'])/sampling_rate_use_acc * 1000
    ## The calculated rri needs to be chop from both sides beacause it was causing errors unpredictable results for some subjects
    filt_rri_e4 = threshold_filter(rri_e4[1:(len(rri_e4) - 2)], threshold='strong', local_median_size=5)
    new_peaks_e4 = rrToPeaks(filt_rri_e4, sampling_rate_use_acc)
    new_peaks_e4 = new_peaks_e4 + e4_peaks['PPG_Peaks'][0]
    e4_peaks['PPG_Peaks'] = new_peaks_e4

    # # Visualize the processing
    plt.show()
    hrv1 = nk.hrv(e4_peaks, sampling_rate=sampling_rate_use_acc, show=True)
    fig = plt.gcf()
    fig.canvas.set_window_title('E4')
    #possible decomposition of a signal. Not very spectecular effect
    #components = nk.signal_decompose(filtered_e4, method='emd')
    #nk.signal_plot(components)  # Visualize components
    #plt.plot(components);
    #plt.show()

    pd.set_option('display.max_rows',52)
    hrv1.to_csv("hrv_e4.txt")
    print(hrv1)

    signals_2, biopac_peaks = nk.ppg_process(filtered_biopac_ppg, sampling_rate=sampling_rate_use_acc)


  #  rri_biopac = np.diff(biopac_peaks['PPG_Peaks'])/32 *1000
    ## The calculated rri needs to be chop from both sides beacause it was causing errors unpredictable results for some subjects
    #print(rri_biopac)
   # filt_rri_biopac = threshold_filter(rri_biopac[1:(len(rri_biopac) - 2)], threshold='strong', local_median_size=5)
   # new_peaks_biopac = rrToPeaks(filt_rri_biopac, 32)
   # new_peaks_biopac = new_peaks_biopac + biopac_peaks['PPG_Peaks'][0]
   # biopac_peaks['PPG_Peaks'] = new_peaks_biopac filter

    hrv2 = nk.hrv(biopac_peaks, sampling_rate=sampling_rate_use_acc, show=True)
    fig = plt.gcf()
    fig.canvas.set_window_title('Biopac PPG')
    print(hrv2)
    hrv2.to_csv("hrv_biopac.txt")

    ecg_biopac_peaks, info = nk.ecg_peaks(ecg[0].to_numpy(), sampling_rate=2000)
    hrv3 = nk.hrv(ecg_biopac_peaks, sampling_rate=2000, show=True)
    fig = plt.gcf()
    fig.canvas.set_window_title('Biopac ECG')
    print(hrv3)

