# -*- coding: UTF-8 -*-
#preprocessing.py
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

## Run settings, set first
run_rest = True; #tbd
run_ert = True; #tbd
plot_hrv = False;
sampling_rate_use_acc=64

path_of_subjects = 'C:\\Users\\Przemek\\Documents\\PILTOR_bysubject\\'

distored_signal_subjects = ['A018','A022','A024','A027','A042','A043','A054','A056','A057','A058','A059','A062','A066','A069','A070','A072'
                            ,'A077' ,'A078','A079','A080','A081','A082','A083','A085','A086','A087','A088','A090','A091',
                           'A093','A097','A104','A105','A108','A121','A122','A123'] #123 lack of acq file 77 - 91 missings in the empatica 58 59 62 70 check a39 - cos dziwnego
subjects = [name for name in listdir(path_of_subjects) if name not in distored_signal_subjects and name[0] == "A"]
subjects = subjects[66:]

print(subjects)


for subject in subjects:

    print('subject ' + subject)
    if run_rest:
        print ("REST")
        bp_ppg_rest = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ppg\\" + subject + "_rest_ppg.csv",
            sep=" ", header=None)
        e4_ppg_rest = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_rest_e4_ppg.csv",
            sep=" ", header=None)
        bp_ecg_rest = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ecg\\" + subject + "_rest_ecg.csv",
            sep=" ", header=None)

        # accelometer
        acc_e4 = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_rest_e4_acc.csv",
            sep=",", header=None)

        acc_x = acc_e4.loc[2:, 0].reset_index(drop=True);
        acc_y = acc_e4.loc[2:, 1].reset_index(drop=True);
        acc_z = acc_e4.loc[2:, 2].reset_index(drop=True);

        acc_x_tmp = acc_x[0]
        acc_y_tmp = acc_y[0]
        acc_z_tmp = acc_z[0]
        avg_sec = 0

        avg = np.zeros(len(acc_x))

        # Empatica algorithm for calculating avarrege of acceleration from X, Y, Z factors
        for n in range(len(acc_x)):
            if (n % 32 == 0):
                avg[n] = avg[n - 1] * 0.9 + (avg_sec / 32) * 0.1
                avg_sec = 0
            elif (n >= 1):
                avg[n] = avg[n - 1]
            avg_sec += max((acc_x[n] - acc_x_tmp), (acc_y[n] - acc_y_tmp), (acc_z[n] - acc_z_tmp), key=abs)
            acc_x_tmp = acc_x[n]
            acc_y_tmp = acc_y[n]
            acc_z_tmp = acc_z[n]

       # plt.plot(avg)
      #  plt.show()
       # print(stat.stdev(avg))


        peaks = filter_signals(bp_ecg_rest, bp_ppg_rest, e4_ppg_rest)

        if plot_hrv:
            hrv1 = nk.hrv(peaks["E4_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc,  show=True)
            plt.show()
        else:
            hrv1 = nk.hrv(peaks["E4_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc,  show=False)

        saveToFile(subject,"rest_e4_ppg",hrv1)

        if plot_hrv:
            hrv2 = nk.hrv(peaks["BP_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=True)
            plt.show()
        else:
            hrv2 = nk.hrv(peaks["BP_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)

        saveToFile(subject,"rest_biopac_ppg",hrv2)
        if plot_hrv:
            hrv3 = nk.hrv(peaks["BP_ECG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=True)
            plt.show()
        else:
            hrv3 = nk.hrv(peaks["BP_ECG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)

        saveToFile(subject,"rest_biopac_ecg",hrv3)

    if run_ert:
        # prefeedback
        print("PREFEEDBACK")
        bp_ppg_prefeedback = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ppg\\" + subject + "_ppg_prefeedback.csv",
            sep=" ", header=None)
        e4_ppg_prefeedback = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_prefeedback_e4_ppg.csv",
            sep=" ", header=None)
        bp_ecg_prefeedback = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ecg\\" + subject + "_ecg_prefeedback_ecg.csv",
            sep=" ", header=None)


        peaks_prefeedback = filter_signals(bp_ecg_prefeedback, bp_ppg_prefeedback, e4_ppg_prefeedback)

        hrv1 = nk.hrv(peaks_prefeedback["E4_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"prefeedback_e4_ppg",hrv1)

        hrv2 = nk.hrv(peaks_prefeedback["BP_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"prefeedback_biopac_ppg",hrv2)

        hrv3 = nk.hrv(peaks_prefeedback["BP_ECG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"prefeedback_biopac_ecg",hrv3)

        # postfeedback
        print("POSTFEEDBACK")
        bp_ppg_postfeedback = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ppg\\" + subject + "_ppg_postfeedback.csv",
            sep=" ", header=None)
        e4_ppg_postfeedback = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_postfeedback_e4_ppg.csv",
            sep=" ", header=None)
        bp_ecg_postfeedback = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ecg\\" + subject + "_ecg_postfeedback_ecg.csv",
            sep=" ", header=None)


        peaks_postfeedback = filter_signals(bp_ecg_prefeedback, bp_ppg_postfeedback, e4_ppg_prefeedback)

        hrv1 = nk.hrv(peaks_postfeedback["E4_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"postfeedback_e4_ppg",hrv1)

        hrv2 = nk.hrv(peaks_postfeedback["BP_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"postfeedback_biopac_ppg",hrv2)

        hrv3 = nk.hrv(peaks_postfeedback["BP_ECG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"postfeedback_biopac_ecg",hrv3)

        #observe1
        print("OBSERVE_1")
        bp_ecg_observe1 = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ecg\\" + subject + "_ecg_observe1.csv",
            sep=" ", header=None)
        bp_ppg_observe1 = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ppg\\" + subject + "_ppg_observe1.csv",
            sep=" ", header=None)
        e4_ppg_observe1 = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_observe1_e4_ppg.csv",
            sep=" ", header=None)


        peaks_observe1 = filter_signals(bp_ecg_observe1, bp_ppg_observe1, e4_ppg_observe1)

        hrv1 = nk.hrv(peaks_observe1["E4_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"observe1_e4_ppg",hrv1)

        hrv2 = nk.hrv(peaks_observe1["BP_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"observe1_biopac_ppg",hrv2)

        hrv3 = nk.hrv(peaks_observe1["BP_ECG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"observe1_biopac_ecg",hrv3)

        #observe2
        print("OBSERVE_2")

        bp_ecg_observe2 = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ecg\\" + subject + "_ecg_observe2.csv",
            sep=" ", header=None)
        bp_ppg_observe2 = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ppg\\" + subject + "_ppg_observe2.csv",
            sep=" ", header=None)
        e4_ppg_observe2 = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_observe2_e4_ppg.csv",
            sep=" ", header=None)


        peaks_observe2 = filter_signals(bp_ecg_observe2, bp_ppg_observe2, e4_ppg_observe2)

        hrv1 = nk.hrv(peaks_observe2["E4_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"observe2_e4_ppg",hrv1)

        hrv2 = nk.hrv(peaks_observe2["BP_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"observe2_biopac_ppg",hrv2)

        hrv3 = nk.hrv(peaks_observe2["BP_ECG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"observe2_biopac_ecg",hrv3)


        #increase
        print("INCREASE")

        bp_ecg_increase = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ecg\\" + subject + "_ecg_increase.csv",
            sep=" ", header=None)
        bp_ppg_increase = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ppg\\" + subject + "_ppg_increase.csv",
            sep=" ", header=None)
        e4_ppg_increase = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_increase_e4_ppg.csv",
            sep=" ", header=None)


        peaks_increase = filter_signals(bp_ecg_increase, bp_ppg_increase, e4_ppg_increase)

        hrv1 = nk.hrv(peaks_increase["E4_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"increase_e4_ppg",hrv1)

        hrv2 = nk.hrv(peaks_increase["BP_PPG"].dropna().drop_duplicates().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"increase_biopac_ppg",hrv2)

        hrv3 = nk.hrv(peaks_increase["BP_ECG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"increase_biopac_ecg",hrv3)


        #decrease
        print("DECREASE")

        bp_ecg_decrease = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ecg\\" + subject + "_ecg_decrease.csv",
            sep=" ", header=None)
        bp_ppg_decrease = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_ppg\\" + subject + "_ppg_decrease.csv",
            sep=" ", header=None)
        e4_ppg_decrease = pd.read_csv(
            r"C:\Users\Przemek\Documents\PILTOR_bysubject\\" + subject + r"\data_processing_e4\\" + subject + "_decrease_e4_ppg.csv",
            sep=" ", header=None)

        peaks_decrease = filter_signals(bp_ecg_decrease, bp_ppg_decrease, e4_ppg_decrease)

        hrv1 = nk.hrv(peaks_decrease["E4_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"decrease_e4_ppg",hrv1)

        hrv2 = nk.hrv(peaks_decrease["BP_PPG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"decrease_biopac_ppg",hrv2)

        hrv3 = nk.hrv(peaks_decrease["BP_ECG"].dropna().drop_duplicates(), sampling_rate=sampling_rate_use_acc, show=False)  # show=True)
        saveToFile(subject,"decrease_biopac_ecg",hrv3)


    #
    # # Acceleromer is sampled by 32 hz so we do downsampling again
    #
    #     acc_e4_df = pd.DataFrame(data=avg, index=range(1, len(avg) + 1), columns=['acc'])
    #     e4_down_acc = downsample(empatica_ppg, acc_e4_df, 64/sampling_rate_use_acc)
    #     ppg_down_acc_biopac = downsample(ppg_down, acc_e4_df, 64/sampling_rate_use_acc)
    #
    #     # Filtering by heartpy algorithm
    #     filtered_e4_3 = hp.filter_signal(e4_down_acc[0].values.tolist(), cutoff=[0.5, 1.5], sample_rate = sampling_rate_use_acc,  order=3,
    #                                            filtertype='bandpass', return_top = False)
    #
    #     filtered_biopac_ppg = hp.filter_signal(ppg_down_acc_biopac[0].values.tolist(), cutoff=[0.5, 1.5], sample_rate = sampling_rate_use_acc, order=3,
    #                                            filtertype='bandpass', return_top = False)
    #
    #     signals, e4_peaks = nk.ppg_process(filtered_e4_3, sampling_rate=sampling_rate_use_acc)
    #
    #     ### Using HRV methotds to filter rri's. Using threshlod filter to replace undetected peaks with median
    #
    #     rri_e4 = np.diff(e4_peaks['PPG_Peaks'])/sampling_rate_use_acc * 1000
    #     ## The calculated rri needs to be chop from both sides beacause it was causing errors unpredictable results for some subjects
    #     filt_rri_e4 = threshold_filter(rri_e4[1:(len(rri_e4) - 2)], threshold = 'strong', local_median_size = 5)
    #     new_peaks_e4 = rrToPeaks(filt_rri_e4, sampling_rate_use_acc)
    #     new_peaks_e4 = new_peaks_e4 + e4_peaks['PPG_Peaks'][0]
    #     e4_peaks['PPG_Peaks'] = new_peaks_e4
    #
    #     # # Visualize the processing
    #  ##   plt.show()

    #
    #     signals_2, biopac_peaks = nk.ppg_process(filtered_biopac_ppg, sampling_rate=sampling_rate_use_acc)


      #  rri_biopac = np.diff(biopac_peaks['PPG_Peaks'])/32 *1000
        ## The calculated rri needs to be chop from both sides beacause it was causing errors unpredictable results for some subjects
        #print(rri_biopac)
       # filt_rri_biopac = threshold_filter(rri_biopac[1:(len(rri_biopac) - 2)], threshold='strong', local_median_size=5)
       # new_peaks_biopac = rrToPeaks(filt_rri_biopac, 32)
       # new_peaks_biopac = new_peaks_biopac + biopac_peaks['PPG_Peaks'][0]
       # biopac_peaks['PPG_Peaks'] = new_peaks_biopac filter
        #fig = plt.gcf()
        #fig.canvas.set_window_title('Biopac ECG')
        #print(hrv3)






# DECOMPOSITION OF SIGNAL
        #   #  fig = plt.gcf()
        #    # fig.canvas.set_window_title('E4')
        #     #possible decomposition of a signal. Not very spectecular effect
        #     #components = nk.signal_decompose(filtered_e4, method='emd')
        #     #nk.signal_plot(components)  # Visualize components
        #     #plt.plot(components);
        #     #plt.show()

