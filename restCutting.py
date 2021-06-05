import os.path
import shutil
import glob
import neurokit
import pandas as pd
import math
from datetime import datetime

path_of_subjects = 'C:\\Users\\Przemek\\Documents\\PILTOR_bysubject\\'

distored_signal_subjects = ['A018','A022','A027','A042','A054','A056','A066','A069','A072','A077','A083','A091','A093','A097','A105','A108','A121','A122']

#subjects = [name for name in listdir(path_of_subjects) if name not in distored_signal_subjects and name[0] == "A"]
subjects = ['A001', 'A002', 'A003']
for subj in subjects:
    path = "C:\\Users\\Przemek\\Documents\\PILTOR_bysubject\\" + subj + "\\"
    os.chdir(path)
    acq_files = glob.glob("*.acq")

    if not os.path.exists('data_processing_ecg'):
        os.mkdir('data_processing_ecg')
    if not os.path.exists('data_processing_ppg'):
        os.mkdir('data_processing_ppg')
    if not os.path.exists('data_processing_eda'):
        os.mkdir('data_processing_eda')
    if not os.path.exists('data_processing_e4'):
        os.mkdir('data_processing_e4')

    print(path + acq_files[0])
    acq_data = neurokit.read_acqknowledge(acq_files[0], path, index='range')

    markers_rest = acq_data[acq_data.iloc[:, 6] == 10].index

    markers_observe = [idx for idx, number in enumerate(acq_data.iloc[:, 6]) if (number in [21, 22, 23, 24])]
    markers_increase = [idx for idx, number in enumerate(acq_data.iloc[:, 6]) if (number in [31, 32])]
    markers_decrease = [idx for idx, number in enumerate(acq_data.iloc[:, 6]) if (number in [41, 42])]

    print("index observe")
    print(markers_observe)

#### Biopac rest

    ecg_rest = acq_data.iloc[markers_rest[19]:markers_rest[20]]['ECG100C']
    ppg_rest = acq_data.iloc[markers_rest[19]:markers_rest[20]]['PPG100C']

    ecg_rest.to_csv("data_processing_ecg\\" + subj + '_rest_ecg.csv', index=False, header=False)
    ppg_rest.to_csv("data_processing_ppg\\" + subj + '_rest_ppg.csv', index=False, header=False)

#### Biopac pre- postfeedback

    ecg_prefeedback = acq_data.iloc[(markers_rest[19] + 600000):markers_rest[20]]['ECG100C'];
    ecg_postfeedback = acq_data.iloc[markers_rest[49]:markers_rest[50]]['ECG100C'];

    ecg_prefeedback.to_csv("data_processing_ecg\\" + subj + '_prefeedback_ecg.csv', index=False, header=False)
    ecg_postfeedback.to_csv("data_processing_ecg\\" + subj + '_postfeedback_ecg.csv', index=False, header=False)


#### Biopac observed

    ecg_observe1 = acq_data.iloc[markers_observe[0]:(markers_observe[math.ceil(len(markers_observe)/2)] + 15000)]['ECG100C']
    ecg_observe2 = acq_data.iloc[(markers_observe[math.ceil(len(markers_observe)/2)] + 15000):(markers_observe[len(markers_observe) - 1] + 15000)]['ECG100C']

    ppg_observe1 = acq_data.iloc[markers_observe[0]:(markers_observe[math.ceil(len(markers_observe)/2)] + 15000)]['PPG100C']
    ppg_observe2 = acq_data.iloc[(markers_observe[math.ceil(len(markers_observe)/2)] + 15000):(markers_observe[len(markers_observe) - 1] + 15000)]['PPG100C']

    ecg_observe1.to_csv("data_processing_ecg\\" + subj + '_ecg_observe1.csv', index=False, header=False)
    ecg_observe2.to_csv("data_processing_ppg\\" + subj + '_ecg_observe2.csv', index=False, header=False)

    ppg_observe1.to_csv("data_processing_ecg\\" + subj + '_ppg_observe1.csv', index=False, header=False)
    ppg_observe2.to_csv("data_processing_ppg\\" + subj + '_ppg_observe2.csv', index=False, header=False)

#### Biopac increased

    ecg_increase = acq_data.iloc[markers_increase[0]:(markers_increase[len(markers_increase) - 1] + 15000)]['ECG100C']
    ppg_increase = acq_data.iloc[markers_increase[0]:(markers_increase[len(markers_increase) - 1] + 15000)]['PPG100C']

    ecg_increase.to_csv("data_processing_ecg\\" + subj + '_ecg_increase.csv', index=False, header=False)
    ppg_increase.to_csv("data_processing_ppg\\" + subj + '_ppg_increase.csv', index=False, header=False)

#### Biopac decrease

    ecg_decrease = acq_data.iloc[markers_decrease[0]:(markers_decrease[len(markers_decrease) - 1] + 15000)]['ECG100C']
    ppg_decrease = acq_data.iloc[markers_decrease[0]:(markers_decrease[len(markers_decrease) - 1] + 15000)]['PPG100C']

    ecg_decrease.to_csv("data_processing_ecg\\" + subj + '_ecg_decrease.csv', index=False, header=False)
    ppg_decrease.to_csv("data_processing_ppg\\" + subj + '_ppg_decrease.csv', index=False, header=False)

#### Empatica rest

    empatica_directory = 'empatica';
    ppg_e4 = pd.read_csv(path + r"\\empatica\\BVP.csv", sep=" ", header=None)
    eda_e4 = pd.read_csv(path + r"\\empatica\\EDA.csv", sep=" ", header=None)
    timestamp_e4 = pd.read_csv(path + r"\\empatica\\tags.csv", sep=" ", header=None)
    acc_e4 = pd.read_csv(path + r"\\empatica\\ACC.csv", sep=",", header=None)

    timestamp_acq = acq_data[acq_data.iloc[:, 6] == 20].index[0]
    seconds_before_stamp = timestamp_e4.iloc[0] - eda_e4.iloc[0];
    print(seconds_before_stamp[0])

    rest_e4_eda_onset = int(round((seconds_before_stamp[0] * 4) + (markers_rest[19] - timestamp_acq) / 500))
    rest_e4_eda_offset = int(round((seconds_before_stamp[0] * 4) + (markers_rest[20] - timestamp_acq) / 500))
    e4_rest_eda = eda_e4.iloc[rest_e4_eda_onset:rest_e4_eda_offset].to_csv("data_processing_e4\\" + subj + '_rest_e4_eda.csv', index=False, header=False)

    rest_e4_ppg_onset = int(round((seconds_before_stamp[0]) * 64 + (markers_rest[19] - timestamp_acq) * 0.032))
    rest_e4_ppg_offset = int(round((seconds_before_stamp[0]) * 64 + (markers_rest[20] - timestamp_acq) * 0.032))
    ppg_e4[rest_e4_ppg_onset:rest_e4_ppg_offset].to_csv("data_processing_e4\\" + subj + '_rest_e4_ppg.csv', index=False, header=False)

    rest_e4_acc_onset = int(round((seconds_before_stamp[0]) * 32 + (markers_rest[19] - timestamp_acq) * 0.016))
    rest_e4_acc_offset = int(round((seconds_before_stamp[0]) * 32 + (markers_rest[20] - timestamp_acq) * 0.016))
    e4_rest_acc = acc_e4[rest_e4_acc_onset:rest_e4_acc_offset].to_csv("data_processing_e4\\" + subj + '_rest_e4_acc.csv', index=False, header=False)

#### E4 pre- postfeedback

    prefeedback_e4_eda_onset = int(round((seconds_before_stamp[0] * 4) + ((markers_rest[19] - timestamp_acq) + 600000) * 0.002))
    prefeedback_e4_ede_offset = int(round((seconds_before_stamp[0] * 4) + (markers_rest[20] - timestamp_acq) * 0.002))
    eda_e4[prefeedback_e4_eda_onset:prefeedback_e4_ede_offset].to_csv("data_processing_e4\\" + subj + '_prefeedback_e4_eda.csv', index=False, header=False)


    prefeedback_e4_ppg_onset = int(round((seconds_before_stamp[0] * 64) + ((markers_rest[19] - timestamp_acq) + 600000) * 0.032))
    prefeedback_e4_ppg_offset = int(round((seconds_before_stamp[0] * 64) + (markers_rest[20] - timestamp_acq) * 0.032))
    ppg_e4[prefeedback_e4_ppg_onset:prefeedback_e4_ppg_offset].to_csv("data_processing_e4\\" + subj + '_prefeedback_e4_ppg.csv', index=False, header=False)


    prefeedback_e4_acc_onset = int(round((seconds_before_stamp[0] * 32) + ((markers_rest[19] - timestamp_acq) + 600000) * 0.016))
    prefeedback_e4_acc_offset = int(round((seconds_before_stamp[0] * 32) + (markers_rest[20] - timestamp_acq) * 0.016))
    acc_e4[prefeedback_e4_acc_onset:prefeedback_e4_acc_offset].to_csv("data_processing_e4\\" + subj + '_prefeedback_e4_acc.csv', index=False, header=False)


    postfeedback_e4_eda_onset = int(round((seconds_before_stamp[0] * 4) + (markers_rest[49] - timestamp_acq) * 0.002))
    postfeedback_e4_ede_offset = int(round((seconds_before_stamp[0] * 4) + (markers_rest[50] - timestamp_acq) * 0.002))
    eda_e4[postfeedback_e4_eda_onset:postfeedback_e4_ede_offset].to_csv("data_processing_e4\\" + subj + '_prefeedback_e4_eda.csv', index=False, header=False)


    postfeedback_e4_ppg_onset = int(round((seconds_before_stamp[0] * 64) + (markers_rest[49] - timestamp_acq) * 0.032))
    postfeedback_e4_ppg_offset = int(round((seconds_before_stamp[0] * 64) + (markers_rest[50] - timestamp_acq) * 0.032))
    ppg_e4[postfeedback_e4_ppg_onset:postfeedback_e4_ppg_offset].to_csv("data_processing_e4\\" + subj + '_prefeedback_e4_ppg.csv', index=False, header=False)


    postfeedback_e4_acc_onset = int(round((seconds_before_stamp[0] * 32) + (markers_rest[49] - timestamp_acq) * 0.016))
    postfeedback_e4_acc_offset = int(round((seconds_before_stamp[0] * 32) + (markers_rest[50] - timestamp_acq) * 0.016))
    acc_e4[postfeedback_e4_acc_onset:postfeedback_e4_acc_offset].to_csv("data_processing_e4\\" + subj + '_prefeedback_e4_acc.csv', index=False, header=False)


### E4 - observe 1 and 2

    observe1_e4_eda_onset = int(round((seconds_before_stamp[0] * 4) + (markers_observe[0] - timestamp_acq) * 0.002))
    observe1_e4_eda_offset = int(round((seconds_before_stamp[0] * 4) + (markers_observe[math.ceil(len(markers_observe)/2)] + 15000 - timestamp_acq) * 0.002))

    observe2_e4_eda_onset = int(round((seconds_before_stamp[0] * 4) + (markers_observe[math.ceil(len(markers_observe) / 2)] + 15000 - timestamp_acq) * 0.002))
    observe2_e4_eda_offset = int(round((seconds_before_stamp[0] * 4) + (markers_observe[len(markers_observe) - 1] + 15000) * 0.002))

    eda_e4[observe1_e4_eda_onset:observe1_e4_eda_offset].to_csv("data_processing_e4\\" + subj + '_observe1_e4_eda.csv', index=False, header=False)
    eda_e4[observe2_e4_eda_onset:observe2_e4_eda_offset].to_csv("data_processing_e4\\" + subj + '_observe2_e4_eda.csv', index=False, header=False)


    observe1_e4_ppg_onset = int(round((seconds_before_stamp[0] * 64) + (markers_observe[0] - timestamp_acq) * 0.032))
    observe1_e4_ppg_offset = int(round((seconds_before_stamp[0] * 64) + (markers_observe[math.ceil(len(markers_observe)/2)] + 15000 - timestamp_acq) * 0.032))

    observe2_e4_ppg_onset = int(round((seconds_before_stamp[0] * 64) + (markers_observe[math.ceil(len(markers_observe) / 2)] + 15000 - timestamp_acq) * 0.032))
    observe2_e4_ppg_offset = int(round((seconds_before_stamp[0] * 64) + (markers_observe[len(markers_observe) - 1] + 15000) * 0.032))

    ppg_e4[observe1_e4_ppg_onset:observe1_e4_ppg_offset].to_csv("data_processing_e4\\" + subj + '_observe1_e4_ppg.csv', index=False, header=False)
    ppg_e4[observe2_e4_ppg_onset:observe2_e4_ppg_offset].to_csv("data_processing_e4\\" + subj + '_observe2_e4_ppg.csv', index=False, header=False)


    observe1_e4_acc_onset = int(round((seconds_before_stamp[0] * 32) + (markers_observe[0] - timestamp_acq) * 0.016))
    observe1_e4_acc_offset = int(round((seconds_before_stamp[0] * 32) + (markers_observe[math.ceil(len(markers_observe)/2)] + 15000 - timestamp_acq) * 0.016))

    observe2_e4_acc_onset = int(round((seconds_before_stamp[0] * 32) + (markers_observe[math.ceil(len(markers_observe) / 2)] + 15000 - timestamp_acq) * 0.016))
    observe2_e4_acc_offset = int(round((seconds_before_stamp[0] * 32) + (markers_observe[len(markers_observe) - 1] + 15000) * 0.016))

    acc_e4[observe1_e4_acc_onset:observe1_e4_acc_offset].to_csv("data_processing_e4\\" + subj + '_observe1_e4_acc.csv', index=False, header=False)
    acc_e4[observe2_e4_acc_onset:observe2_e4_acc_offset].to_csv("data_processing_e4\\" + subj + '_observe2_e4_acc.csv', index=False, header=False)

#### E4 - increase

    increase_e4_eda_onset = int(round((seconds_before_stamp[0] * 4) + (markers_increase[0] - timestamp_acq) * 0.002))
    increase_e4_eda_offset = int(round((seconds_before_stamp[0] * 4) + (markers_increase[len(markers_increase) - 1] + 15000 - timestamp_acq) * 0.002))
    eda_e4[increase_e4_eda_onset:increase_e4_eda_offset].to_csv("data_processing_e4\\" + subj + '_increase_e4_eda.csv', index=False, header=False)

    increase_e4_ppg_onset = int(round((seconds_before_stamp[0] * 64) + (markers_increase[0] - timestamp_acq) * 0.032))
    increase_e4_ppg_offset = int(round((seconds_before_stamp[0] * 64) + (markers_increase[len(markers_increase) - 1] + 15000 - timestamp_acq) * 0.032))
    ppg_e4[increase_e4_ppg_onset:increase_e4_ppg_offset].to_csv("data_processing_e4\\" + subj + '_increase_e4_ppg.csv', index=False, header=False)

    increase_e4_acc_onset = int(round((seconds_before_stamp[0] * 32) + (markers_increase[0] - timestamp_acq) * 0.016))
    increase_e4_acc_offset = int(round((seconds_before_stamp[0] * 32) + (markers_increase[len(markers_increase) - 1] + 15000 - timestamp_acq) * 0.016))
    acc_e4[increase_e4_acc_onset:increase_e4_acc_offset].to_csv("data_processing_e4\\" + subj + '_increase_e4_acc.csv', index=False, header=False)

#### E4 - decrease

    decrease_e4_eda_onset = int(round((seconds_before_stamp[0] * 4) + (markers_decrease[0] - timestamp_acq) * 0.002))
    decrease_e4_eda_offset = int(round((seconds_before_stamp[0] * 4) + (markers_decrease[len(markers_decrease) - 1] + 15000 - timestamp_acq) * 0.002))
    eda_e4[decrease_e4_eda_onset:decrease_e4_eda_offset].to_csv("data_processing_e4\\" + subj + '_decrease_e4_eda.csv', index=False, header=False)

    decrease_e4_ppg_onset = int(round((seconds_before_stamp[0] * 64) + (markers_decrease[0] - timestamp_acq) * 0.032))
    decrease_e4_ppg_offset = int(round((seconds_before_stamp[0] * 64) + (markers_decrease[len(markers_decrease) - 1] + 15000 - timestamp_acq) * 0.032))
    ppg_e4[decrease_e4_ppg_onset:decrease_e4_ppg_offset].to_csv("data_processing_e4\\" + subj + '_decrease_e4_ppg.csv', index=False, header=False)

    decrease_e4_acc_onset = int(round((seconds_before_stamp[0] * 32) + (markers_decrease[0] - timestamp_acq) * 0.016))
    decrease_e4_acc_offset = int(round((seconds_before_stamp[0] * 32) + (markers_decrease[len(markers_decrease) - 1] + 15000 - timestamp_acq) * 0.016))
    acc_e4[decrease_e4_acc_onset:decrease_e4_acc_offset].to_csv("data_processing_e4\\" + subj + '_decrease_e4_acc.csv', index=False, header=False)
