# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 23:16:51 2018

Script to compare the individual courves of the differenct frequencies

@author: Roman
"""

from __future__ import division
import os
#from pyface.qt import QtGui, QtCore
import mne # import the processing tools
#import os.path as op
import numpy as np
import matplotlib.pyplot as plt
import math
import shutil
from itertools import repeat
from mne.time_frequency import tfr_morlet, psd_multitaper # for freq and time-freq analyses

import seaborn as sns
sns.axes_style('darkgrid')

parent_dir = 'D:/CHF'
base_dir = 'D:/CHF/processing_new_all'
data_dir = base_dir + '/group/data_arrays/'

os.chdir(base_dir)

stimdurations=[0.546,0.278,0.144,0.078]
stimdurations_ms=[int(i*1000) for i in stimdurations]
#############################
### load arrays #############
#############################

full2hz = np.load(data_dir + os.sep + 'zoom_2_Hz_full_electrode_Pz_reference [TP9, TP10].npy')
full4hz = np.load(data_dir + os.sep + 'zoom_4_Hz_full_electrode_Pz_reference [TP9, TP10].npy')
full8hz = np.load(data_dir + os.sep + 'zoom_8_Hz_full_electrode_Pz_reference [TP9, TP10].npy')
full15hz = np.load(data_dir + os.sep + 'zoom_15_Hz_full_electrode_Pz_reference [TP9, TP10].npy')

upper2hz = np.load(data_dir + os.sep + 'zoom_2_Hz_upper_electrode_Pz_reference [TP9, TP10].npy')
upper4hz = np.load(data_dir + os.sep + 'zoom_4_Hz_upper_electrode_Pz_reference [TP9, TP10].npy')
upper8hz = np.load(data_dir + os.sep + 'zoom_8_Hz_upper_electrode_Pz_reference [TP9, TP10].npy')
upper15hz = np.load(data_dir + os.sep + 'zoom_15_Hz_upper_electrode_Pz_reference [TP9, TP10].npy')

lower2hz = np.load(data_dir + os.sep + 'zoom_2_Hz_lower_electrode_Pz_reference [TP9, TP10].npy')
lower4hz = np.load(data_dir + os.sep + 'zoom_4_Hz_lower_electrode_Pz_reference [TP9, TP10].npy')
lower8hz = np.load(data_dir + os.sep + 'zoom_8_Hz_lower_electrode_Pz_reference [TP9, TP10].npy')
lower15hz = np.load(data_dir + os.sep + 'zoom_15_Hz_lower_electrode_Pz_reference [TP9, TP10].npy')

sum2hz = np.load(data_dir + os.sep + 'zoom_2_Hz_sum_electrode_Pz_reference [TP9, TP10].npy')
sum4hz = np.load(data_dir + os.sep + 'zoom_4_Hz_sum_electrode_Pz_reference [TP9, TP10].npy')
sum8hz = np.load(data_dir + os.sep + 'zoom_8_Hz_sum_electrode_Pz_reference [TP9, TP10].npy')
sum15hz = np.load(data_dir + os.sep + 'zoom_15_Hz_sum_electrode_Pz_reference [TP9, TP10].npy')

diff2hz = np.load(data_dir + os.sep + 'zoom_DIFF_2_Hz_electrode_Pz_reference [TP9, TP10].npy')
diff4hz = np.load(data_dir + os.sep + 'zoom_DIFF_4_Hz_electrode_Pz_reference [TP9, TP10].npy')
diff8hz = np.load(data_dir + os.sep + 'zoom_DIFF_8_Hz_electrode_Pz_reference [TP9, TP10].npy')
diff15hz = np.load(data_dir + os.sep + 'zoom_DIFF_15_Hz_electrode_Pz_reference [TP9, TP10].npy')

######## DIFF CURVES #######

### 8 Hz
hypo_diff_8hz_pred_by_2hz = np.ndarray((4,546)) # 10 ist einfach so
for i in range(4):
    hypo_diff_8hz_pred_by_2hz[i][:] = diff2hz
    # shift along the 2nd axis
    for j in range(i):
        hypo_diff_8hz_pred_by_2hz[i] = np.roll(hypo_diff_8hz_pred_by_2hz[i], stimdurations_ms[2], axis=0)
    
summe = np.sum(hypo_diff_8hz_pred_by_2hz, axis = 0)

f, ax = plt.subplots()
sns.lineplot(x = range(len(summe)), y = summe, ax=ax, color='gray')
sns.lineplot(x = range(len(np.tile(diff8hz,4))), y = np.tile(diff8hz,4), ax=ax, color='black')
plt.legend(['hypothetical summation','real time course'])
plt.show()

### 15 Hz
hypo_diff_15hz_pred_by_2hz = np.ndarray((7,546)) # 10 ist einfach so
for i in range(7):
    hypo_diff_15hz_pred_by_2hz[i][:] = diff2hz
    # shift along the 2nd axis
    for j in range(i):
        hypo_diff_15hz_pred_by_2hz[i] = np.roll(hypo_diff_15hz_pred_by_2hz[i], stimdurations_ms[3], axis=0)
    
summe = np.sum(hypo_diff_15hz_pred_by_2hz, axis = 0)

f, ax = plt.subplots()
sns.lineplot(summe, ax=ax, color='gray')
sns.lineplot(np.tile(diff15hz,7), ax=ax, color='black')
plt.legend(['hypothetical summation','real time course'])



######## COMMON TIME CURVES #######

for VF, field in zip(['FVF','UVF','LVF'],['full','upper','lower']):
    
    if field == 'full':
        v2hz = full2hz
        v4hz = full4hz
        v8hz = full8hz
        v15hz = full15hz
    elif field == 'upper':
        v2hz = upper2hz
        v4hz = upper4hz
        v8hz = upper8hz
        v15hz = upper15hz
    elif field == 'lower':
        v2hz = lower2hz
        v4hz = lower4hz
        v8hz = lower8hz
        v15hz = lower15hz
    
    ### 8 Hz
    hypo_8hz_pred_by_2hz = np.ndarray((4,546)) # 10 ist einfach so
    for i in range(4):
        hypo_8hz_pred_by_2hz[i][:] = v2hz
        # shift along the 2nd axis
        for j in range(i):
            hypo_8hz_pred_by_2hz[i] = np.roll(hypo_8hz_pred_by_2hz[i], stimdurations_ms[2], axis=0)
        
    summe = np.sum(hypo_8hz_pred_by_2hz, axis = 0)
    
    f, ax = plt.subplots()
    sns.lineplot(summe, ax=ax, color='gray')
    sns.lineplot(np.tile(v8hz,4), ax=ax, color='black')
    plt.legend(['hypothetical summation','real time course'])
    plt.title(VF + '_8_Hz_predicted_by_2_Hz')

    ### 15 Hz
    hypo_15hz_pred_by_2hz = np.ndarray((7,546)) # 10 ist einfach so
    for i in range(7):
        hypo_15hz_pred_by_2hz[i][:] = v2hz
        # shift along the 2nd axis
        for j in range(i):
            hypo_15hz_pred_by_2hz[i] = np.roll(hypo_15hz_pred_by_2hz[i], stimdurations_ms[3], axis=0)
        
    summe = np.sum(hypo_15hz_pred_by_2hz, axis = 0)
    
    f, ax = plt.subplots()
    sns.lineplot(summe, ax=ax, color='gray')
    sns.lineplot(np.tile(v15hz,7), ax=ax, color='black')
    plt.legend(['hypothetical summation','real time course'])
    plt.title(VF + '_15_Hz_predicted_by_2_Hz')
