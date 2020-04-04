#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pipeline for usage

@author: cth
"""

# define working directory

import os
import mne # import the processing tools
import numpy as np
import matplotlib.pyplot as plt
import math
import shutil
from itertools import repeat

parent_dir = 'D:/CHF/polarityReversal/'
base_dir = 'D:/CHF/polarityReversal/eeg/'
os.chdir(base_dir)

# import personal scriptsimport sys
#import sys
#sys.path.append(base_dir + '/functions')
#from rk_calculate_timings import *
#from rk_chf_processing_functions import *

# define the trigger lengths
screenfreq=90.0 #Hz Monitor frame rate
time_per_frame=1.0/screenfreq
Ncond=12 # number of conditions
stimfreq=[1.875,3.75,7.5,15]
stimfreq2=[i for i in stimfreq for _ in range(3)]

#Nframes_per_stim = [ screenfreq / x / 2 for x in stimfreq ] # pure stimulus duration, same duration counts for the ISI to come
#Nframes_per_ISI = [ x * 2 - min(Nframes_per_stim) for x in Nframes_per_stim ]
#time_per_stim_and_ISI = [ ((x + min(Nframes_per_stim)) / screenfreq) for x in Nframes_per_ISI ]
#time_per_stim_and_ISI2 = [x for item in time_per_stim_and_ISI for x in repeat(item, 3)]

stimdurations = open(parent_dir + "parameters/stimdurations.txt").readlines()
stimdurations=[0.546,0.278,0.144,0.078]
stimdurations2=[i for i in stimdurations for _ in range(3)]
hemifields=['full','lower','upper'] # this was upside down in the initial processing, but now changed upper to lower and vice versa, because wrong definition in paradigm found
hemifields2=[i for _ in range(3) for i in hemifields]

resolution_eeg=0.002 # temporal resolution of the eeg, equals to 1 / EEG-Frequency
min_triggerlength=[]
max_triggerlength=[]
for i in range(2,14): # den Trigger, der den einzelnen Stimulus onset definiert (1 frame), mal weggelassen
    min_triggerlength.append(math.floor((i*time_per_frame-resolution_eeg*1.5)*1000)/1000) #abrunden
    max_triggerlength.append(math.ceil((i*time_per_frame+resolution_eeg*1.5)*1000)/1000) #aufrunden

# these are valid for subjects 4-10
single_event_timings=[[322, 868, 1413, 1959],
                      [322, 601, 879, 1157, 1435, 1714, 1992, 2270],
                      [322, 467, 612, 756, 901, 1046, 1190, 1335, 1480, 1625, 1769, 1914, 2060, 2204, 2349],
                      [322, 400, 478, 556, 634, 712, 790, 868, 946, 1023, 1101, 1179, 1257, 1335, 1413, 1491, 1569, 1647, 1725, 1804, 1882, 1960, 2038, 2116, 2194, 2272, 2349, 2427, 2505, 2583]]
single_period_times= [546,279,145,78]

# turn information status messages on, if off, type: mne.set_log_level('WARNING')
mne.set_log_level('INFO')

colormapQ3a=[[166,206,227], [31,120,180], [178,223,138]]
colormapQ3a=[[j/255 for j in i] for i in colormapQ3a];

colormapQ4a=[(0,0,0),(51/255,160/255,44/255),(31/255,120/255,180/255),(0.5,0.5,0.5)]
#colormapQ4a=[[j/255 for j in i] for i in colormapQ4a];

cm       = colormapQ4a
tmin     = 0.
tmax     = 3.

for sub in range(4,11): # 4, 11
#for sub in range(1,11):
    if sub<10:
        subs = '0' + str(sub)
    else:
        subs = str(sub)
    subject='RK_200' + subs
    
    # delete old files
    shutil.rmtree(subject + os.sep + 'plots')
    shutil.rmtree(subject + os.sep + 'data')
    
    if not os.path.isdir(subject + os.sep + 'data'):
        os.makedirs(subject + os.sep + 'data')
    if not os.path.isdir(subject + os.sep + 'plots'):
        os.makedirs(subject + os.sep + 'plots')
    
    #for eoi in ['Oz', 'O1', 'O2', 'Pz']:
    for eoi in ['Pz']:    
        print(eoi)
                
        #for ref in [[],['Fz'],['Cz'],['TP9','TP10']]:
        for ref in [['TP9','TP10']]: # [['TP9','TP10']]:
            print(ref)
            #ref = ['Fz'] #Fz
            
            ###### PIPELINE - Time Analysis ############
            logfile = open('n_blinks.txt', 'a')
            logfile.write('electrode ' + eoi + ' ... reference [%s]'  % ', '.join(map(str, ref)) + ' ... \n')
            
            sparsen_eeg_logfile(os, subject)
            #serial_to_parallel_codes(max_triggerlength, min_triggerlength, np, os, subject)
                 
            eoi_idx, raw =                          pp_raw(eoi, mne, os, ref, subject)

            if eoi in raw.info['bads']:
                break            
            #if ref in raw.info['bads']:
            #    continue
            
            raw =                                   pp_lpf(raw)
            eog_events, n_blinks, raw =             pp_blinks(mne, np, raw, subs)
            raw =                                   pp_reference(raw, ref)
            events, events_old, stim_events_single = pp_events(max_triggerlength, min_triggerlength, mne, np, raw)
            
            
            baseline, epochs, event_id, evoked, picks, raw, tmax, tmin = \
                                                pp_epochs_evoked(events, mne, os, plt, raw, subject)
            
            liste, liste_all =                  calculate_timings(math, np.mean, np, os, plt, stim_events_single, subject, subs)                                    
            


            
            
            #evoked_oz_cond = pp_evoked_oz_seperate(epochs, hemifields2, Ncond, np, os, picks, plt, screenfreq, stimfreq2, subject, time_per_stim_and_ISI2, tmin, tmax)
            #pp_evoked_oz_merged(division, evoked_oz_cond, math, np, os, plt, stimfreq, subject, tmin, tmax)
            
            #time_seperate_freq_merged_fields(cm, eoi, eoi_idx, epochs, np, os, picks, plt, ref, single_event_timings, stimdurations, stimfreq, subject, subs, tmin, tmax)
            #plt.close('all')
            time_seperate_freq_merged_fields_summation(cm, eoi, eoi_idx, epochs, np, os, picks, plt, ref, single_event_timings, stimdurations, stimfreq, subject, subs, tmin, tmax)
            plt.close('all')
            #epochs, montage, raw =              time_topomaps(epochs, hemifields2, mne, np, os, parent_dir, plt, raw, ref, stimdurations, stimfreq2, subject, subs)
            #plt.close('all')
            
            '''
            ###### PIPELINE - Frequency Analysis ############
            epochs, montage, raw =              freq_analysis(eoi, eoi_idx, epochs, hemifields2, mne, np, os, parent_dir, plt, raw, ref, stimfreq, stimfreq2, subject, subs)
            
            ###########
            
            sparsen_eeg_logfile(os, subject)
            eoi_idx, raw =                                                  pp_raw(eoi, mne, os, ref, subject)
            #raw =                                                           pp_lpf(raw)
            eog_events, n_blinks, raw =                                     pp_blinks(mne, np, raw, subs)
            events, events_old, stim_events_single  =                       pp_events(max_triggerlength, min_triggerlength, mne, np, raw)
            baseline, epochs, event_id, evoked, picks, raw, tmax, tmin =    pp_epochs_evoked(events, mne, os, plt, raw, subject)
            evoked_oz_cond =                                                pp_evoked_oz_seperate(eoi_idx, epochs, hemifields2, Ncond, np, os, picks, plt, screenfreq, stimfreq2, subject, time_per_stim_and_ISI2, tmin, tmax)
            pp_evoked_oz_merged(division, evoked_oz_cond, math, np, os, plt, stimfreq, subject, tmin, tmax)
          
            
'''            
###### Group analysis

#for eoi in ['Oz','O1', 'O2', 'Pz']:
for eoi in ['Pz']: # 'Oz'
    
    for ref in [['TP9','TP10']]:
    #for ref in [['TP9','TP10']]:    # ['Cz'],['Fp1','Fp2']
        # decide, if you want to plot the summation of UVF/LVF, too
        summation = 1 # 1 = true, 0 = false
        Nsub = 7 # additionally change the subjects in the loop of the subsequent function
        num_period = 3 # single integer, start counting at zero
        # for displaying the whole time snipped
        GROUP_TIME_seperate_freq_merged_fields(cm, eoi, np, Nsub, os, plt, ref, single_event_timings, stimdurations, stimfreq, summation, tmin, tmax)
        
        # if you want only 1 period displayed
        #GROUP_TIME_seperate_freq_merged_fields_zoom_period(cm, eoi, np, Nsub, num_period, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax)
        
        # if you want an average of all periods per condition
        #GROUP_TIME_seperate_freq_merged_fields_zoom_period_all(cm, eoi, np, Nsub, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax)
        GROUP_TIME_seperate_freq_merged_fields_zoom_period_all_sameX(cm, eoi, hemifields, np, Nsub, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax)
        # plot the difference wave (DIFF) in the same fashion
        mse = False
        absolute = True
        GROUP_TIME_seperate_freq_merged_fields_zoom_period_all_sameX_DIFF(absolute,cm, eoi, np, mse, Nsub, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax)
        # do the same for all single subjects
        # dies auch noch machen
        SINGLE_SUBJECTS_TIME_seperate_freq_merged_fields_zoom_period_all_sameX_DIFF(cm, eoi, np, Nsub, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax)

        bin_length=10 # the size of the bins
        #group_additivity_over_time(bin_length, eoi,hemifields, mne, np, os, plt, ref, stimfreq)
        group_similarity_over_time(bin_length, eoi,hemifields, mne, np, os, plt, ref, stimfreq)


'''
from rk_chf_processing_functions import *
Nsub = 7         
summation = 1
mse = False
absolute = True
testarr = GROUP_TIME_seperate_freq_merged_fields_zoom_period_all_sameX_DIFF(absolute,cm, eoi, np, mse, Nsub, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax)



#########
array1 = np.asarray([1,2,3,4,5,6])
array2 = np.asarray([2,3,4,4,3,3])
mse1   = (array1 - array2) ** 2
'''