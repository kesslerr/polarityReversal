#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 17:57:50 2016
behavioral logfiles analysis of the CHF paradigm
@author: cth
"""

from __future__ import division
import os
from pyface.qt import QtGui, QtCore
import numpy as np
import matplotlib.pyplot as plt
import math
import shutil
from itertools import repeat

os.chdir('/media/cth/Samsung USB/CHF/processing_new')

for sub in range(4,4): ## later 11
    if sub<10:
        subs = '0' + str(sub)
    else:
        subs = str(sub)
    subject='RK_200' + subs
    
    os.chdir(subject + os.sep + 'log')
    event = open('log_event.txt', 'r').readlines()
    stimfreq = open('log_stimfreq.txt', 'r').readlines()
    hemifield = open('log_hemifield.txt', 'r') .readlines()   
    time_exp = open('log_time_exp.txt', 'r').readlines()
    time_block = open('log_time_block.txt', 'r').readlines()
    act_color = open('log_act_color.txt', 'r').readlines()
    
    # convert the variables to int/floats and delete \n in strings, so that you can work with them
    for i in range(0,len(event)):
        event[i]=event[i].strip()
        stimfreq[i]=float(stimfreq[i].strip())
        hemifield[i]=hemifield[i].strip()
        time_exp[i]=float(time_exp[i].strip())
        time_block[i]=float(time_block[i].strip())
        act_color[i]=int(act_color[i].strip())
        
    
    resp_cor = 0
    resp_incor = 0
    rts = []
    actcolor = 0
    color_onset = []
    stim_onset_1 = []
    block_onset_tmp = 0
    block_on = False
    # iterate through lines of the files
    for i in range(0,len(event)):
        
        # responses
        if 'keypress_1_incorrect' in event[i]:
            resp_incor += 1
        elif 'keypress_1_correct' in event[i]:
            resp_cor += 1
            rts.append(float(time_exp[i]) - color_onset)
        elif 'keypress_2_incorrect' in event[i]:
            resp_incor += 1
        elif 'keypress_2_correct' in event[i]:
            resp_cor += 1
            rts.append(float(time_exp[i]) - color_onset)
        # response times
        elif 'colorchange_color_0'  in event[i]:
            actcolor=0
        elif 'colorchange_color_1'  in event[i]:
            actcolor=1
            color_onset = float(time_exp[i])
        elif 'colorchange_color_2'  in event[i]:
            actcolor=2
            color_onset = float(time_exp[i])
        
        # stimulus timings
        elif 'onset_stim'  in event[i]:
            block_act_freq = float(stimfreq[i])
            block_act_VF = hemifield[i]
            onset_block = float(time_exp[i])
        elif 'IBI_begin'
            block_on = False
            

    hitrate = resp_cor / (resp_cor + resp_incor)
    