#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 17:29:29 2017

single subject analysis of CHF paradigm

topics: fourier transform, peak identification, ...

@author: kesslerr@med.uni-marburg.de
"""


# coding: utf-8

# # CHF - Statistics
# 
# This script is the follow-up of "rk_find_peaks4" for the CHF data.
# 
# author: *roman kessler* (kesslerr@med.uni-marburg.de)
# 
# date: *Sep 24, 2017* (election day)
# 
# What it does:
# 1. it reads the evoked responses from a text file
# 2. finding the maximum/minimum average values in each condition
# 3. doing a fourier transformation of the files to extract
#     a. amplitude information (not yet corretly converted)
#     b. phase shift
# 4. plots the results
# 

# In[6]:

#get_ipython().magic(u'matplotlib inline')
#get_ipython().magic(u'pylab inline --no-import-all')
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import plotly.plotly as py
import csv
import math
from pylab import *
import os 


# # group responses

# # for each subject

for subjectno in range(4,11): # load subjects # 4 to subject # 10
    if subjectno < 10:
        subjectstr = 'RK_2000' + str(subjectno)
    else:
        subjectstr = 'RK_200'  + str(subjectno)


pathExp    = '/media/cth/Samsung_T3/CHF/processing_new/'
pathBase   = '/media/cth/Samsung_T3/CHF/processing_new/group/data/'
pathFiles1 = 'evoked_group_'
pathFiles1z = 'evoked_zoom_avg_'
pathCond   = ['cond_1','cond_2','cond_3','sum_1','cond_4','cond_5','cond_6','sum_2','cond_7','cond_8','cond_9','sum_3','cond_10','cond_11','cond_12','sum_4']
pathFiles2 = '_eoi_Oz_ref_[TP9, TP10].txt'

dataMin  = [];     dataMax  = [];      dataAbs   = []
dataMinz = [];     dataMaxz = [];      dataAbsz  = []

evoked_data = []
for icond in pathCond:  
    
    data      = np.loadtxt(pathBase + pathFiles1 + str(icond) + pathFiles2)                 # load ERPs
    dataMin.append(np.amin(data[322:-1]))                                                   # find maximum value
    dataMax.append(np.amax(data[322:-1]))                                                   # find minimum value
    dataAbs.append(np.amax([abs(dataMax[-1]),abs(dataMin[-1])]))                      # absolute value
    evoked_data.append(data)
    
    data      = np.loadtxt(pathBase + pathFiles1z + str(icond) + pathFiles2)                # load zoomed ERPs
    dataMinz.append(np.amin(data))                                                          # find maximum value
    dataMaxz.append(np.amax(data))                                                          # find minimum value    
    dataAbsz.append(np.amax([abs(dataMaxz[-1]),abs(dataMinz[-1])]))                   # absolute value
    


# # phase of the signal

# ## fourier transformation

# In[8]:

periodDur            = [546,279,145,78]
Nstim                = [4,8,15,30]
stimFreq             = [2,4,8,15]
visualFields         = ['FVF','UVF','LVF','SUM']
periodDur            = [periodDur[i//4] for i in range(16)]
Nstim                = [Nstim[i//4] for i in range(16)]

def fourier(coi, evoked_data, plot, save):
    y = evoked_data[coi]                                        # change here to get a full period length 
    y = y[ 322 : 322 + Nstim[coi] * periodDur[coi] ]            # eliminate leakage effect, ensure periodicity of the signal

    Fs = len(y);                                                # sampling rate
    Ts = 1;                                                     # sampling interval
    t = np.arange(0,Fs,1)                                       # time vector

    n = len(y)                                                  # length of the signal
    k = np.arange(n)
    T = n/Fs
    frq = k/2 #T                                                # two sides frequency range
    frq = frq[range(int(n/2))]                                  # one side frequency range

    Y = np.fft.fft(y)/n                                         # fft computing and normalization
    Y = Y[range(int(n/2))]

    border=40
    
    if plot == True:
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(t,y)
        ax[0].set_xlabel('Time')
        ax[0].set_ylabel('Amplitude')
        ax[1].plot(frq,abs(Y),'r')                              # plotting the spectrum
        xlim(0,border)
        ax[1].set_xlabel('Freq (Hz)')
        ax[1].set_ylabel('|Y(freq)|')
    
    if save == True:
        plt.savefig(pathExp + 'group' + os.sep + 'plots' + os.sep + 'fourier_1000dpi_' + str(stimFreq[int(coi/4)]) + visualFields[int(coi%4)], dpi=1000)
        plt.savefig(pathExp + 'group' + os.sep + 'plots' + os.sep + 'fourier_50dpi_'   + str(stimFreq[int(coi/4)]) + visualFields[int(coi%4)], dpi=50)
    
    angle = np.angle(Y)
    angledeg = [math.degrees(i) for i in angle]; #angledeg

    idx, = np.where( abs(Y) == amax(abs(Y)) )
    #print('index of greatest value: ' + str(idx[0]) )
    frequency = frq[idx[0]]
    #print('frequency with greatest power: ' + str(frequency) )
    shiftDeg = angledeg[idx[0]]
    #print('angle of greatest value: ' + str(shift))
    shiftRad = angle[idx[0]]
    
    return shiftDeg, shiftRad, frequency, y


# In[9]:

shiftDeg               = []
shiftRad               = []
frequency              = []
evoked_data_sparse     = []                            # the full-period data, without start and end-time

plot = False                                           # plot the results?
save = False                                           # save the figures?

for i in range(len(evoked_data)):                      # loop through conditions
    a, b, c, y = fourier(i, evoked_data, plot, save)
    shiftDeg.append(a)
    shiftRad.append(b)
    frequency.append(c)
    evoked_data_sparse.append(y)


# In[10]:

len(evoked_data)


# ### normalizing angles
# Use the first of 4 values (full visual field stimulation) as 0 degree shift. Proceed with the other angles in relation to this value.

# In[11]:

shiftRadNorm = []                      # write normalized angles (norm. to FVF stimulation) in a new list
shiftDegNorm = []
for icond in range(16):                # 16 conditions
    shiftRadNorm.append(shiftRad[icond] - shiftRad[int(floor(icond/4)*4)])
    shiftDegNorm.append(shiftDeg[icond] - shiftDeg[int(floor(icond/4)*4)])
#    print(shiftRad[int(floor(icond/4)*4)])


# ### converting angles to millisecons delay
# 
# lag in ms = radial shift  x  length of full cycle in ms / ( 2 * pi )

# In[12]:

shiftTimeNorm = []                                                 # write ms delay in a new list
for icond in range(16):                                            # 16 conditions
    timeCycle = periodDur[icond]
    # calculate time shift
    timeShift = timeCycle * shiftRadNorm[icond] / ( 2 * np.pi )    # formula to calculate times
    shiftTimeNorm.append(timeShift)


# This is the average delay (in mm) of the respective Hz response, normalized to the FVF response (which is set to 0 mm)

# In[13]:

#shiftTimeNorm


# This is the average delay in degree. Normalized to FVF response (which is set to 0 degree). Note, that the angle is with respect to the frequency, thus a particular angle in the 2 Hz conditions means a higher time difference than in the 4 Hz condition (and so on).

# In[14]:

#shiftDegNorm


# ### polar plot of responses
# example fot this polar bar plot: https://matplotlib.org/examples/pie_and_polar_charts/polar_bar_demo.html
# or a polar scatter plot: https://matplotlib.org/examples/pie_and_polar_charts/polar_scatter_demo.html

# In[129]:

for i in range(4):                                          # compute for each frequency separately
    N      = 4                                              # N of conditions per frequency: here: FVF, UVF, LVF, SUM
    theta  = shiftRadNorm[i*4:i*4+4]                        # position on the circle (like position on the x-axis)


    radii  = 1000000 * np.array(dataAbsz[i*4:i*4+4])        # heigh of the bars
    width  = np.pi / 25                                     # width of the bars
    colors = [(0,0,0),(51/255,160/255,44/255),(31/255,120/255,180/255),(0.5,0.5,0.5)]    # here black, green, blue, gray

    ax = plt.subplot(111, projection='polar')
    ax.set_title(str(stimFreq[i]) + ' Hz stimulation\n')
    bars = ax.bar(theta, radii, width=width, bottom=0.0)
    
    # Use custom colors and opacity
    for r, bar, color in zip(radii, bars, colors):
        bar.set_facecolor(color)
        bar.set_alpha(0.7)
    
    ax.legend(bars, visualFields, loc='best')
    
    # save figures
    plt.savefig(pathExp + 'group' + os.sep + 'plots' + os.sep + 'polar_1000dpi_' + str(stimFreq[i]), dpi=1000)
    plt.savefig(pathExp + 'group' + os.sep + 'plots' + os.sep + 'polar_50dpi_' + str(stimFreq[i]), dpi=50)
    
    # show plot
    plt.show()


# ## cross-correlation as a measure of phase

# to do ? or is a FFT sufficient ?
# 

# # similarity between curves

# ## sum of squares (here: mean squared error)
# For each x value, compute delta_yi = y1,i - y2,i and accumulate delta_yi2. This metric is the basis for a least square optimization, where the goal is to minimize the sum of the squares of the errors. This is a widely used approach because oftentimes it is fairly easy to implement.
# https://stackoverflow.com/questions/6723157/how-to-compare-two-curves-arrays-of-points
# 
# MSE = 1 / n * sum ( ( Y(pred) - Y(obs) ) ^ 2 )

# In[63]:

Ncond = len(evoked_data_sparse)                              # total # conditions
MSE = zeros((4,4))                                # empty array to store all comparisons
for i1, i1i in zip( [0,4,8,12] , range(4) ):           # go in 1st dimension
    for i2, i2i in zip( range(i1,i1+4), range(4) ):  # go in 2nd dimension
        SE = []
        Ntp   = len(evoked_data_sparse[i1])                  # total # time points
        
        for j in range(Ntp):                                 # go through the data
            Yp = evoked_data_sparse[i1][j]
            Yo = evoked_data_sparse[i2][j]
            SE.append( ( ( Yp - Yo ) ** 2 ) )
            
        MSE[i1i][i2i] = 1 / Ntp * sum(SE)

fig = plt.imshow(MSE, cmap='gray', interpolation='nearest')  # plot
plt.colorbar()
xlabels = ['FVF','UVF','LVF','SUM']
ylabels = ['2 Hz','4 Hz','8 Hz','15 Hz']
plt.xticks(range(len(xlabels)), xlabels, fontsize=12)
plt.yticks(range(len(ylabels)), ylabels, fontsize=12)
plt.savefig(pathExp + 'group' + os.sep + 'plots' + os.sep + 'MSE_50dpi',   dpi=50)
plt.savefig(pathExp + 'group' + os.sep + 'plots' + os.sep + 'MSE_1000dpi', dpi=1000)
plt.show()

print(MSE)

# only for the SUM condition compared to FVF
sumMSE = MSE[:][-1] # In your specific case you got a 1D array, so you need to add a dimension with np.expand_dims()
sumMSE = np.asarray(sumMSE)
sumMSE = np.expand_dims(sumMSE, axis=0)
sumMSE = np.transpose(sumMSE)

fig = plt.imshow(sumMSE, cmap='gray', interpolation='nearest')  # plot
plt.colorbar()
xlabels = ['SUM']
ylabels = ['2 Hz','4 Hz','8 Hz','15 Hz']
plt.xticks(range(len(xlabels)), xlabels, fontsize=12)
plt.yticks(range(len(ylabels)), ylabels, fontsize=12)
plt.savefig(pathExp + 'group' + os.sep + 'plots' + os.sep + 'MSEsum_50dpi',   dpi=50)
plt.savefig(pathExp + 'group' + os.sep + 'plots' + os.sep + 'MSEsum_1000dpi', dpi=1000)
plt.show()


# In[59]:

get_ipython().magic(u'pinfo np.asarray')


# ## maximum deviation
# Find the abs_delta_yi = |y1,i - y2,i| that maximizes the |y1,i - y2,i| for all x values. This metric is the basis for a lot of the implementations of the functions in the math library, where the goal is to minimize the maximum error. These math library implementations are approximations of the true function. As a consumer of such an approximation, I typically care more about the worst thing that the approximation is going to do to my application than I care about how that approximation is going to behave on average. https://stackoverflow.com/questions/6723157/how-to-compare-two-curves-arrays-of-points

# ## weighted least squares
# 'In my research I compared two curves (data vs model) and I used the  weighted least squares (WLS)
# For the WLS method it was defined as the sum of the absolute differences between observed and expected values, divided by the observed values (Lika et al. 2011).'
# https://www.researchgate.net/post/Is_there_any_statistical_method_to_compare_two_curves_in_a_graph_eg_quantifying_how_similar_different_they_are
# 
# WLS = 1 / n * sum ( | 1 - ( expected value / observed value) | )
Ncond = len(evoked_data)                                     # total # conditions
Ntp   = len(evoked_data[0])                                  # total # time points
WLS = zeros( (Ncond, Ncond) )                             # empty array to store all comparisons
for i1 in [0,4,8,12]:                                      # go in 1st dimension
    for i2 in range(i1+1,i1+4):                                  # go in 2nd dimension
        LS=[]
        for j in range(Ntp):                                 # go through the data
            ev = evoked_data[i1][j]   # put here FVF ?
            ov = evoked_data[i2][j]   # put here the other VFs ?
            LS.append( abs( 1 - ( ev / ov ) ) )
        WLS[i1][i2] = 1 / Ntp * sum(LS)

plt.imshow(WLS, cmap='gray', interpolation='nearest')
plt.colorbar()
plt.show()
# Ncond = len(evoked_data)                                     # total # conditions
# Ntp   = len(evoked_data[0])                                  # total # time points
# WLS2 = zeros( (Ncond, Ncond) )                             # empty array to store all comparisons
# for i1 in range(Ncond):                                      # go in 1st dimension
#     for i2 in range(Ncond):                                  # go in 2nd dimension
#         LS=[]
#         for j in range(Ntp):                                 # go through the data
#             ev = evoked_data[i1][j]
#             ov = evoked_data[i2][j]
#             LS.append( abs( ( 1 - ev / ov ) ) )
#         WLS2[i1][i2] = 1 / Ntp * sum(LS)
# fig = plt.imshow(WLS2, cmap='hot', interpolation='nearest')
# plt.colorbar()
# plt.show()