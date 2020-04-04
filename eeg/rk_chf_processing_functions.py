# -*- coding: utf-8 -*-
"""
Processing of the rk_CHF paradigm.
Under construction.
"""


def sparsen_eeg_logfile(os, subject):
    if not os.path.isfile(subject + os.sep + 'eeg' + os.sep + 'old_' + subject + '.vmrk'):
        os.chdir(subject + os.sep + 'eeg')
        file = open(subject + '.vmrk', 'r') #'w'
        new_file = open('new_' + subject + '.vmrk', 'w')
        for line in file:
            if "S  3" not in line:  # substring
                new_file.write(line)
                continue    
        file.close()
        new_file.close()
        # rename
        os.rename(subject + '.vmrk','old_' + subject + '.vmrk')
        os.rename('new_' + subject + '.vmrk' , subject + '.vmrk')
        os.chdir('../..')
    return

# necessary for parallel logfile
def serial_to_parallel_codes(max_triggerlength, min_triggerlength, np, os, subject):
    os.chdir(subject + os.sep + 'eeg')    
    file = open(subject + '.vmrk', 'r')
    new_file = open('tmp_' + subject + '.vmrk', 'w')
    for line in file:
        if "S  1" in line or "S  2" in line:  # only choose relevant lines
            line=line.replace('S  ', '') # substitute "S__x" string with simple number, for easier handling in the future
            new_file.write(line) # write to new file
            continue    
    file.close()
    new_file.close()
    
    data = np.genfromtxt('tmp_' + subject + '.vmrk',delimiter=',') # import newly created data table
    onoff=data[:,1]
    onofftime=data[:,2]
    
    onofftime=[i*2 for i in onofftime] # convert onsets/offsets to real ms timing, from 500 Hz to 1000 Hz
    
    # loop throuth the single events/lines
    onset_tmp0=0.0
#    onset_tmp1=0.0
    event_type=[]
    event_time=[]
    trashtriggercounter=0
    for i in range(0,len(onoff)):
        if onoff[i]==1:
            onset_tmp0=onofftime[i]
            #onset_tmp1=onofftime[i]
        elif onoff[i]==2:
            dist=onofftime[i]-onset_tmp0
            if dist < (min_triggerlength[0]*1000):
                # hier etwas mit den einzelnen events machen spaeter?
                trashtriggercounter += 1
            else:
                for j in range(0,len(min_triggerlength)):
                    if (min_triggerlength[j]*1000) < dist < (max_triggerlength[j]*1000):
                        event_type.append(j+1)
                        event_time.append(onset_tmp0)
    
    event_time = [i/2 for i in event_time] # back from ms to sampling

    file = open(subject + '.vmrk', 'r')
    file_old = open('old_' + subject + '.vmrk', 'r')
    new_file = open('par_' + subject + '.vmrk', 'w')
    counter = 0
    for line_old in file_old: # copy the first header lines from the very old file to the new file
        new_file.write(line_old) # write to new file
        counter += 1
        if counter > 12:
            break
        
    for line in file:
        if "S  1" in line:  # only choose relevant lines
            for i in range(0,len(event_time)):
                if ',' + str(int(event_time[i])) + ',' in line:
                    line=line.replace('S  1', 'S  ' + str(event_type[i+1]))# replace the S**1 with the new S**1-12
                    new_file.write(line) # write to new file
    
    file.close()
    new_file.close()
                        
    os.rename(subject + '.vmrk','old2_' + subject + '.vmrk')
    os.rename('par_' + subject + '.vmrk' , subject + '.vmrk')
    os.chdir('../..')
    return
    
    

####################################################################
############ PREPROCESSING #########################################
####################################################################

def pp_raw(eoi, mne, os, ref, subject):
    # read raw data
    raw=mne.io.read_raw_brainvision(subject + os.sep + 'eeg' + os.sep + subject + '.vhdr',
                            #event_id={'R  64':64;'R 128':128},
                            preload=True
                            )

    eoi_idx=raw.ch_names.index(eoi)
    # if necessary, exclude some channels because of bad data quality
    if os.path.isfile(subject + os.sep + 'eeg' + os.sep + 'bad_channels.txt'):
        file = open(subject + os.sep + 'eeg' + os.sep + 'bad_channels.txt', 'r') #'w'
        for line in file:
            raw.info['bads'].append(line)
#            content = [x.strip('\n') for x in content] 
            raw.info['bads'] = [x.strip('\n') for x in raw.info['bads']] # delete the \n for every element
        file.close()    
    else:
        print("no bad channels to be rejected")
    return eoi_idx, raw
    
def pp_lpf(raw): # low-pass filter
    #raw_unf=raw
    raw.filter(None, 47.5, h_trans_bandwidth='auto', filter_length='auto',
               phase='zero')
    return raw#, raw_unf

def pp_blinks(mne, np, raw, subs): # eye blink detection
    # Fp1
    eog_events1 = mne.preprocessing.find_eog_events(raw, ch_name='Fp1')
    eog_events2 = mne.preprocessing.find_eog_events(raw, ch_name='Fp2')
    eog_events = np.concatenate((eog_events1, eog_events2))
    
    n_blinks = len(eog_events)
    # Center to cover the whole blink with full duration of 0.5s:
    onsets = eog_events[:, 0] / raw.info['sfreq'] - 0.25
    durations = np.repeat(0.5, n_blinks)
    descriptions = ['bad blink'] * n_blinks
    
    raw.annotations = mne.Annotations(onsets, durations, descriptions,
                                      orig_time=raw.info['meas_date'])
    raw.plot(events=eog_events)  # To see the annotated segments.

    # save number of detected blinks per subject
    logfile = open('n_blinks.txt', 'a')
    logfile.write('subject ' + subs + ' has ' + str(n_blinks) + ' blinks \n')

    return eog_events, n_blinks, raw

def pp_reference(raw, ref):
    # setting EEG references
    #mne.set_eeg_reference(raw, ref_channels = ref, copy = False)
    if ref:
        raw.set_eeg_reference(ref_channels=ref)
    print('channels')
    print(raw.ch_names)
    print('reference')
    print(ref)
    return raw
    
#############################
### TIME DOMAIN ANALYSIS ####
#############################    
    
def pp_events(max_triggerlength, min_triggerlength, mne, np, raw): # get serial port signal
    events = mne.find_events(raw, stim_channel='STI 014')
    print('events')
    print(events)
    # create here appropriate trigger values for the conditions
    eventsman=[]
    trashtriggercounter=0
#    stim_events_all=[]
    stim_events_single=[]
    last_cond=0
    onset_cond_tmp=0.0
    
    onset_tmp0=0.0
    onset_tmp1=0.0
    
    for i in range(0,len(events)):
        if events[i][2]==64:
            eventsman.append([events[i][0],events[i][1],101])
        elif events[i][2]==128:
            eventsman.append([events[i][0],events[i][1],102])
        elif events[i][2]==1:
            onset_tmp0=events[i][0]
            onset_tmp1=events[i][1]
            
        elif events[i][2]==2:
            dist=(events[i][0]-onset_tmp0)/500.0 # wenn in 500 Hz geloggt
            if dist<min_triggerlength[0]:
                stim_events_single.append([last_cond,(events[i][0]-onset_cond_tmp)*2]) # somit gleich in [ms]
                trashtriggercounter+=1
                
            else:
                for j in range(0,len(min_triggerlength)):
                    if min_triggerlength[j] < dist < max_triggerlength[j]:
                        eventsman.append([onset_tmp0,onset_tmp1,j+1])
                        last_cond = j+1
                        onset_cond_tmp=onset_tmp0
                    
    eventsman2=np.ndarray( (len(eventsman),len(eventsman[0])) ,dtype = 'int')
    for i in range(0,len(eventsman)):
        for j in range(0,len(eventsman[0])):
            eventsman2[i][j]=eventsman[i][j]
    events_old=events
    events=eventsman2
    
    return events, events_old, stim_events_single

def pp_epochs_evoked(events, mne, os, plt, raw, subject):
    # weise den IDs ein Event zu
    event_id = dict(cond_1=1, cond_2=2, cond_3=3, cond_4=4, cond_5=5, cond_6=6, 
                    cond_7=7, cond_8=8, cond_9=9, cond_10=10, cond_11=11, cond_12=12)
                    #resp_1=101, resp_2=102)  # event trigger and conditions
    tmin = 0.0  # start of each epoch (0ms before the trigger)
    tmax = 3.0  # end of each epoch (2500ms after the trigger)
    picks = mne.pick_types(raw.info, meg=False, eeg=True,
                           #eog=True,
                           stim=False, exclude='bads'
                           )
    baseline = (0, 0.322)  # means from the first instant to t = 0
# original    baseline = (None, None)  # means from the first instant to t = 0
    reject = dict(#grad=4000e-13, # T / m (gradiometers)
                  #mag=4e-12, # T (magnetometers)
                  eeg=40e-6, # V (EEG channels) ####### war 40, dies ist eher empirisch ausgewaehlt
    #              eeg=40e-6, # V (EEG channels)
                  #eog=250e-6 # V (EOG channels)
                  )
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, proj=True, picks=picks,
                        baseline=baseline, preload=False,
    #                    reject=reject,
                        reject_by_annotation = True,
                        #add_eeg_ref=False,
                        detrend=0) # detrend = 0 bedeuted Mittelwert auf 0, detrend = 1 waere lineares detrenden 
    
    epochs.drop_bad()
    # plot drop log... optional
    print(epochs.drop_log) #[40:45])  # only a subset
    epochs.plot_drop_log(show=False) 
    plt.savefig(subject + os.sep + 'plots' + os.sep + 'drop_log.png', dpi=300)#1000)
    plt.close('all')
    
    print(epochs)
    evoked = epochs.average()
    
    '''     new: save the epochs and evoked in mne fashion, so you can use them in future analyses '''
    if not os.path.exists(subject + os.sep + 'data_mne'):
        os.makedirs(subject + os.sep + 'data_mne')
    epochs.save(subject + os.sep + 'data_mne' + os.sep + 'epochs.fif')
    evoked.save(subject + os.sep + 'data_mne' + os.sep + 'evoked.fif')  # save evoked data to disk
    
    
    # plot condition and response distribution over time
    mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp,
                        event_id=event_id, show=False)
    plt.savefig(subject + os.sep + 'plots' + os.sep + 'condition_order.png', dpi=300)#1000)
    plt.close('all')
    return baseline, epochs, event_id, evoked, picks, raw, tmax, tmin

def time_seperate_freq_merged_fields(cm, eoi, eoi_idx, epochs, np, os, picks, plt, ref, single_event_timings, stimdurations, stimfreq, subject, subs, tmin, tmax):
    # plot evoked of each frequency combined
    for i in range(1,5):
        plt.figure()
        colors = [ cm.Paired(x) for x in np.linspace(0.1, 0.4, 3) ]
        yrange = []
        for j in range(1,4):
            # resampling = copying every data point
            evoked=epochs['cond_' + str((i-1)*3 + j)].average(picks=[eoi_idx]).data[0]
            evoked2=[]
            for k in range(0,len(evoked)*2):
                evoked2.append(evoked[np.floor(k/2)])
            yrange.append(np.ceil(max(max(evoked2),abs(min(evoked2)))*100000))
            plt.plot([ x * 100000 for x in evoked2 ],dth=3, color = colors[j-1])
            #save evoked2 to file for group analyses
            np.savetxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt',
                       evoked2 ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='')
            #np.savetxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + '.txt',
            #           evoked2 ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='')

            
        plt.legend(['FVF','LVF','UVF'])
        plt.title('subject ' + subs + ' / ' + str(stimfreq[i-1]) + ' Hz' + ' / electrode ' + eoi + ' / reference [%s]' % ', '.join(map(str, ref)))
        plt.ylabel('$\mu$V')
        plt.xlabel('ms')
        
        yrange2=max(yrange)
        plt.ylim(-yrange2,yrange2)
        plt.xlim(tmin*1000,tmax*1000)
        plt.hlines(y=0, xmin=tmin*1000, xmax=tmax*1000)
        # plot stimulus onset vertical lines
        for k in range(0,len(single_event_timings[i-1])):
            plt.vlines(single_event_timings[i-1][k], -yrange2, yrange2, color='red', linestyles='solid', label='onset')
        #k=0
        #while (0.334 + k * stimdurations[i-1]) < 3.0:
        #    plt.vlines(334 + k * stimdurations[i-1]*1000, -yrange2, yrange2, color='red', linestyles='solid', label='onset')
        #    k += 1
        # plot vertical lines for all 100ms
        k=0
        while k <= 3000:
            plt.vlines(k, -yrange2, yrange2, color='gray', linestyles='dashed', label='')
            k += 100
        
        #plt.savefig(subject + os.sep + 'plots' + os.sep + 'evoked_' + str(stimfreq[i-1]) + '.png', dpi=300)#1000)
        plt.savefig(subject + os.sep + 'plots' + os.sep + 'evoked_subject_' + subs + '_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)))
        plt.close('all')
        #plt.show()
    return # evoked_all

def time_seperate_freq_merged_fields_summation(cm, eoi, eoi_idx, epochs, np, os, picks, plt, ref, single_event_timings, stimdurations, stimfreq, subject, subs, tmin, tmax):
    # plot evoked of each frequency combined
    for i in range(1,5):
        plt.figure()
        colors = cm
        #old colors = [ cm.Paired(x) for x in np.linspace(0.1, 0.4, 3) ]
        #old version: colors = [ cm.Paired(x) for x in np.linspace(0.1, 0.9, 4) ]
        yrange = []
        evoked_summation = []
        for j in range(1,5):
            if j < 4:
                # resampling = copying every data point
                evoked=epochs['cond_' + str((i-1)*3 + j)].average(picks=[eoi_idx]).data[0]
                evoked2=[]
                for k in range(0,len(evoked)*2):
                    evoked2.append(evoked[int(np.floor(k/2))])
                yrange.append(np.ceil(max(max(evoked2),abs(min(evoked2)))*100000))
                plt.plot([ x * 100000 for x in evoked2 ],linewidth=2, color = colors[j-1])
                #save evoked2 to file for group analyses
                np.savetxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt',
                           evoked2 ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='')
                #np.savetxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + '.txt',
                #           evoked2 ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='')
                
                # add values of upper and lower hemifield for a summation effect
                if j == 2: # if upper hemifield, save values for next iteration of loop to sum with lower hemifield          
                    evoked_summation = evoked2
                elif j == 3:
                    for k in range(0,len(evoked2)):
                        evoked_summation[k] = evoked_summation[k] + evoked2[k]
            elif j == 4:
                yrange.append(np.ceil(max(max(evoked_summation),abs(min(evoked_summation)))*100000))
                plt.plot([ x * 100000 for x in evoked_summation ],linewidth=2,linestyle='--', color='k')
                np.savetxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + 'summation__eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt',
                               evoked_summation ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='')
                            
        plt.legend(['FVF','LVF','UVF','SUM'])
        plt.title('subject ' + subs + ' / ' + str(stimfreq[i-1]) + ' Hz' + ' / electrode ' + eoi + ' / reference [%s]' % ', '.join(map(str, ref)))
        plt.ylabel('$\mu$V')
        plt.xlabel('ms')
        
        yrange2=max(yrange)
        plt.ylim(-yrange2,yrange2)
        plt.xlim(tmin*1000,tmax*1000)
        plt.hlines(y=0, xmin=tmin*1000, xmax=tmax*1000)
        # plot stimulus onset vertical lines
        for k in range(0,len(single_event_timings[i-1])):
            plt.vlines(single_event_timings[i-1][k], -yrange2, yrange2, color='red', linestyles='solid', label='onset')
        
        #k=0
        #while (0.334 + k * stimdurations[i-1]) < 3.0:
        #    plt.vlines(334 + k * stimdurations[i-1]*1000, -yrange2, yrange2, color='red', linestyles='solid', label='onset')
        #    k += 1
        # plot vertical lines for all 100ms
        k=0
        while k <= 3000:
            plt.vlines(k, -yrange2, yrange2, color='gray', linestyles='dashed', label='')
            k += 100
        
        #plt.savefig(subject + os.sep + 'plots' + os.sep + 'evoked_' + str(stimfreq[i-1]) + '.png', dpi=300)#1000)
        plt.savefig(subject + os.sep + 'plots' + os.sep + 'evoked_subject_' + subs + '_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)))
        plt.close('all')
        #plt.show()
    return # evoked_all

def time_topomaps(epochs, hemifields2, mne, np, os, parent_dir, plt, raw, ref, stimdurations, stimfreq2, subject, subs):
    montage = mne.channels.read_montage(parent_dir + '/electrode_positions/M1_ThetaPhi_32_electrodes.txt', ch_names=None, path=None, unit='m', transform=False)    
    raw.set_montage(montage)
    epochs.set_montage(montage)
    for i in range(0,12):
        # coarser temporal resolution
        times=[0.322+0.1*i for i in np.linspace(0,1,num=19)]
        evoked_tmp = epochs['cond_'+str(i+1)].average()
        evoked_tmp.plot_topomap(times = times)
        plt.savefig(subject + os.sep + 'plots' + os.sep + 'topomap_time_coarse_subject_' + subs + '_' + str(int(np.round(stimfreq2[i]))) + '_Hz_' + hemifields2[i] + '_visual_field__reference [%s]' % ', '.join(map(str, ref)))
        plt.close('all')
        # finer temporal resolution
        times=[0.332+0.1*i for i in np.linspace(0,5,num=19)]
        evoked_tmp.plot_topomap(times = times)
        plt.savefig(subject + os.sep + 'plots' + os.sep + 'topomap_time_fine_subject_' + subs + '_' + str(int(np.round(stimfreq2[i]))) + '_Hz_' + hemifields2[i] + '_visual_field__reference [%s]' % ', '.join(map(str, ref)))
        plt.close('all')
    return epochs, montage, raw

################################
##### GROUP ANALYSES ###########
################################

def GROUP_TIME_seperate_freq_merged_fields(cm, eoi, np, Nsub, os, plt, ref, single_event_timings, stimdurations, stimfreq, summation, tmin, tmax):
    # plot evoked of each frequency combined FOR ALL SUBJECTS AVERAGED
    #colors = [ cm.Paired(x) for x in np.linspace(0.1, 0.4, 3) ]
    #colors.append((0.5, 0.5, 0.5, 1.0)) # gray for the summation of lower und upper hemifield, if summation is selected
    #Nsub = 10
    colors = cm
    
    jrange=4+summation
    for i in range(1,5): # Frequencies
        plt.figure()
        yrange=[]

        for j in range(1,jrange): # hemifields + summation of upper and lower hemifield
        
            evoked_array = np.ndarray((0,0))
            for sub in range(4,11): # subjects ######## CHANGE HERE
                if sub<10:
                    subs = '0' + str(sub)
                else:
                    subs = str(sub)
                subject='RK_200' + subs
                if j<4: #for the full, upper and lower hemifield, not for the summation
                    tmp = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ') 
                else:
                    tmp = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + 'summation__eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ') 
                evoked_array = np.append(evoked_array, tmp)    
            
            evoked_array = evoked_array.reshape((Nsub,3002))
            evoked_group = np.average(evoked_array, axis=0)
                
        
            yrange.append(np.ceil(max(max(evoked_group),abs(min(evoked_group)))*100000))
            plt.plot([ x * 100000 for x in evoked_group ],linewidth=2, color = colors[j-1])
            
            # save data
            # for normal data
            if j < 4:
                np.savetxt('group' + os.sep + 'data' + os.sep + 'evoked_group_cond_' + str((i-1)*3 + j) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt',
                           evoked_group ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='') 
            # for summation data
            else:
                np.savetxt('group' + os.sep + 'data' + os.sep + 'evoked_group_sum_' + str(i) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt',
                           evoked_group ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='') 
        
   
        if summation:
            plt.legend(['FVF','LVF','UVF','SUM'])
        else:
            plt.legend(['FVF','LVF','UVF'])
        plt.title('group_' + str(stimfreq[i-1]) + '_Hz' + '__electrode_' + eoi + '__reference [%s]' % ', '.join(map(str, ref)))
        plt.ylabel('$\mu$V')
        plt.xlabel('ms')
        
        yrange2=max(yrange)
        plt.ylim(-yrange2,yrange2)
        plt.xlim(tmin*1000,tmax*1000)
        plt.hlines(y=0, xmin=tmin*1000, xmax=tmax*1000)
        # plot stimulus onset vertical lines
        k=0
        #while (0.334 + k * stimdurations[i-1]) < 3.0:
        for k in range(0,len(single_event_timings[i-1])):
            plt.vlines(single_event_timings[i-1][k], -yrange2, yrange2, color='gray', linestyles='solid', label='onset') # color='red', linestyles='solid'
        #       plt.vlines(0.312 + k * stimdurations[i-1], -2, +2)
            #k += 1
        '''
        # plot vertical lines for all 100ms
        k=0
        while k <= 3000:
            plt.vlines(k, -yrange2, yrange2, color='gray', linestyles='dashed', label='')
        #       plt.vlines(0.312 + k * stimdurations[i-1], -2, +2)
            k += 100
        '''    
        if not os.path.isdir('group' + os.sep + 'plots'):
            os.makedirs('group' + os.sep + 'plots')
        plt.savefig('group' + os.sep + 'plots' + os.sep + 'group_1000dpi_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=1000)
        plt.savefig('group' + os.sep + 'plots' + os.sep + 'group_50dpi_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=50)
        plt.close('all')
        
def GROUP_TIME_seperate_freq_merged_fields_zoom_period(cm, eoi, np, Nsub, num_period, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax):
    # plot evoked of each frequency combined FOR ALL SUBJECTS AVERAGED
    colors = [ cm.Paired(x) for x in np.linspace(0.1, 0.4, 3) ]
    colors.append((0.5, 0.5, 0.5, 1.0)) # gray for the summation of lower und upper hemifield, if summation is selected
    tmin=322+num_period*single_period_times[num_period]
    for i in range(1,5): # Frequencies
        plt.figure()
        yrange=[]
        tmax=single_period_times[i-1]+tmin

        for j in range(1,4+summation): # hemifields + summation of upper and lower hemifield
        
            evoked_array = np.ndarray((0,0))
            for sub in range(4,11): # subjects ######## CHANGE HERE
                if sub<10:
                    subs = '0' + str(sub)
                else:
                    subs = str(sub)
                subject='RK_200' + subs
                if j<4: #for the full, upper and lower hemifield, not for the summation
                    tmp = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ') 
                else: # for the summation
                    tmp = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + 'summation__eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ') 
                evoked_array = np.append(evoked_array, tmp)    
            
            evoked_array = evoked_array.reshape((Nsub,3002))
            evoked_group = np.average(evoked_array, axis=0)
                
        
            yrange.append(np.ceil(max(max(evoked_group),abs(min(evoked_group)))*100000))
            
            # plot only one period, beginning at data point 322 until the first cycle (single_period_times[FreqBed])
            plt.plot([ x * 100000 for x in evoked_group[322:(322+single_period_times[i-1])] ],linewidth=2, color = colors[j-1])
   
        if summation:
            plt.legend(['FVF','LVF','UVF','SUM'])
        else:
            plt.legend(['FVF','LVF','UVF'])
        plt.title('group_' + str(stimfreq[i-1]) + '_Hz' + '__electrode_' + eoi + '__reference [%s]' % ', '.join(map(str, ref)))
        plt.ylabel('$\mu$V')
        plt.xlabel('ms')
        
        yrange2=max(yrange)
        plt.ylim(-yrange2,yrange2)
        plt.xlim(0,single_period_times[i-1])
        plt.hlines(y=0, xmin=0, xmax=single_period_times[i-1])
        
        # plot vertical lines for all 100ms
        k=0
        while k <= tmax:
            plt.vlines(k, -yrange2, yrange2, color='gray', linestyles='dashed', label='')
            k += 50
            
        if not os.path.isdir('group' + os.sep + 'plots'):
            os.makedirs('group' + os.sep + 'plots')
        plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_period_' + str(num_period) + '_evoked_group_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)))
        plt.close('all')

def GROUP_TIME_seperate_freq_merged_fields_zoom_period_all(cm, eoi, np, Nsub, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax):
    # plot evoked of each frequency combined FOR ALL SUBJECTS AVERAGED
    #colors = [ cm.Paired(x) for x in np.linspace(0.1, 0.4, 3) ]
    #colors.append((0.5, 0.5, 0.5, 1.0)) # gray for the summation of lower und upper hemifield, if summation is selected    
    colors = cm
    
    for i in range(1,5): # Frequencies
    
        # calculate mean of all periods
        if i==1: # if 2Hz, there are 4 periods
            Nperiods=4
        elif i==2: # if 4 Hz, there are 8 periods
            Nperiods=8
        elif i==3: # if 8 Hz, there are 15 periods
            Nperiods=15
        elif i ==4: # if 15 Hz, there are 30 periods
            Nperiods=30
    
        plt.figure()
        yrange=[]

        for j in range(1,4+summation): # hemifields + summation of upper and lower hemifield
        
            evoked_array = np.ndarray((0,0))
            for sub in range(4,11): # subjects ######## CHANGE HERE
                if sub<10:
                    subs = '0' + str(sub)
                else:
                    subs = str(sub)
                subject='RK_200' + subs
                if j<4: #for the full, upper and lower hemifield, not for the summation
                    tmp = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ') 
                else: # for the summation
                    tmp = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + 'summation__eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ') 
                evoked_array = np.append(evoked_array, tmp)    
            
            evoked_array = evoked_array.reshape((Nsub,3002)) # reshape to have not only one long line but data points from each subject in a line
            evoked_group = np.average(evoked_array, axis=0)
            
            # calculate average of all periods
            evoked_group_avg_period = np.ndarray((0,0))
            for k in range(0,Nperiods):
                evoked_group_avg_period = np.append(evoked_group_avg_period, evoked_group[(322+single_period_times[i-1]*(k)):(322+single_period_times[i-1]*(k+1))]) 
            evoked_group_avg_period = evoked_group_avg_period.reshape((Nperiods,single_period_times[i-1])) # reshape to have each period in a line so you can create an average in the next step
            evoked_group_avg_period = np.average(evoked_group_avg_period, axis=0)  # average over periods
        
            yrange.append(np.ceil(max(max(evoked_group_avg_period),abs(min(evoked_group_avg_period)))*100000))
 
            # plot mean of all periods
            plt.plot([ x * 100000 for x in evoked_group_avg_period ],linewidth=2, color = colors[j-1])
   
            # save data
            # normal data
            if j < 4:
                np.savetxt('group' + os.sep + 'data' + os.sep + 'evoked_zoom_avg_cond_' + str((i-1)*3 + j) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt',
                           evoked_group_avg_period ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='')
            else:
            # for summation data
                np.savetxt('group' + os.sep + 'data' + os.sep + 'evoked_zoom_avg_sum_' + str(i) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt',
                           evoked_group_avg_period ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='')
        
    
        if summation:
            plt.legend(['FVF','LVF','UVF','SUM'])
        else:
            plt.legend(['FVF','LVF','UVF'])
        plt.title('group_' + str(stimfreq[i-1]) + '_Hz' + '__electrode_' + eoi + '__reference [%s]' % ', '.join(map(str, ref)))
        plt.ylabel('$\mu$V')
        plt.xlabel('ms')
        
        yrange2=max(yrange)
        plt.ylim(-yrange2,yrange2)
        plt.xlim(0,single_period_times[i-1])
        plt.hlines(y=0, xmin=0, xmax=single_period_times[i-1])
        
        # plot vertical lines for all 100ms
        k=0
        while k <= single_period_times[i-1]:
            plt.vlines(k, -yrange2, yrange2, color='gray', linestyles='dashed', label='')
            k += 50
            
        if not os.path.isdir('group' + os.sep + 'plots'):
            os.makedirs('group' + os.sep + 'plots')
        plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_1000dpi_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=1000)
        plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_50dpi_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=50)
        plt.close('all')

''' group zoom but with same x-axis --> 1 cycle for 2 Hz, 2 Cycles for 4 Hz, 4 cycles for 8 Hz, and 8 cycles for 15 Hz '''
def GROUP_TIME_seperate_freq_merged_fields_zoom_period_all_sameX(cm, eoi, hemifields, np, Nsub, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax):
    # plot evoked of each frequency combined FOR ALL SUBJECTS AVERAGED
    #colors = [ cm.Paired(x) for x in np.linspace(0.1, 0.4, 3) ]
    #colors.append((0.5, 0.5, 0.5, 1.0)) # gray for the summation of lower und upper hemifield, if summation is selected    
    colors = cm
    
    for i in range(1,5): # Frequencies
    
        # calculate mean of all periods
        if i==1: # if 2Hz, there are 4 periods
            Nperiods=4
            Nrepeats_per_plot = 1
        elif i==2: # if 4 Hz, there are 8 periods
            Nperiods=8
            Nrepeats_per_plot = 2
        elif i==3: # if 8 Hz, there are 15 periods
            Nperiods=15
            Nrepeats_per_plot = 4
        elif i==4: # if 15 Hz, there are 30 periods
            Nperiods=30
            Nrepeats_per_plot = 8
        
        plt.figure()
        yrange=[]

        for j in range(1,4+summation): # hemifields + summation of upper and lower hemifield
        
            evoked_array = np.ndarray((0,0))
            for sub in range(4,11): # subjects ######## CHANGE HERE
                if sub<10:
                    subs = '0' + str(sub)
                else:
                    subs = str(sub)
                subject='RK_200' + subs
                if j<4: #for the full, upper and lower hemifield, not for the summation
                    tmp = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ') 
                else: # for the summation
                    tmp = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j) + 'summation__eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ') 
                evoked_array = np.append(evoked_array, tmp)    
            
            evoked_array = evoked_array.reshape((Nsub,3002)) # reshape to have not only one long line but data points from each subject in a line
            evoked_group = np.average(evoked_array, axis=0)
            
            # calculate average of all periods
            evoked_group_avg_period = np.ndarray((0,0))
            for k in range(0,Nperiods):
                evoked_group_avg_period = np.append(evoked_group_avg_period, evoked_group[(322+single_period_times[i-1]*(k)):(322+single_period_times[i-1]*(k+1))]) 
            evoked_group_avg_period = evoked_group_avg_period.reshape((Nperiods,single_period_times[i-1])) # reshape to have each period in a line so you can create an average in the next step
            evoked_group_avg_period = np.average(evoked_group_avg_period, axis=0)  # average over periods
        
            # save one avg period in np array
            if j<4:
                np.save('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_' + hemifields[j-1] + '_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), evoked_group_avg_period)
            else:
                np.save('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_sum_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), evoked_group_avg_period)
                
            yrange.append(np.ceil(max(max(evoked_group_avg_period),abs(min(evoked_group_avg_period)))*100000))
            
            # plot mean of all periods
            plt.plot([ x * 100000 for x in evoked_group_avg_period ]*Nrepeats_per_plot,linewidth=2, color = colors[j-1])
   
            # save data
            # normal data
            #if j < 4:
            #    np.savetxt('group' + os.sep + 'data' + os.sep + 'evoked_zoom_avg_cond_' + str((i-1)*3 + j) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt',
            #               evoked_group_avg_period ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='')
            #else:
            ## for summation data
            #    np.savetxt('group' + os.sep + 'data' + os.sep + 'evoked_zoom_avg_sum_' + str(i) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt',
            #               evoked_group_avg_period ,fmt='%.18e', delimiter=' ', newline='\n', header='', footer='')
        
    
        if summation:
            leg = plt.legend(['FVF','LVF','UVF','SUM'], loc='upper right')
        else:
            leg = plt.legend(['FVF','LVF','UVF'], loc='upper right')
        leg.set_zorder(99.)
        
        plt.title('group_' + str(stimfreq[i-1]) + '_Hz' + '__electrode_' + eoi + '__reference [%s]' % ', '.join(map(str, ref)))
        plt.ylabel('$\mu$V')
        plt.xlabel('ms')
        
        yrange2=max(yrange)
        plt.ylim(-yrange2,yrange2)
        plt.xlim(0,single_period_times[0]) # special for the zoom with same x axes, else it was [i-1] 
        plt.hlines(y=0, xmin=0, xmax=single_period_times[0]) # was also [i-1]
        
        # plot vertical lines for all 100ms
        k=0
        while k <= single_period_times[0]: # special for the zoom with same x axes, else it was [i-1] 
            plt.vlines(k, -yrange2, yrange2, color='gray', linestyles='dashed', label='')
            k += 50
        
        # special for sameX zoom: plot vertical bar at period transition to mask the offset and indicate transition
        for k in range(Nrepeats_per_plot-1):
            plt.vlines(single_period_times[i-1]*(k+1), -yrange2, yrange2, color='red', linestyles='solid', label='onset', zorder=98.0) # color='red', linestyles='solid'

        
        if not os.path.isdir('group' + os.sep + 'plots'):
            os.makedirs('group' + os.sep + 'plots')
        
        #plt.show() #plt.show makes it impossible to 
        plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_1000dpi_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=1000)
        plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_50dpi_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=50)
        plt.close('all')
    

''' group zoom but with same x-axis, plot only DIFFerence wave here !! '''
def GROUP_TIME_seperate_freq_merged_fields_zoom_period_all_sameX_DIFF(absolute, cm, eoi, np, mse, Nsub, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax):
    xranges = [275,275,275,275]
    for i in range(1,5): # Frequencies
    
        # calculate mean of all periods
        if i==1: # if 2Hz, there are 4 periods
            Nperiods=4
            Nrepeats_per_plot = 1
        elif i==2: # if 4 Hz, there are 8 periods
            Nperiods=8
            Nrepeats_per_plot = 2
        elif i==3: # if 8 Hz, there are 15 periods
            Nperiods=15
            Nrepeats_per_plot = 4
        elif i==4: # if 15 Hz, there are 30 periods
            Nperiods=30
            Nrepeats_per_plot = 8
        
        w, h = plt.figaspect(.3)
        fig = plt.figure(figsize=(w, h))
        
        yrange=[]

#        for j in range(1,4+summation): # hemifields + summation of upper and lower hemifield
        
        j1 = 2 #LVF
        j2 = 3 #UVF
            
        evoked_array1 = np.ndarray((0,0)) #LVF
        evoked_array2 = np.ndarray((0,0)) #UVF
        
        for sub in range(4,11): # subjects ######## CHANGE HERE
            if sub<10:
                subs = '0' + str(sub)
            else:
                subs = str(sub)
            subject='RK_200' + subs
            #if j<4: #for the full, upper and lower hemifield, not for the summation
            
            tmp1 = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j1) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ')  #UVF
            tmp2 = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j2) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ')  #LVF
            
            evoked_array1 = np.append(evoked_array1, tmp1)    
            evoked_array2 = np.append(evoked_array2, tmp2)  
            
        evoked_array1 = evoked_array1.reshape((Nsub,3002)) # reshape to have not only one long line but data points from each subject in a line
        evoked_array2 = evoked_array2.reshape((Nsub,3002)) # reshape to have not only one long line but data points from each subject in a line

        evoked_group1 = np.average(evoked_array1, axis=0)
        evoked_group2 = np.average(evoked_array2, axis=0)

        if mse == False:
            if absolute == False:
                evoked_diff  = np.subtract(evoked_group1,evoked_group2)
            else:
                evoked_diff  = abs(np.subtract(evoked_group1,evoked_group2))
        else: # squared error instead of substraction
            evoked_diff   = np.subtract(evoked_group1,evoked_group2) ** 2 # .mean(axis=0) # calculaes mse between two arrays along a given axis
        
        # calculate average of all periods
        evoked_group_avg_period = np.ndarray((0,0))
        for k in range(0,Nperiods):
            evoked_group_avg_period = np.append(evoked_group_avg_period, evoked_diff[(322+single_period_times[i-1]*(k)):(322+single_period_times[i-1]*(k+1))]) 
        evoked_group_avg_period = evoked_group_avg_period.reshape((Nperiods,single_period_times[i-1])) # reshape to have each period in a line so you can create an average in the next step
        evoked_group_avg_period = np.average(evoked_group_avg_period, axis=0)  # average over periods

        
        yrange.append(np.ceil(max(max(evoked_group_avg_period),abs(min(evoked_group_avg_period)))*100000))
        
        # save one avg period in np array
        np.save('group' + os.sep + 'data_arrays' + os.sep + 'zoom_DIFF_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), evoked_group_avg_period)
           
        # plot 
        plt.plot([ x * 100000 for x in evoked_group_avg_period ]*Nrepeats_per_plot,linewidth=2, color='black', linestyle='--')
        plt.axvspan(70, 90, alpha=0.5, color=[251/255,154/255,153/255])
        plt.axvspan(100, 120, alpha=0.5, color=[202/255,178/255,214/255])

        leg = plt.legend(['DIFF'], loc='upper right')
        leg.set_zorder(99.)
        
        plt.title('group_' + str(stimfreq[i-1]) + '_Hz' + '__electrode_' + eoi + '__reference [%s]' % ', '.join(map(str, ref)))
        plt.ylabel('$\mu$V')
        plt.xlabel('ms')
        
        yrange2=max(yrange)
        if absolute == False:
            plt.ylim(-yrange2,yrange2)
        else:
            plt.ylim(0,yrange2)
        plt.xlim(0,xranges[0]) # special for the zoom with same x axes, else it was [i-1] 
        plt.hlines(y=0, xmin=0, xmax=xranges[0]) # was also [i-1]
        
        # plot vertical lines for all 100ms
        k=0
        while k <= xranges[0]: # special for the zoom with same x axes, else it was [i-1] 
            #if absolute == False:
            plt.vlines(k, 0, yrange2, color='gray', linestyles='dashed', label='')
            #else:
            #    plt.vlines(k, 0, yrange2, color='gray', linestyles='dashed', label='')
            
            k += 50
        
        # special for sameX zoom: plot vertical bar at period transition to mask the offset and indicate transition
        for k in range(Nrepeats_per_plot-1):
            #if absolute==False:
            plt.vlines(single_period_times[i-1]*(k+1), -yrange2, yrange2, color='red', linestyles='solid', label='onset', zorder=98.0) # color='red', linestyles='solid'
            #else:
            #    plt.vlines(single_period_times[i-1]*(k+1), 0, yrange2, color='red', linestyles='solid', label='onset', zorder=98.0) # color='red', linestyles='solid'
            
        
        if not os.path.isdir('group' + os.sep + 'plots'):
            os.makedirs('group' + os.sep + 'plots')
        
        plt.show() #plt.show makes it impossible to 
        #plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_DIFF_1000dpi_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=1000)
        #plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_DIFF_50dpi_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=50)
        plt.close('all')
    
''' Single subject zoom but with same x-axis, plot only DIFFerence wave here !! '''
def SINGLE_SUBJECTS_TIME_seperate_freq_merged_fields_zoom_period_all_sameX_DIFF(cm, eoi, np, Nsub, os, plt, ref, single_event_timings, single_period_times, stimdurations, stimfreq, summation, tmin, tmax):
    
    for i in range(1,5): # Frequencies
    
        # calculate mean of all periods
        if i==1: # if 2Hz, there are 4 periods
            Nperiods=4
            Nrepeats_per_plot = 1
        elif i==2: # if 4 Hz, there are 8 periods
            Nperiods=8
            Nrepeats_per_plot = 2
        elif i==3: # if 8 Hz, there are 15 periods
            Nperiods=15
            Nrepeats_per_plot = 4
        elif i==4: # if 15 Hz, there are 30 periods
            Nperiods=30
            Nrepeats_per_plot = 8
        
        j1 = 2 #LVF
        j2 = 3 #UVF            

        for sub in range(4,11): # subjects ######## CHANGE HERE
            
            w, h = plt.figaspect(.3)
            fig = plt.figure(figsize=(w, h))        
            yrange=[]
            evoked_array1 = np.ndarray((0,0)) #LVF
            evoked_array2 = np.ndarray((0,0)) #UVF
            
            if sub<10:
                subs = '0' + str(sub)
            else:
                subs = str(sub)
            subject='RK_200' + subs
            
            tmp1 = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j1) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ')  #LVF
            tmp2 = np.loadtxt(subject + os.sep + 'data' + os.sep + 'evoked_cond_' + str((i-1)*3 + j2) + '_eoi_' + eoi + '_ref_[%s]' % ', '.join(map(str, ref)) + '.txt', delimiter=' ')  #UVF
            
            evoked_array1 = np.append(evoked_array1, tmp1)    
            evoked_array2 = np.append(evoked_array2, tmp2)  
            
#        evoked_array1 = evoked_array1.reshape((Nsub,3002)) # reshape to have not only one long line but data points from each subject in a line
#        evoked_array2 = evoked_array2.reshape((Nsub,3002)) # reshape to have not only one long line but data points from each subject in a line

#        evoked_group1 = np.average(evoked_array1, axis=0)
#        evoked_group2 = np.average(evoked_array2, axis=0)

            evoked_diff  = np.subtract(evoked_array1,evoked_array2)
    
            
            # calculate average of all periods
            evoked_avg_period = np.ndarray((0,0))
            for k in range(0,Nperiods):
                evoked_avg_period = np.append(evoked_avg_period, evoked_diff[(322+single_period_times[i-1]*(k)):(322+single_period_times[i-1]*(k+1))]) 
            evoked_avg_period = evoked_avg_period.reshape((Nperiods,single_period_times[i-1])) # reshape to have each period in a line so you can create an average in the next step
            evoked_avg_period = np.average(evoked_avg_period, axis=0)  # average over periods
    
            
            yrange.append(np.ceil(max(max(evoked_avg_period),abs(min(evoked_avg_period)))*100000))
            
            # plot 
            plt.plot([ x * 100000 for x in evoked_avg_period ]*Nrepeats_per_plot,linewidth=2, color='black', linestyle='--')
       
            leg = plt.legend(['DIFF'], loc='upper right')
            leg.set_zorder(99.)
            
            plt.title('subject_' + subs + '_' + str(stimfreq[i-1]) + '_Hz' + '__electrode_' + eoi + '__reference [%s]' % ', '.join(map(str, ref)))
            plt.ylabel('$\mu$V')
            plt.xlabel('ms')
            
            yrange2=max(yrange)
            plt.ylim(-yrange2,yrange2)
            plt.xlim(0,single_period_times[0]) # special for the zoom with same x axes, else it was [i-1] 
            plt.hlines(y=0, xmin=0, xmax=single_period_times[0]) # was also [i-1]
            
            # plot vertical lines for all 100ms
            k=0
            while k <= single_period_times[0]: # special for the zoom with same x axes, else it was [i-1] 
                plt.vlines(k, -yrange2, yrange2, color='gray', linestyles='dashed', label='')
                k += 50
            
            # special for sameX zoom: plot vertical bar at period transition to mask the offset and indicate transition
            for k in range(Nrepeats_per_plot-1):
                plt.vlines(single_period_times[i-1]*(k+1), -yrange2, yrange2, color='red', linestyles='solid', label='onset', zorder=98.0) # color='red', linestyles='solid'
    
            
            if not os.path.isdir('group' + os.sep + 'single'):
                os.makedirs('group' + os.sep + 'single')
            
            #plt.show() #plt.show makes it impossible to 
            plt.savefig('group' + os.sep + 'single' + os.sep + 'zoom_sameX_DIFF_1000dpi_subject_' + subs + '_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=1000)
            plt.savefig('group' + os.sep + 'single' + os.sep + 'zoom_sameX_DIFF_50dpi_subject_' + subs + '_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=50)
            plt.close('all')
    
def group_additivity_over_time(bin_length, eoi,hemifields, mne, np, os, plt, ref, stimfreq):
    # uses the zoomed_sameX data, and calculates how well LVF and UVF sum up to FVF over small time bins
    lenarray = 80
    data = np.ndarray((4,lenarray)) # 80 ist empriscih
    lengths = []
    for i in range(1,5): # Frequencies
        # load all hemifields + SUM
        full1  = np.load('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_' + hemifields[0] + '_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)) + '.npy')
#        lower = np.load('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_' + hemifields[1] + '_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)) + '.npy')
#        upper = np.load('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_' + hemifields[2] + '_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)) + '.npy')
        summe1 = np.load('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_' + 'sum'         + '_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)) + '.npy')
        if i==1: # if 2Hz, there are 4 periods
            Nrepeats_per_plot = 1
        elif i==2: # if 4 Hz, there are 8 periods
            Nrepeats_per_plot = 2
        elif i==3: # if 8 Hz, there are 15 periods
            Nrepeats_per_plot = 4
        elif i==4: # if 15 Hz, there are 30 periods
            Nrepeats_per_plot = 8
        summe = np.ndarray((0))
        full  = np.ndarray((0))
        for r in range(Nrepeats_per_plot):
            summe = np.append(summe, summe1)        
            full  = np.append(full, full1)

        # calculate distance bt. sum and full for small time bins
        diff = abs(full-summe)
        diff_ds = []
        # bin it
        steps = len(full) / bin_length
        for step in range(lenarray):
            if step <= int(np.floor(steps)-1):
                diff_ds.append(np.mean(diff[(step*bin_length):((step+1)*bin_length)])) #downsampled = ds
            else:
                if step == int(np.floor(steps)):
                    lengths.append(step-1)
                diff_ds.append(0)
        data[i-1] = np.asarray(diff_ds)
    data = data * 100000 # to bring it to muV
    data_sparse = np.delete(data, range(min(lengths),lenarray), axis=1)

    # define xticklabels
    xlabels = range(10,min(lengths)*bin_length+10,bin_length)
    xlabels = [str(x) for x in xlabels]
    for i in range(len(xlabels)): # replace every but 5 element with '' to make it sparse
        if i % 5 != 4:
            xlabels[i] = ''
    # plot results
    import seaborn as sns
    fig = sns.heatmap(data_sparse, cmap='gray', yticklabels=['2 Hz','4 Hz','8 Hz','15 Hz'], xticklabels = xlabels) 
    plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_additivity_1000dpi_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=1000)
    plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_additivity_50dpi__electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=50)
    
    return data, data_sparse


def group_similarity_over_time(bin_length, eoi,hemifields, mne, np, os, plt, ref, stimfreq):
    # uses the zoomed_sameX data, and calculates how well LVF and UVF and SUM correspond to FVF
    # similar to the group_additicity_over_time, but with all conditions compared to FVF, not just SUM
    lenarray = 80
    data_sum = np.ndarray((4,lenarray)) # 80 ist empriscih
    data_lvf = np.ndarray((4,lenarray)) # 80 ist empriscih
    data_uvf = np.ndarray((4,lenarray)) # 80 ist empriscih
    lengths = []
    for i in range(1,5): # Frequencies
        # load all hemifields + SUM
        full1  = np.load('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_' + hemifields[0] + '_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)) + '.npy')
        lower1 = np.load('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_' + hemifields[1] + '_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)) + '.npy')
        upper1 = np.load('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_' + hemifields[2] + '_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)) + '.npy')
        summe1 = np.load('group' + os.sep + 'data_arrays' + os.sep + 'zoom_' + str(int(np.round(stimfreq[i-1]))) + '_Hz_' + 'sum'         + '_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)) + '.npy')
        if i==1: # if 2Hz, there are 4 periods
            Nrepeats_per_plot = 1
        elif i==2: # if 4 Hz, there are 8 periods
            Nrepeats_per_plot = 2
        elif i==3: # if 8 Hz, there are 15 periods
            Nrepeats_per_plot = 4
        elif i==4: # if 15 Hz, there are 30 periods
            Nrepeats_per_plot = 8
        summe = np.ndarray((0))
        full  = np.ndarray((0))
        lower = np.ndarray((0))
        upper = np.ndarray((0))
        for r in range(Nrepeats_per_plot):
            summe = np.append(summe, summe1)        
            full  = np.append(full, full1)
            lower = np.append(lower, lower1)
            upper = np.append(upper, upper1)

        # calculate distance bt. sum and full for small time bins
        # for sum, lvf, and uvf, repectively
        diff_sum = abs(full-summe)
        diff_ds_sum = [] # downsampled
        diff_lvf = abs(full-lower)
        diff_ds_lvf = [] # downsampled
        diff_uvf = abs(full-upper)
        diff_ds_uvf = [] # downsampled
        
        # bin it
        steps = len(full) / bin_length
        for step in range(lenarray):
            if step <= int(np.floor(steps)-1):
                diff_ds_sum.append(np.mean(diff_sum[(step*bin_length):((step+1)*bin_length)])) #downsampled = ds
                diff_ds_lvf.append(np.mean(diff_lvf[(step*bin_length):((step+1)*bin_length)])) #downsampled = ds
                diff_ds_uvf.append(np.mean(diff_uvf[(step*bin_length):((step+1)*bin_length)])) #downsampled = ds
            else:
                if step == int(np.floor(steps)):
                    lengths.append(step-1)
                diff_ds_sum.append(0)
                diff_ds_lvf.append(0)
                diff_ds_uvf.append(0)
        data_sum[i-1] = np.asarray(diff_ds_sum)
        data_lvf[i-1] = np.asarray(diff_ds_lvf)
        data_uvf[i-1] = np.asarray(diff_ds_uvf)
    data_sum = data_sum * 100000 # to bring it to muV
    data_sparse_sum = np.delete(data_sum, range(min(lengths),lenarray), axis=1)
    data_lvf = data_lvf * 100000 # to bring it to muV
    data_sparse_lvf = np.delete(data_lvf, range(min(lengths),lenarray), axis=1)
    data_uvf = data_uvf * 100000 # to bring it to muV
    data_sparse_uvf = np.delete(data_uvf, range(min(lengths),lenarray), axis=1)

    # define xticklabels
    xlabels = range(10,min(lengths)*bin_length+10,bin_length)
    xlabels = [str(x) for x in xlabels]
    for i in range(len(xlabels)): # replace every but 5 element with '' to make it sparse
        if i % 5 != 4:
            xlabels[i] = ''
    # plot results
    import seaborn as sns
    fig = sns.heatmap(data_sparse_sum, cmap='gray', yticklabels=['2 Hz','4 Hz','8 Hz','15 Hz'], xticklabels = xlabels) 
    plt.title('SUM vs. FVF')
    plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_similarity_sum_1000dpi_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=1000)
    plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_similarity_sum_50dpi__electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=50)
    plt.show()
    plt.close()
    
    fig = sns.heatmap(data_sparse_lvf, cmap='gray', yticklabels=['2 Hz','4 Hz','8 Hz','15 Hz'], xticklabels = xlabels) 
    plt.title('LVF vs. FVF')
    plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_similarity_lvf_1000dpi_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=1000)
    plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_similarity_lvf_50dpi__electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=50)
    plt.show()
    plt.close()
    
    fig = sns.heatmap(data_sparse_uvf, cmap='gray', yticklabels=['2 Hz','4 Hz','8 Hz','15 Hz'], xticklabels = xlabels) 
    plt.title('UVF vs. FVF')
    plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_similarity_uvf_1000dpi_electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=1000)
    plt.savefig('group' + os.sep + 'plots' + os.sep + 'zoom_sameX_similarity_uvf_50dpi__electrode_' + eoi + '_reference [%s]' % ', '.join(map(str, ref)), dpi=50)
    plt.show()
    plt.close()
    
    #return data, data_sparse

##########################
### FREQUENCY ANALYSIS ###
##########################

def freq_analysis(eoi, eoi_idx, epochs, hemifields2, mne, np, os, parent_dir, plt, raw, ref, stimfreq, stimfreq2, subject, subs):
    #epochs.plot_psd(fmin=2., fmax=40.) #Let’s first check out all channel types by averaging across epochs.
    montage = mne.channels.read_montage(parent_dir + '/electrode_positions/M1_ThetaPhi_32_electrodes.txt', ch_names=None, path=None, unit='m', transform=False)  
    
    raw.set_montage(montage)
    epochs.set_montage(montage)
    
#    epochs.plot_psd_topomap(bands = [(4, 8, 'Theta')], tmin = 0.3, tmax = 2.3,
#                                     ch_type='eeg', normalize=True)
    
    # time-frequency representations (TFRs) from Epochs.
    # We’ll look at power and intertrial coherence (ITC)
    for i in range(0,12):

        freqs = np.arange(1.875, 48, 1.875)
        n_cycles = freqs / 2.  # different number of cycle per frequency
        power, itc =  mne.time_frequency.tfr_morlet(epochs['cond_'+str(i+1)], freqs=freqs,
                            n_cycles=n_cycles, use_fft=True, #picks = eoi_idx, 
                            return_itc=True, decim=3, n_jobs=1)
        power.plot([eoi_idx], baseline=(-0.5, 0), mode='logratio')
        plt.title('TFR_subject_' + subs + '_' + str(int(np.round(stimfreq2[i]))) + '_Hz_' + hemifields2[i] + '_visual_field__electrode' + eoi + '_reference [%s]' % ', '.join(map(str, ref)))
        plt.savefig(subject + os.sep + 'plots' + os.sep + 'TFR_subject_' + subs + '_' + str(int(np.round(stimfreq2[i]))) + '_Hz_' + hemifields2[i] + '_visual_field__electrode' + eoi + '_reference [%s]' % ', '.join(map(str, ref)))
        plt.show()
        plt.close('all')
    return epochs, montage, raw
        
    '''
        if np.isnan(power.data).any():
            power.data = np.nan_to_num(power.data)
            
        # workaround for the bug
        threshold_for_bug = 0.0001#00001 # could be any value, ex numpy.min
        power.data[power.data < threshold_for_bug] = threshold_for_bug
        power.plot_topomap(ch_type='eeg', tmin=0.312, tmax=2.3, fmin=8, fmax=12,
                   baseline=(-0.5, 0), mode='logratio', axes=axis[0],
                   title='Alpha', vmax=0.45, show=True)
        
        power.plot_topomap(ch_type='eeg', tmin=0.312, tmax=2.3, fmin=8, fmax=12,
                   #baseline=(-0.5, 0), mode='logratio', axes=axis[0],
                   title='Alpha', vmax=0.45, show=True)

        power.plot_topomap(ch_type='eeg', tmin=0.312, tmax=2.3, fmin=8, fmax=12,
                   baseline=(-0.5, 0), #mode='logratio', axes=axis[0],
                   title='Alpha', vmax=0.45, show=False)
    
    fig, axis = plt.subplots(1, 2, figsize=(7, 4))
    power.plot_topomap(ch_type='eeg', tmin=0.5, tmax=1.5, fmin=8, fmax=12,
                       baseline=(-0.5, 0), mode='logratio', axes=axis[0],
                       title='Alpha', vmax=0.45, show=False)
    power.plot_topomap(ch_type='eeg', tmin=0.5, tmax=1.5, fmin=13, fmax=25,
                       baseline=(-0.5, 0), mode='logratio', axes=axis[1],
                       title='Beta', vmax=0.45, show=False)
    mne.viz.tight_layout()
    plt.show()
    '''

    '''
    evoked1 = epochs['cond_1'].average()
    evoked2 = epochs['cond_2'].average()
    evoked3 = epochs['cond_3'].average()
    evoked4 = epochs['cond_4'].average()
    evoked5 = epochs['cond_5'].average()
    evoked6 = epochs['cond_6'].average()
    evoked7 = epochs['cond_7'].average()
    evoked8 = epochs['cond_8'].average()
    evoked9 = epochs['cond_9'].average()
    evoked10 = epochs['cond_10'].average()
    evoked11 = epochs['cond_11'].average()
    evoked12 = epochs['cond_12'].average()
    '''
    
######## MISC

def pp_evoked_oz_merged(division, evoked_oz_cond, math, np, os, plt, stimfreq, subject, tmin, tmax):
    ########### THIS is a 3x4 plot !!! ######################
    #plt.close('all')
    colorQ12=[[166,206,227], [31,120,180], [178,223,138], [51,160,44], [251,154,153], [227,26,28], [253,191,111], [255,127,0], [202,178,214], [106,61,154], [255,255,153], [177,89,40]]
    colorQ12=[[j/255 for j in i] for i in colorQ12];
    x = np.linspace(tmin, tmax, len(evoked_oz_cond[0][0]))
    # Four axes, returned as a 2-d array
    counter=0
    f, axarr = plt.subplots(3, 4, figsize=(20, 60))
    #plt.title('Evoked responses')
    #erst horizontal, dann vertical zaehlen
    for i in range(0,3):
        for j in range(0,4):
            axarr[i, j].plot(x, evoked_oz_cond[counter][0]*1000000, color='k')#color=colorQ12[counter])
            #plt.hlines(y=0, xmin=tmin, xmax=tmax, hold=True)
            
            # axarr[i, j].
            axarr[i, j].set_ylim([-50, 50])
            axarr[i, j].axhline(y=0.0,xmin=tmin,xmax=tmax,c="blue",linewidth=0.5,zorder=0)
            axarr[i, j].axvline(x=0.3,c="blue",linewidth=0.5,zorder=0)
            # condition specific vertical lines for stim onsets
            n=1
            while (0.3+1/stimfreq[j]*n) < tmax:
                axarr[i, j].axvline(x=0.3+1/stimfreq[j]*n,c="blue",linewidth=0.5,zorder=0)
                n += 1
            #axarr[j, j].ylim([-0.001,0.001])
            if counter<4: #erste reihe ueberschriften einfuegen
                axarr[i, j].set_title(str(stimfreq[j])+' Hz')
            #    axarr[i, j].
            #if counter%4==0: #erste Spalte wieder ueberschriften der Bedingungen einfuehren
            #    axarr[i, j].
            counter += 1
    # change figure size !!
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    
    # shared y label
    f.text(0.06, 0.5, '$\mu$V', ha='center', va='center', rotation = 'vertical')
    # y label for conditions:
    f.text(0.95, 0.5, 'lower hemifield - upper hemifield - full visual field', ha='center', va='center', rotation = 'vertical')
    plt.savefig(subject + os.sep + 'plots' + os.sep + 'evoked_3x4.png', dpi=300)#1000)
    return



def pp_evoked_oz_seperate(eoi_idx, epochs, hemifields2, Ncond, np, os, picks, plt, screenfreq, stimfreq2, subject, time_per_stim_and_ISI2, tmin, tmax):
    # create epochs for each condition
    evoked_oz_cond=[]
    #evoked_oz_cond=np.array((12,len(epochs['cond_1'].average(picks=[29]).data)))
    for i in range(0,Ncond):
        evoked_oz_cond.append(epochs[i].average(picks=[eoi_idx]).data)
    for ncond in range(0,12):
        evoked_oz_cond_tmp=epochs['cond_'+str(ncond+1)].average(picks=[eoi_idx]) # 29='Oz'
        ### plotting
        ########### single conditions - time analysis ###########
        #plt.close('all')
        fig1 = evoked_oz_cond_tmp.plot()
        plt.axhline(y=0.0,xmin=tmin,xmax=tmax,c="black",linewidth=0.5,zorder=0)
        plt.axvline(x=300,c="black",linewidth=0.5,zorder=0)
        n = 0
    #    while (0.3 + 1/stimfreq2[ncond]*n ) < tmax:
        while (np.ceil(0.3*screenfreq)/screenfreq + np.ceil(time_per_stim_and_ISI2[ncond]*screenfreq)/screenfreq * n) < tmax:
            #plt.axvline(x=300+1000/stimfreq2[ncond]*n,c="black",linewidth=0.5,zorder=0)
            plt.axvline(x= (np.ceil(0.3*screenfreq)/screenfreq*1000 + np.ceil(time_per_stim_and_ISI2[ncond]*screenfreq)/screenfreq * n * 1000) ,c="black",linewidth=0.5,zorder=0)
    #        plt.axvline(x= (np.ceil(0.3*screenfreq)/screenfreq + np.ceil(time_per_stim_and_ISI2[0]*screenfreq)/screenfreq * n) ,c="black",linewidth=0.5,zorder=0)
            n += 1
        plt.title('Oz electrode, ' + str(stimfreq2[ncond]) + ' Hz, ' + hemifields2[ncond] + ' hemifield' )    
        plt.show()
        plt.savefig(subject + os.sep + 'plots' + os.sep + 'evoked3_cond' + str(ncond+1) + '.png', dpi=300)#1000)
    return evoked_oz_cond

