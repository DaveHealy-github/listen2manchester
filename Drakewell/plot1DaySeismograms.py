#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 10:28:36 2022

@author: davidhealy
"""
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import pandas as pd 
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*60*24     # length in seconds of data to plot after origin time
Flow = 0.7
Fhigh = 1.   

#   now get some station meta data, esp. lat & long 
namStationList = ['R4DD7', 'R3FEA', 'R6C8A', 'R36CB', 'LBWR']
labelStationList = ['Firs', 'Chorlton', 'MMU', 'Cheadle Hulme', 'Ladybower']
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client1 = Client('RASPISHAKE')

EQ_TIME = "2022-06-26T00:00:00"
dtEvent = UTCDateTime(EQ_TIME)
t1start = dtEvent + T_START
t1end = dtEvent + T_END

# Now plot the waveform data
#fig, axs = plt.subplots(len(namStationList), 1, figsize=(3,6))
fig, axs = plt.subplots(figsize=(4,5))

i = 0 
for namStation in namStationList:

    if "LBWR" in namStation:
        client1=Client('http://eida.bgs.ac.uk')
        netStation = 'GB'
        locStation = '00' 
        chaStation = 'HHZ'
        
    # Download and filter data for station1
    st1 = client1.get_waveforms(netStation, namStation, locStation, chaStation,
                          starttime=t1start, endtime=t1end, attach_response=True)
    st1.merge(method=0, fill_value='latest')
    st1.detrend(type="demean")
    st1.remove_response()
    st1f = st1.copy()
#    st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
    st1f.filter("highpass", freq=Fhigh, corners=2)
    st1.trim(t1start, t1end)
    st1f.trim(t1start, t1end)
    st1f.decimate(5)

    #   build a dataframe of velocity amplitudes and their datetimes 
    #   this helps to improve the x-axis formatting...         
    index = pd.DatetimeIndex([(dtEvent + ns).datetime for ns in st1f[0].times()])
    vampl = pd.DataFrame(st1f[0].data*1000, index=index)
    vampl_clip = np.copy(vampl[0])
    vampl_clip[abs(vampl_clip) > 0.05] = 0.
    
    # axs[i].plot(st1f[0].times(reftime=dtEvent), st1f[0].data*1000, linewidth=1,
    #         color="darkred")
    offset = float(i * 0.08) 
    axs.plot(vampl.index, vampl_clip + offset, linewidth=0.5, label=namStation) 
    axs.grid(True)
    axs.set_xlabel("Hour of day on 26 June 2022")
#    axs[i].set_title("{:}.{:}.{:}.{:} - highpass filter: {:} Hz - {:}".format(
#        st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
#        st1f[0].stats.channel, Fhigh, labelStationList[i]))
#    axs[i].set_ylabel("Ground velocity (mm/s)")
    axs.set_ylim(0.38, -0.04)
#    axs.legend()

    i = i + 1

plt.yticks([0.0, 0.08, 0.16, 0.24, 0.32], namStationList)
date_form = DateFormatter("%H")
axs.xaxis.set_major_formatter(date_form)
axs.autoscale(enable=True, axis='x', tight=True)
plt.tight_layout()
fn = "1DaySeismograms_v2.png" 
plt.savefig(fn, dpi=300)
