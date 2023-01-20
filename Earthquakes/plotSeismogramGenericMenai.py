#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:58:30 2021

@author: davidhealy
"""
import matplotlib.pyplot as plt

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

T_START = -10     # length in seconds of data to plot before origin time
T_END = 60*2      # length in seconds of data to plot after origin time
Flow = 1.
Fhigh = 20.   

#   now get some station meta data, esp. lat & long 
#namStationList = ['R3FEA', 'R6C8A', 'RAEED', 'R36CB', 'RE28A', 'R6055']
namStationList = ['RE28A']
#namStationList = ['R6055', 'R72C7', 'RE28A']
#labelStationList = ['Chorlton', 'MMU', 'UoM', 'Cheadle Hulme', 'Menai', 'Amlwch']
labelStationList = ['Menai Bridge']
#labelStationList = ['Amlwch', 'Aberdeen', 'Menai']
netStation = 'AM'
locStation = '00'
#chaStation = 'EHZ'      #   vertical component of geophone 
chaStationList = ['EHZ', 'EHE', 'EHN']      #   vertical component of geophone 
client1 = Client('RASPISHAKE')

#   event datetime, e.g. from USGS or BGS
# EQ_TIME = "2022-03-16T14:36:33"
# EQ_NAME = "M7.3 Honshu 2022-03-16"
# EQ_TIME = "2022-03-16T23:15:45"
# EQ_NAME = "M5.8 Iran 2022-03-16"
# EQ_TIME = "2022-03-22T17:41:38"
# EQ_NAME = "M6.7 Taiwan 2022-03-22"
# EQ_TIME = "2022-03-30T20:56:58"
# EQ_NAME = "M6.9 New Caledonia 2022-03-30"
# EQ_TIME = "2022-05-29T20:40:37.9"
# EQ_NAME = "M2.3 Sale, GM 2022-05-29"
EQ_TIME = "2022-05-30T14:36:57.7"
EQ_NAME = "M3.8 Shropshire, GM 2022-05-30"
dtEvent = UTCDateTime(EQ_TIME)
t1start = dtEvent + T_START
t1end = dtEvent + T_END

# Now plot the waveform data
fig, axs = plt.subplots(len(chaStationList), 1, figsize=(12,12))

i = 0 
for namStation in namStationList:
    
    for chaStation in chaStationList:
        
        # Download and filter data for station1
        st1 = client1.get_waveforms(netStation, namStation, locStation, chaStation,
                              starttime=t1start, endtime=t1end, attach_response=True)
        st1.merge(method=0, fill_value='latest')
        st1.detrend(type="demean")
        st1.remove_response()
        st1f = st1.copy()
        st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
        st1.trim(t1start, t1end)
        st1f.trim(t1start, t1end)
    
        axs[i].plot(st1f[0].times(reftime=dtEvent+1), st1f[0].data*1000, linewidth=1,
                color="darkred")
        axs[i].grid(True)
        axs[i].set_xlabel("Time after earthquake (s)")
        axs[i].set_title("{:} {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz - {:}".format(
            EQ_NAME, st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
            st1f[0].stats.channel, Flow, Fhigh, labelStationList[0]))
        axs[i].set_ylabel("Ground velocity (mm/s)")
        # if i >= 4: 
        #     axs[i].set_xlim(15, 45)
        #     axs[i].set_ylim(-0.0025, 0.0025)
        # else: 
        #     axs[i].set_xlim(0, 25)
        #     axs[i].set_ylim(-0.02, 0.02)
        axs[i].set_xlim(10, 100)
        axs[i].set_ylim(-0.05, 0.05)
    
        i = i + 1

plt.tight_layout()
fn = EQ_NAME + ".png" 
plt.savefig(fn, dpi=300)
