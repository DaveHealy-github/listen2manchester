#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:58:30 2021

@author: davidhealy
"""
import matplotlib.pyplot as plt
import pandas as pd 

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*1.5     # length in seconds of data to plot after origin time
Flow = 3.
Fhigh = 20.   

#   now get some station meta data, esp. lat & long 
namStationList = ['R4DD7', 'RE28A', 'LBWR']
#namStationList = ['RE28A']
#namStationList = ['R6055', 'R72C7', 'RE28A']
labelStationList = ['Fallowfield', 'Menai Bridge', 'BGS Ladybower']
#labelStationList = ['Menai Bridge']
#labelStationList = ['Amlwch', 'Aberdeen', 'Menai']
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
#chaStationList = ['EHZ', 'EHE', 'EHN']      #   vertical component of geophone 
client1 = Client('RASPISHAKE')

#   event datetime, e.g. from USGS or BGS
#EQ_TIME = "2022-06-21T20:54:36"
#EQ_NAME = "M5.9 Khost, Afghanistan, 2022-06-21"
#EQ_TIME = "2022-06-08T00:55:47"
#EQ_NAME = "M6.5 Tarauaca, Brazil 2022-06-08"
#EQ_TIME = "2022-05-29T20:40:37.9"
#EQ_NAME = "M2.3 Sale, GM 2022-05-29"
EQ_TIME = "2022-07-14T14:43:33.3"
EQ_NAME = "M2.4 Leigh, GM 2022-07-14"
dtEvent = UTCDateTime(EQ_TIME)
t1start = dtEvent + T_START
t1end = dtEvent + T_END

# Now plot the waveform data
fig, axs = plt.subplots(len(namStationList), 1, figsize=(10,10))

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
    st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
    st1.trim(t1start+1, t1end)
    st1f.trim(t1start+1, t1end)

    #   build a dataframe of velocity amplitudes and their datetimes 
    #   this helps to improve the x-axis formatting...         
    index = pd.DatetimeIndex([(dtEvent + ns).datetime for ns in st1f[0].times()])
    vampl = pd.DataFrame(st1f[0].data*1000, index=index)

    if i == 2:
        # axs[i].plot(st1f[0].times(reftime=dtEvent+1), st1f[0].data*1000, linewidth=1,
        #             color="blue")
        axs[i].plot(vampl.index, vampl[0], linewidth=1,
                    color="blue")
    else:  
        # axs[i].plot(st1f[0].times(reftime=dtEvent+1), st1f[0].data*1000, linewidth=1,
        #             color="darkred")
        axs[i].plot(vampl.index, vampl[0], linewidth=1,
                    color="darkred")
        
    axs[i].grid(True)
    axs[i].set_xlabel("Time (HH:MM:SS)")
    axs[i].set_title("{:} {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz - {:}".format(
        EQ_NAME, st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
        st1f[0].stats.channel, Flow, Fhigh, labelStationList[i]))
    axs[i].set_ylabel("Ground velocity (mm/s)")
#    axs[i].set_xlim(0, 100)
    if i == 2:
        axs[i].set_ylim(-0.0015, +0.0015)
    # else:
    #     axs[i].set_ylim(-0.015, +0.015)

    i = i + 1

plt.tight_layout()
fn = "QuakeSeismogram_Paper.png" 
plt.savefig(fn, dpi=300)
