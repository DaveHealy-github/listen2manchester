#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:58:30 2021

@author: davidhealy
"""
import numpy as np  
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

Flow = 0.1
Fhigh = 1.5  

client=Client("IRIS") # set IRIS as the place you're going to download data from
t_eq=UTCDateTime("2021-10-12T09:24:03") # set a date and time for the earthquake
tstart=t_eq+60*4 #start downloading from 300s after the earthquake
tend=t_eq+60*13 #stop downloading data this number of seconds after the earthquake

st=client.get_waveforms("GB","STNC", "*", "*HZ", tstart, tend, attach_response=True) 
st.merge(method=0, fill_value='latest')

std = st.copy()
std.detrend(type="demean")

strem = std.copy()
strem.remove_response()

stf = strem.copy() 
stf.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
stf.trim(tstart, tend)

# Now plot the waveform data
fig = plt.figure(figsize=(12,5)) 
plt.plot(st[0].times(reftime=tstart+1), st[0].data, linewidth=1, color="darkblue" )
plt.grid(True)
plt.xlabel("Time after earthquake (s)")
plt.xlim(0., np.max(st[0].times(reftime=tstart+1)))
plt.ylabel("Counts - raw")

fig = plt.figure(figsize=(12,5)) 
plt.plot(std[0].times(reftime=tstart+1), std[0].data, linewidth=1, color="darkblue" )
plt.grid(True)
plt.xlabel("Time after earthquake (s)")
plt.xlim(0., np.max(std[0].times(reftime=tstart+1)))
plt.ylabel("Velocity, m/s - detrend(demean)")

fig = plt.figure(figsize=(12,5)) 
plt.plot(strem[0].times(reftime=tstart+1), strem[0].data, linewidth=1, color="darkblue" )
plt.grid(True)
plt.xlabel("Time after earthquake (s)")
plt.xlim(0., np.max(strem[0].times(reftime=tstart+1)))
plt.ylabel("Velocity, m/s - remove_response()")

fig = plt.figure(figsize=(12,5)) 
plt.plot(stf[0].times(reftime=tstart+1), stf[0].data*1.e6, linewidth=1, color="darkblue" )
plt.grid(True)
plt.xlabel("Time after earthquake (s)")
plt.xlim(0., np.max(stf[0].times(reftime=tstart+1)))
plt.ylabel(r"Velocity, $\mu$m/s - filtered (bandpass)")
