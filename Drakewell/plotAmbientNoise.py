#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 13:42:35 2021

@author: davidhealy
"""
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd 

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy import signal 
import obspy as obspy 

#   now get some station meta data, esp. lat & long 
namStation1 = 'R72C7'    #   in my office at Aberdeen 
netStation1 = 'AM'
locStation1 = '00'
chaStation1 = 'EHZ'      #   vertical component of geophone 
client1 = Client('RASPISHAKE')
metadata1 = client1.get_stations(network=netStation1, station=namStation1, location=locStation1, channel=chaStation1, level='resp')
idSEED1 = netStation1 + '.' + namStation1 + '.' + locStation1 + '.' + chaStation1 
latStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['latitude']
lonStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['longitude']
eleStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['elevation']

#   time window of interest  
T_END = 2*24*60*60      # length in seconds of data to plot after origin time
Flow = 4.
Fhigh = 14.  
t1string = "2021-09-08"
dt1 = UTCDateTime(t1string)
dt2 = dt1 + T_END
sEventName1 = "Last few weeks from: " + t1string 
datelist = pd.date_range(dt1.datetime, dt2.datetime, freq="D")

#   Download and filter data for station1, by day and then splice 
#   Otherwise, memory overload for more than 3-4 days of data in one go  
istream = 0 
for day in datelist:
    istream += 1 
    st = client1.get_waveforms(netStation1, namStation1, 
                                locStation1, chaStation1,
                                UTCDateTime(day)-1801, 
                                UTCDateTime(day)+86400+1801, 
                                attach_response=True)
    print(st) 
    if istream == 1:
        stall = st.copy()
    else:
        stall += st 
    
stall.merge(method=0, fill_value='latest')
stall.detrend(type="demean")
stall.remove_response()
stall.resample(1000)
npts = stall[0].stats.npts
samprate = stall[0].stats.sampling_rate
stall.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
data_envelope = signal.filter.envelope(stall[0].data)
t = np.arange(0, npts / samprate, 1 / samprate)

# Now plot the waveform data
fig, axs = plt.subplots(figsize=(12,6))
axs.plot(t, stall[0].data, linewidth=1,
        color="darkred")
axs.plot(t, data_envelope, ':r')
axs.grid(True)
axs.set_xlabel("Date & time since start")
axs.set_title("{:} - RaspberryShake {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    sEventName1, stall[0].stats.network, stall[0].stats.station, stall[0].stats.location,
    stall[0].stats.channel, Flow, Fhigh))
axs.set_ylabel("Ground velocity (mm/s)")
#axs.set_xlim(T_START+1, T_END)
