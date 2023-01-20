#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:58:30 2021

@author: davidhealy
"""
import matplotlib.pyplot as plt

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*30    # length in seconds of data to plot after origin time
Flow = 0.7
Fhigh = 2.   

#   now get some station meta data, esp. lat & long 
namStation1 = 'R3FEA'    #   Oswald Road school, Chorlton   
netStation1 = 'AM'
locStation1 = '00'
chaStation1 = 'EHZ'      #   vertical component of geophone 
client1 = Client('RASPISHAKE')
metadata1 = client1.get_stations(network=netStation1, station=namStation1, location=locStation1, channel=chaStation1, level='resp')
idSEED1 = netStation1 + '.' + namStation1 + '.' + locStation1 + '.' + chaStation1 
latStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['latitude']
lonStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['longitude']
eleStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['elevation']

#   event datetime, e.g. from USGS or BGS
#   M6.8 Alaska
EQ_TIME1 = "2022-01-11T11:35:46"
dtEvent1 = UTCDateTime(EQ_TIME1)
t1start = dtEvent1 - T_START
t1end = dtEvent1 + T_END

#   M6.6 Cyprus
EQ_TIME2 = "2022-01-11T01:07:47"
dtEvent2 = UTCDateTime(EQ_TIME2)
t2start = dtEvent2 - T_START
t2end = dtEvent2 + T_END

# Download and filter data for station1
st1 = client1.get_waveforms(netStation1, namStation1, locStation1, chaStation1,
                          starttime=t1start, endtime=t1end, attach_response=True)
st1.merge(method=0, fill_value='latest')
st1.detrend(type="demean")
st1.remove_response()
st1f = st1.copy()
st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
st1.trim(t1start, t1end)
st1f.trim(t1start, t1end)

st2 = client1.get_waveforms(netStation1, namStation1, locStation1, chaStation1,
                          starttime=t2start, endtime=t2end, attach_response=True)
st2.merge(method=0, fill_value='latest')
st2.detrend(type="demean")
st2.remove_response()
st2f = st2.copy()
st2f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
st2.trim(t2start, t2end)
st2f.trim(t2start, t2end)

# Now plot the waveform data
fig, axs = plt.subplots(2, 1, figsize=(12,8))

axs[0].plot(st1f[0].times(reftime=dtEvent1+1), st1f[0].data*1000, linewidth=1,
        color="darkred")
axs[0].grid(True)
axs[0].set_xlabel("Time after earthquake (s)")
axs[0].set_title("M6.8 Alaska 2022-01-11 11:35 {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
    st1f[0].stats.channel, Flow, Fhigh))
axs[0].set_ylabel("Ground velocity (mm/s)")
axs[0].set_xlim(500, 1250)

axs[1].plot(st2f[0].times(reftime=dtEvent2+1), st2f[0].data*1000, linewidth=1,
        color="darkblue")
#ymins2f, ymaxs2f = axs[1].get_ylim()
axs[1].grid(True)
axs[1].set_xlabel("Time after earthquake (s)")
axs[1].set_title("M6.6 Cyprus 2022-01-11 01:07 {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    st2f[0].stats.network, st2f[0].stats.station, st2f[0].stats.location,
    st2f[0].stats.channel, Flow, Fhigh))
axs[1].set_ylabel("Ground velocity (mm/s)")
axs[1].set_xlim(10, 1250)

plt.tight_layout() 
plt.savefig('plotR3FEA.png', dpi=300)
