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
T_END = 60*20    # length in seconds of data to plot after origin time
Flow = 0.5
Fhigh = 2.   

#   now get some station meta data, esp. lat & long 
namStation1 = 'RE28A'    #   Oswald Road school, Chorlton   
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
EQ_TIME1 = "2022-03-16T14:36:33"
dtEvent1 = UTCDateTime(EQ_TIME1)
t1start = dtEvent1 - T_START
t1end = dtEvent1 + T_END

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

# Now plot the waveform data
fig, axs = plt.subplots(1, 1, figsize=(12,5))

axs.plot(st1f[0].times(reftime=dtEvent1+1), st1f[0].data*1000, linewidth=1,
        color="darkred")
axs.grid(True)
axs.set_xlabel("Time after earthquake (s)")
axs.set_title("M7.3 Honshu 2022-03-16 {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
    st1f[0].stats.channel, Flow, Fhigh))
axs.set_ylabel("Ground velocity (mm/s)")
axs.set_xlim(8*60, 18*60)

plt.tight_layout() 
plt.savefig('plotRE28A Menai - M7.3 Honshu 2022-03-16.png', dpi=300)
