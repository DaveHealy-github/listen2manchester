#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 13:30:59 2021

@author: davidhealy
"""
import numpy as np 
import matplotlib.pyplot as plt

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*24   # length in seconds of data to plot after origin time
Flow = 4.
Fhigh = 20.  

#   now get some station meta data, esp. lat & long 
#namStation1 = 'R9E71'    #   RHS Bridgewater, Worsley   
#namStation1 = 'R3FEA'    #   Oswald Road, Chorlton   
#namStation1 = 'R72C7'    #   Oswald Road, Chorlton   
namStation1 = 'R6055'    #   GeoMon, Amlwch   
netStation1 = 'AM'
locStation1 = '00'
chaStation1 = 'EHZ'      #   vertical component of geophone 
client1 = Client('RASPISHAKE')
metadata1 = client1.get_stations(network=netStation1, station=namStation1, location=locStation1, channel=chaStation1, level='resp')
idSEED1 = netStation1 + '.' + namStation1 + '.' + locStation1 + '.' + chaStation1 
latStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['latitude']
lonStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['longitude']
eleStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['elevation']

#   start of 'event' 
#EQ_TIME1 = "2021-12-05T12:33:00"
EQ_TIME1 = "2021-12-16T10:35:00"
dtEvent1 = UTCDateTime(EQ_TIME1)
sEventName1 = 'Raw data (response removed, demeaned)'
sEventName2 = 'Filtered data'

#   get the data
#   plot seismograms  
t1 = dtEvent1 - T_START
t2 = dtEvent1 + T_END

# Download and filter data for station1
st1 = client1.get_waveforms(netStation1, namStation1, locStation1, chaStation1,
                          starttime=t1, endtime=t2, attach_response=True)
st1.merge(method=0, fill_value='latest')
st1.detrend(type="demean")
st1.remove_response()
st1.trim(t1, t2)

# Now plot the waveform data
fig = plt.figure(figsize=(12,4))
axs = fig.gca()
st1[0].spectrogram(log=False, wlen=4, axes=axs, dbscale=True, 
                   cmap='plasma', show=True, samp_rate=st1[0].stats.sampling_rate)
mappable = axs.images[0]
#mappable = axs.collections[0]
axs.set_ylabel('Frequency, Hz')
axs.set_xlabel('Time from start, s')
#axs.set_ylim(0.1, 2.)
for im in fig.gca().get_images():
    im.set_clim(-160,-100)
plt.colorbar(mappable=mappable, ax=axs) 

plt.tight_layout() 
plt.savefig('plotDopplerSpectrogramR6055_16Dec21_microseism.png', dpi=300)
