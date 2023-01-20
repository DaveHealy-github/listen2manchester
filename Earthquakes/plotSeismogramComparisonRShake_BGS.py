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
T_END = 60*2.5    # length in seconds of data to plot after origin time
Flow = 0.7
Fhigh = 2.   

#   now get some station meta data, esp. lat & long 
#namStation1 = 'R9E71'    #   RHS Bridgewater 
#namStation1 = 'R6055'    #   GeoMon, Amlwch 
#namStation1 = 'R72C7'    #   DH office, Aberdeen 
namStation1 = 'RE28A'    #   Menai Bridge, Anglesey 
netStation1 = 'AM'
locStation1 = '00'
chaStation1 = 'EHZ'      #   vertical component of geophone 
client1 = Client('RASPISHAKE')
metadata1 = client1.get_stations(network=netStation1, station=namStation1, location=locStation1, channel=chaStation1, level='resp')
idSEED1 = netStation1 + '.' + namStation1 + '.' + locStation1 + '.' + chaStation1 
latStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['latitude']
lonStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['longitude']
eleStation1 = metadata1[0].get_coordinates(idSEED1, UTCDateTime())['elevation']

#   now get some station meta data, esp. lat & long 
namStation2 = 'WPS'    #   BGS Ladybower 
netStation2 = 'GB'
locStation2 = '00'
chaStation2 = 'HHZ'      #   vertical component of geophone 
client2 = Client('http://eida.bgs.ac.uk')
metadata2 = client2.get_stations(network=netStation2, station=namStation2, location=locStation2, channel=chaStation2, level='resp')
idSEED2 = netStation2 + '.' + namStation2 + '.' + locStation2 + '.' + chaStation2 
latStation2 = metadata2[0].get_coordinates(idSEED2, UTCDateTime())['latitude']
lonStation2 = metadata2[0].get_coordinates(idSEED2, UTCDateTime())['longitude']
eleStation2 = metadata2[0].get_coordinates(idSEED2, UTCDateTime())['elevation']

#   event datetime, e.g. from USGS or BGS
EQ_TIME1 = "2021-12-16T12:48:30"
dtEvent1 = UTCDateTime(EQ_TIME1)

#   get the data
#   plot seismograms  
t1 = dtEvent1 - T_START
t2 = dtEvent1 + T_END

# Download and filter data for station1
st1 = client1.get_waveforms(netStation1, namStation1, locStation1, chaStation1,
                          starttime=t1, endtime=t2, attach_response=True)
#st1.merge(method=0, fill_value='latest')
#st1.detrend(type="demean")
#st1.remove_response()

st1f = st1.copy()
st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
st1.trim(t1, t2)
st1f.trim(t1, t2)

# Download and filter data for station2
st2 = client2.get_waveforms(netStation2, namStation2, locStation2, chaStation2,
                          starttime=t1, endtime=t2, attach_response=True)
st2.merge(method=0, fill_value='latest')
st2.detrend(type="demean")
st2.remove_response()
st2f = st2.copy()
st2f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
st2.trim(t1, t2)
st2f.trim(t1, t2)

# Now plot the waveform data
fig, axs = plt.subplots(2, 1, figsize=(12,8))

axs[0].plot(st1[0].times(reftime=dtEvent1+1), st1[0].data, linewidth=1,
        color="darkred")
axs[0].grid(True)
axs[0].set_xlabel("Time after earthquake (s)")
axs[0].set_title("RaspberryShake {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
    st1f[0].stats.channel, Flow, Fhigh))
axs[0].set_ylabel("Ground velocity (mm/s)")
#axs[0].set_xlim(T_START+1, T_END)
#axs[0].set_xlim(10, 300)
#axs[0].set_ylim(yminf, ymaxf)
#axs[0].set_ylim(-0.003, 0.003)

axs[1].plot(st2f[0].times(reftime=dtEvent1+1), st2[0].data*1000, linewidth=1,
        color="darkblue")
#ymins2f, ymaxs2f = axs[1].get_ylim()
axs[1].grid(True)
axs[1].set_xlabel("Time after earthquake (s)")
axs[1].set_title("BGS station {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    st2f[0].stats.network, st2f[0].stats.station, st2f[0].stats.location,
    st2f[0].stats.channel, Flow, Fhigh))
axs[1].set_ylabel("Ground velocity (mm/s)")
#axs[1].set_xlim(T_START+1, T_END)
#axs[1].set_xlim(10, 300)
#axs[1].set_ylim(yminf, ymaxf)
#axs[1].set_ylim(-0.00025, 0.00025)

plt.tight_layout() 
plt.savefig('plotCompareR6055_v_BGS.png', dpi=300)
