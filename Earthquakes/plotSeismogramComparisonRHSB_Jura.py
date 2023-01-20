#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:58:30 2021

@author: davidhealy
"""
import numpy as np 
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fftshift

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from obspy.geodetics import gps2dist_azimuth

T_START = 0     # length in seconds of data to plot before origin time
T_END = 600    # length in seconds of data to plot after origin time
Flow = 3.
Fhigh = 8.  

#   some standard phases 
PHASES = ["P", "pP", "PP", "PKP"]   # list of phases to compute theoretical times for
MODEL = 'iasp91'  # Velocity model to predict travel-times through
model = TauPyModel(model=MODEL)

#   now get some station meta data, esp. lat & long 
namStation1 = 'R9E71'    #   RHS Bridgewater 
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
namStation2 = 'LBWR'    #   BGS Ladybower 
netStation2 = 'GB'
locStation2 = '00'
chaStation2 = 'HHZ'      #   vertical component of geophone 
client2 = Client('http://eida.bgs.ac.uk')
metadata2 = client2.get_stations(network=netStation2, station=namStation2, location=locStation2, channel=chaStation2, level='resp')
idSEED2 = netStation2 + '.' + namStation2 + '.' + locStation2 + '.' + chaStation2 
latStation2 = metadata2[0].get_coordinates(idSEED2, UTCDateTime())['latitude']
lonStation2 = metadata2[0].get_coordinates(idSEED2, UTCDateTime())['longitude']
eleStation2 = metadata2[0].get_coordinates(idSEED2, UTCDateTime())['elevation']

#   now get some station meta data, esp. lat & long 
namStation3 = 'R72C7'    #   BGS DRUM 
netStation3 = 'AM'
locStation3 = '00'
chaStation3 = 'EHZ'      #   vertical component of geophone 
client3 = Client('RASPISHAKE')
metadata3 = client3.get_stations(network=netStation3, station=namStation3, location=locStation3, channel=chaStation3, level='resp')
idSEED3 = netStation3 + '.' + namStation3 + '.' + locStation3 + '.' + chaStation3 
latStation3 = metadata3[0].get_coordinates(idSEED3, UTCDateTime())['latitude']
lonStation3 = metadata3[0].get_coordinates(idSEED3, UTCDateTime())['longitude']
eleStation3 = metadata3[0].get_coordinates(idSEED3, UTCDateTime())['elevation']

#   event datetime, e.g. from USGS or BGS
EQ_TIME1 = "2021-11-16T01:44:07"
latEvent1 = 56.033 
lonEvent1 = -5.537 
dtEvent1 = UTCDateTime(EQ_TIME1)
depEvent1 = 12 
sEventName1 = "Jura M3.4 " + EQ_TIME1 

#   work out epicentral distance from station, degrees and km 
distDeg1 = locations2degrees(latEvent1, lonEvent1, latStation1, lonStation1)
distKm1, _, _ = gps2dist_azimuth(latStation1, lonStation1, latEvent1, lonEvent1)   
print("Epicentral distance 1 (degrees): %3.1f" % distDeg1) 

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

# Download and filter data for station3
st3 = client3.get_waveforms(netStation3, namStation3, locStation3, chaStation3,
                          starttime=t1, endtime=t2, attach_response=True)
st3.merge(method=0, fill_value='latest')
st3.detrend(type="demean")
st3.remove_response()
st3f = st3.copy()
st3f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
st3.trim(t1, t2)
st3f.trim(t1, t2)

ymin = np.min([st1[0].data, st2[0].data, st3[0].data])*1000.*1.1
yminf = np.min([st1f[0].data, st2f[0].data, st3f[0].data])*1000.*1.1
ymax = np.max([st1[0].data, st2[0].data, st3[0].data])*1000.*1.1 
ymaxf = np.max([st1f[0].data, st2f[0].data, st3[0].data])*1000.*1.1

# Now plot the waveform data
fig, axs = plt.subplots(2, 1, figsize=(12,8))

axs[0].plot(st1f[0].times(reftime=dtEvent1+1), st1f[0].data*1000, linewidth=1,
        color="darkred")
axs[0].grid(True)
axs[0].set_xlabel("Time after earthquake (s)")
axs[0].set_title("{:} - RaspberryShake {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    sEventName1, st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
    st1f[0].stats.channel, Flow, Fhigh))
axs[0].set_ylabel("Ground velocity (mm/s)")
#axs[0].set_xlim(T_START+1, T_END)
axs[0].set_xlim(10, 300)
#axs[0].set_ylim(yminf, ymaxf)
axs[0].set_ylim(-0.003, 0.003)

axs[1].plot(st2f[0].times(reftime=dtEvent1+1), st2f[0].data*1000, linewidth=1,
        color="darkblue")
#ymins2f, ymaxs2f = axs[1].get_ylim()
axs[1].grid(True)
axs[1].set_xlabel("Time after earthquake (s)")
axs[1].set_title("{:} - BGS station {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    sEventName1, st2f[0].stats.network, st2f[0].stats.station, st2f[0].stats.location,
    st2f[0].stats.channel, Flow, Fhigh))
axs[1].set_ylabel("Ground velocity (mm/s)")
#axs[1].set_xlim(T_START+1, T_END)
axs[1].set_xlim(10, 300)
#axs[1].set_ylim(yminf, ymaxf)
axs[1].set_ylim(-0.00025, 0.00025)

# axs[2].plot(st3f[0].times(reftime=dtEvent1+1), st3f[0].data*1000, linewidth=1,
#         color="darkblue")
# #ymins2f, ymaxs2f = axs[1].get_ylim()
# axs[2].grid(True)
# axs[2].set_xlabel("Time after earthquake (s)")
# axs[2].set_title("{:} - Raspberry Shake {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
#     sEventName1, st3f[0].stats.network, st3f[0].stats.station, st3f[0].stats.location,
#     st3f[0].stats.channel, Flow, Fhigh))
# axs[2].set_ylabel("Ground velocity (mm/s)")
# #axs[2].set_xlim(T_START+1, T_END)
# axs[2].set_xlim(250, 1000)
# #axs[2].set_ylim(yminf, ymaxf)
# #axs[2].set_ylim(-0.025, 0.025)

plt.tight_layout() 
plt.savefig('plotCompareRHSvBGS_Jura.png', dpi=300)
