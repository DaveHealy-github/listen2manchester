#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:58:30 2021

@author: davidhealy
"""
import numpy as np 
import matplotlib.pyplot as plt
from scipy import signal

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*120    # length in seconds of data to plot after origin time
Flow = 0.5
Fhigh = 1.5  

#   some standard phases 
PHASES = ["P", "pP", "PP", "PKP"]   # list of phases to compute theoretical times for
MODEL = 'iasp91'  # Velocity model to predict travel-times through
model = TauPyModel(model=MODEL)

# #   now get some station meta data, esp. lat & long 
# namStation2 = 'CWF'    #   BGS CWF 
# netStation2 = 'GB'
# locStation2 = '??'
# chaStation2 = 'BHZ'      #   vertical component 
# client2 = Client('http://eida.bgs.ac.uk')
# metadata2 = client2.get_stations(network=netStation2, station=namStation2, location=locStation2, channel=chaStation2, level='resp')
# idSEED2 = netStation2 + '.' + namStation2 + '.' + locStation2 + '.' + chaStation2 
# latStation2 = metadata2[0].get_coordinates(idSEED2, UTCDateTime())['latitude']
# lonStation2 = metadata2[0].get_coordinates(idSEED2, UTCDateTime())['longitude']
# eleStation2 = metadata2[0].get_coordinates(idSEED2, UTCDateTime())['elevation']

# #   event - Tohoku, M9.1
# EQ_TIME1 = "2011-03-11T05:46:24"
# latEvent1 = 38.297 
# lonEvent1 = 142.373 
# dtEvent1 = UTCDateTime(EQ_TIME1)
# depEvent1 = 19.7 
# sEventName1 = "Tohoku M9.1 " + EQ_TIME1 

# #   work out epicentral distance from station, degrees and km 
# distDeg2 = locations2degrees(latEvent1, lonEvent1, latStation2, lonStation2)
# distKm2, _, _ = gps2dist_azimuth(latStation2, lonStation2, latEvent1, lonEvent1)   
# print("Epicentral distance (degrees): %3.1f" % distDeg2) 

# #   get the data
# #   plot seismograms  
# t1 = dtEvent1 - T_START
# t2 = dtEvent1 + T_END

# # Download and filter data for station2
# st2 = client2.get_waveforms(netStation2, namStation2, locStation2, chaStation2,
#                           starttime=t1, endtime=t2, attach_response=True)

client=Client("IRIS") # set IRIS as the place you're going to download data from

t_eq=UTCDateTime("2011-03-11T05:46:24.000") # set a date and time for the earthquake
tstart=t_eq #start downloading from 300s after the earthquake
tend=t_eq+60*32 #stop downloading data 1500s after the earthquake

st2=client.get_waveforms("GB","CWF", "*", "BHZ", tstart, tend, attach_response=True) 

st2.merge(method=0, fill_value='latest')
st2.detrend(type="spline", order=2, dspline=1000)
st2.remove_response()
st2f = st2.copy()
st2f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
st2.trim(tstart, tend)
st2f.trim(tstart, tend)

# Now plot the waveform data
fig, axs = plt.subplots(1, 1, figsize=(12,5))
axs.plot(st2f[0].times(reftime=tstart+1), st2f[0].data*1000, linewidth=1.5,
        color=(0.85,0.16,0.11))
axs.set_xlabel("Time after earthquake (s)")
axs.set_ylabel("Ground velocity (mm/s)")
axs.set_xlim(600, 1600)
plt.tight_layout() 
plt.savefig('plotTohokuLogoRed.png', dpi=600)

# Now plot the waveform data
fig, axs = plt.subplots(1, 1, figsize=(12,5))
axs.plot(st2f[0].times(reftime=tstart+1), st2f[0].data*1000, linewidth=1.5,
        color=(0.42,0.67,0.87))
axs.set_xlabel("Time after earthquake (s)")
axs.set_ylabel("Ground velocity (mm/s)")
axs.set_xlim(400, 1400)
plt.tight_layout() 
plt.savefig('plotTohokuLogoBlue.png', dpi=600)
