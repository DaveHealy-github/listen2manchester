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
T_END = 60*60    # length in seconds of data to plot after origin time
Flow = 4.
Fhigh = 14.  

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

#client=Client("IRIS") # set IRIS as the place you're going to download data from
client=Client('http://eida.bgs.ac.uk')
station = "HLM1"

t_eq=UTCDateTime("2022-05-30T14:36:57.7") # set a date and time for the earthquake
tstart=t_eq #start downloading from 300s after the earthquake
tend=t_eq+60*3 #stop downloading data 1500s after the earthquake

stz=client.get_waveforms("GB", station, "00", "HHZ", tstart, tend, attach_response=True) 
stz.merge(method=0, fill_value='latest')
stz.detrend(type="spline", order=2, dspline=1000)
stz.detrend(type="demean")
stz.remove_response()

stn=client.get_waveforms("GB", station, "00", "HHN", tstart, tend, attach_response=True) 
stn.merge(method=0, fill_value='latest')
stn.detrend(type="spline", order=2, dspline=1000)
stn.detrend(type="demean")
stn.remove_response()

ste=client.get_waveforms("GB", station, "00", "HHE", tstart, tend, attach_response=True) 
ste.merge(method=0, fill_value='latest')
ste.detrend(type="spline", order=2, dspline=1000)
ste.detrend(type="demean")
ste.remove_response()

# Now plot the waveform data
fig, axs = plt.subplots(3, 1, figsize=(12,6))

axs[0].plot(stz[0].times(reftime=tstart+1), stz[0].data*1000, linewidth=1,
        color="darkblue")
axs[0].grid(True)
axs[0].set_xlabel("Time after earthquake (s)")
axs[0].set_ylabel("BHZ Ground velocity (mm/s)")

axs[1].plot(stn[0].times(reftime=tstart+1), stn[0].data*1000, linewidth=1,
        color="darkblue")
axs[1].grid(True)
axs[1].set_xlabel("Time after earthquake (s)")
axs[1].set_ylabel("BHN Ground velocity (mm/s)")

axs[2].plot(ste[0].times(reftime=tstart+1), ste[0].data*1000, linewidth=1,
        color="darkblue")
axs[2].grid(True)
axs[2].set_xlabel("Time after earthquake (s)")
axs[2].set_ylabel("BHE Ground velocity (mm/s)")

plt.tight_layout() 
plt.savefig('plotBGSstation.png', dpi=300)
