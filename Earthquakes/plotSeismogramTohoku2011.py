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
tend=t_eq+60*120 #stop downloading data 1500s after the earthquake

st2=client.get_waveforms("GB","CWF", "*", "BHZ", tstart, tend, attach_response=True) 

st2.merge(method=0, fill_value='latest')
st2.detrend(type="spline", order=2, dspline=1000)
st2.remove_response()
st2f = st2.copy()
st2f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
st2.trim(tstart, tend)
st2f.trim(tstart, tend)

f, Pxx_den = signal.periodogram(st2[0].data*1000., 1)
ff, Pxx_denf = signal.periodogram(st2f[0].data*1000., 1)

#st2.write("Tohoku.txt", format='SLIST')

ymin = np.min(st2[0].data)*1000.*1.1
yminf = np.min(st2f[0].data)*1000.*1.1
ymax = np.max(st2[0].data)*1000.*1.1 
ymaxf = np.max(st2f[0].data)*1000.*1.1

# Now plot the waveform data
fig, axs = plt.subplots(3, 1, figsize=(12,12))

axs[0].plot(st2[0].times(reftime=tstart+1), st2[0].data*1000, linewidth=1,
        color="darkblue")
#ymins2f, ymaxs2f = axs[1].get_ylim()
axs[0].grid(True)
axs[0].set_xlabel("Time after earthquake (s)")
# axs.set_title("{:} - BGS station {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
#     sEventName1, st2f[0].stats.network, st2f[0].stats.station, st2f[0].stats.location,
#     st2f[0].stats.channel, Flow, Fhigh))
axs[0].set_ylabel("Ground velocity (mm/s)")
#axs[1].set_xlim(T_START+1, T_END)
#axs.set_xlim(1050, 1400)
#axs[1].set_ylim(yminf, ymaxf)
#axs.set_ylim(-0.025, 0.025)
# # Now plot the theoretical arrival times
# for phase in PHASES:
#     phase = [phase]
#     tt = model.get_travel_times(source_depth_in_km=depEvent1,
#                                 distance_in_degree=distDeg1,
#                                 phase_list=phase)
#     if len(tt) > 0:
#         axs[1].vlines(tt[0].time, yminf, ymaxf, color="red",
#               linewidth=1.2, zorder=3, linestyle="--", alpha=0.5)
#         axs[1].text(tt[0].time*1.005, 0.02, phase[0], fontsize=12,
#             horizontalalignment="left", verticalalignment="top")

axs[1].plot(st2f[0].times(reftime=tstart+1), st2f[0].data*1000, linewidth=1,
        color="darkblue")
#ymins2f, ymaxs2f = axs[1].get_ylim()
axs[1].grid(True)
axs[1].set_xlabel("Time after earthquake (s)")
# axs.set_title("{:} - BGS station {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
#     sEventName1, st2f[0].stats.network, st2f[0].stats.station, st2f[0].stats.location,
#     st2f[0].stats.channel, Flow, Fhigh))
axs[1].set_ylabel("Ground velocity (mm/s)")
#axs[1].set_xlim(T_START+1, T_END)
#axs.set_xlim(1050, 1400)
#axs[1].set_ylim(yminf, ymaxf)
#axs.set_ylim(-0.025, 0.025)
# # Now plot the theoretical arrival times
# for phase in PHASES:
#     phase = [phase]
#     tt = model.get_travel_times(source_depth_in_km=depEvent1,
#                                 distance_in_degree=distDeg1,
#                                 phase_list=phase)
#     if len(tt) > 0:
#         axs[1].vlines(tt[0].time, yminf, ymaxf, color="red",
#               linewidth=1.2, zorder=3, linestyle="--", alpha=0.5)
#         axs[1].text(tt[0].time*1.005, 0.02, phase[0], fontsize=12,
#             horizontalalignment="left", verticalalignment="top")

axs[2].loglog(f, Pxx_den, '-b')
axs[2].loglog(ff, Pxx_denf, '-r')
axs[2].grid(True)
axs[2].set_xlabel("Frequency, Hz")
axs[2].set_ylabel("Power spectral density ((mm/s)^2/Hz)")
axs[2].set_xlim(1e-5, 1e0)
axs[2].set_ylim(1e-23, 1e3)

plt.tight_layout() 
plt.savefig('plotTohoku.png', dpi=300)
