#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:58:30 2021

@author: davidhealy
"""
import matplotlib.pyplot as plt
import numpy as np 
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.taup import plot_travel_times
from obspy.geodetics.base import locations2degrees
from obspy.geodetics import gps2dist_azimuth

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*40    # length in seconds of data to plot after origin time
Flow = 0.5
Fhigh = 2.   

#   some initialisation... 
PHASES = ["P", "pP", "S"]   # list of phases to compute theoretical times for
MODEL = 'iasp91'  # Velocity model to predict travel-times through
model = TauPyModel(model=MODEL)

#   now get some station meta data, esp. lat & long 
namStation1 = 'RCFCD'    #   Runshaw College, Leyland nr Preston    
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
EQ_TIME1 = "2023-10-16T11:35:31"
dtEvent1 = UTCDateTime(EQ_TIME1)
t1start = dtEvent1 - T_START
t1end = dtEvent1 + T_END
latEvent = 52.445 
lonEvent = 176.850 
depEvent = 187.4    #   in km 

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

#   work out epicentral distance from station, degrees and km 
distDeg = locations2degrees(latEvent, lonEvent, latStation1, lonStation1)

# Now plot the waveform data
fig, axs = plt.subplots(1, 1, figsize=(12,5))

axs.plot(st1f[0].times(reftime=dtEvent1+1), st1f[0].data*1000, linewidth=1,
        color="darkred")

#   now get y limits for this plotted data window...         
ymin, ymax = axs.get_ylim()
ylim = np.max([abs(ymin), abs(ymax)])

#   plot the theoretical arrival times
for phase in PHASES:
    phase = [phase]
    tt = model.get_travel_times(source_depth_in_km=depEvent,
                                distance_in_degree=distDeg,
                                phase_list=phase)
    if len(tt) > 0:
        axs.vlines(tt[0].time, -ylim, ylim, color="blue",
              linewidth=1.2, zorder=3, linestyle="--", alpha=0.5)
        axs.text(tt[0].time*1.005, ylim*.95, phase[0], fontsize=12,
            horizontalalignment="left", verticalalignment="top")

axs.grid(True)
axs.set_xlabel("Time after earthquake (s)")
axs.set_title("M6.4 Aleutian Islands 2023-10-16 11:35 {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
    st1f[0].stats.channel, Flow, Fhigh))
axs.set_ylabel("Ground velocity (mm/s)")
axs.set_xlim(400, 1500)

plt.tight_layout() 
plt.savefig('plotRCFCD Runshaw College - M6.4 Aleutian Islands 2023-10-16.png', dpi=300)
