#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:58:30 2021

@author: davidhealy
"""
import matplotlib.pyplot as plt

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*10     # length in seconds of data to plot after origin time
Flow = 0.7
Fhigh = 2.   
PHASES = ["P"]   # list of phases to compute theoretical times for
#PHASES = []   # list of phases to compute theoretical times for
MODEL = "iasp91"  # Velocity model to predict travel-times through
model = TauPyModel(model=MODEL)

#   now get some station meta data, esp. lat & long 
namStationList = ['R4DD7', 'R3FEA', 'R6C8A', 'RAEED', 'R36CB', 'RE28A', 'S53C6']
#namStationList = ['RE28A']
#namStationList = ['R6055', 'R72C7', 'RE28A']
labelStationList = ['Firs', 'Chorlton', 'MMU', 'UoM', 'Cheadle Hulme', 'Menai', 'Macclesfield']
#labelStationList = ['Menai Bridge']
#labelStationList = ['Amlwch', 'Aberdeen', 'Menai']
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
#chaStationList = ['EHZ', 'EHE', 'EHN']      #   vertical component of geophone 
client1 = Client('RASPISHAKE')

#   event datetime, e.g. from USGS or BGS
# EQ_TIME = "2022-03-16T14:36:33"
# EQ_NAME = "M7.3 Honshu 2022-03-16"
# EQ_TIME = "2022-03-16T23:15:45"
# EQ_NAME = "M5.8 Iran 2022-03-16"
# EQ_TIME = "2022-03-22T17:41:38"
# EQ_NAME = "M6.7 Taiwan 2022-03-22"
# EQ_TIME = "2022-03-30T20:56:58"
# EQ_NAME = "M6.9 New Caledonia 2022-03-30"
# EQ_TIME = "2022-05-29T20:40:37.9"
# EQ_NAME = "M2.3 Sale, GM 2022-05-29"
#EQ_TIME = "2022-05-30T07:59:17.3"
#EQ_NAME = "M2.1 Arran 2022-05-30"
# EQ_TIME = "2022-06-01T01:51:22.9"
# EQ_NAME = "M2.8 Laneshawbridge, Lancs 2022-06-01"
#EQ_TIME = "2022-06-08T00:55:47"
#EQ_NAME = "M6.5 Tarauaca, Brazil 2022-06-08"
#EQ_TIME = "2022-06-15T06:06:03"
#EQ_NAME = "M5.5 Southern Iran, 2022-06-15"
EQ_TIME = "2022-11-09T06:07:27"
EQ_NAME = "M5.6 Marotta, Italy, 2022-11-09"
EVT_LAT = 43.930
EVT_LON = 13.310
EVT_Z = 10. 
dtEvent = UTCDateTime(EQ_TIME)
t1start = dtEvent + T_START
t1end = dtEvent + T_END

# Now plot the waveform data
fig, axs = plt.subplots(len(namStationList), 1, figsize=(10,20))

i = 0 
for namStation in namStationList:

    metadata = client1.get_stations(network=netStation, station=namStation, location=locStation, channel=chaStation, level='resp')
    idSEED = netStation + '.' + namStation + '.' + locStation + '.' + chaStation 
    latStation = metadata[0].get_coordinates(idSEED, UTCDateTime())['latitude']
    lonStation = metadata[0].get_coordinates(idSEED, UTCDateTime())['longitude']
    eleStation = metadata[0].get_coordinates(idSEED, UTCDateTime())['elevation']

    dist = locations2degrees(EVT_LAT, EVT_LON, latStation, lonStation)
    print(dist)
    
    # Download and filter data for station1
    st1 = client1.get_waveforms(netStation, namStation, locStation, chaStation,
                          starttime=t1start, endtime=t1end, attach_response=True)
    st1.merge(method=0, fill_value='latest')
    st1.detrend(type="demean")
    st1.remove_response()
    st1f = st1.copy()
    st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
    st1.trim(t1start+20, t1end)
    st1f.trim(t1start+20, t1end)

    axs[i].plot(st1f[0].times(reftime=dtEvent), st1f[0].data*1000, linewidth=1,
            color="darkred")
    ymin, ymax = axs[i].get_ylim()
    
    # # Now plot the theoretical arrival times
    for phase in PHASES:
        
        phase = [phase]
        tt = model.get_travel_times(source_depth_in_km=EVT_Z,
                                    distance_in_degree=dist,
                                    phase_list=phase)
        axs[i].vlines(tt[0].time, ymin, ymax, color="blue",
                  linewidth=1.2, zorder=3, linestyle="--", alpha=0.5)
        axs[i].text(tt[0].time*1.005, ymax*1.05, phase[0], fontsize=12,
                horizontalalignment="left", verticalalignment="top")

    axs[i].grid(True)
    axs[i].set_xlabel("Time after earthquake (s)")
    axs[i].set_title("{:} {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz - {:}".format(
        EQ_NAME, st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
        st1f[0].stats.channel, Flow, Fhigh, labelStationList[i]))
    axs[i].set_ylabel("Ground velocity (mm/s)")
#     if i >= 4: 
#         axs[i].set_xlim(10, 80)
# #        axs[i].set_ylim(-0.0025, 0.0025)
#     else: 
#         axs[i].set_xlim(0, 50)
#         axs[i].set_ylim(-0.02, 0.02)
#    axs[i].autoscale(enable=True, axis='y', tight=True)
#    axs[i].set_ylim(ymin, ymax)
#    axs[i].set_xlim(400, 1000)

    i = i + 1

plt.tight_layout()
fn = EQ_NAME + ".png" 
plt.savefig(fn, dpi=300)
