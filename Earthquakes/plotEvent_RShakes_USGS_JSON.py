#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 11:02:25 2021

script to read the USGS GeoJSON feed and plot RShake data for each quake

@author: davidhealy
"""

import urllib
import json 
import matplotlib.pyplot as plt

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from obspy.geodetics import gps2dist_azimuth

#   some initialisation... 
T_START = 0      # length in seconds of data before origin time
T_END = 60*40    # length in seconds of data after origin time
Flow1 = 0.7
Fhigh1 = 2. 
magThreshold = 6.0

PHASES = ["P", "PKiKP", "PP"]   # list of phases to compute theoretical times for
MODEL = 'iasp91'  # Velocity model to predict travel-times through
model = TauPyModel(model=MODEL)

#   now get some station meta data, esp. lat & long 
namStations = ['R72C7',     #   in my office at Aberdeen 
               'R6055',     #   GeoMon Amlwch 
               'RE28A',     #   Menai Bridge
               'R3FEA']     #   Chorlton 
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

#   get the data for the last 30 days from USGS; GeoJSON summary feed 
url = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/significant_month.geojson"
uh = urllib.request.urlopen(url)
datas = uh.read()
data = json.loads(datas)

#   loop through the data, stripping out lat, long, depth and name 
nEvents = data['metadata']['count']
iPlotCount = len(namStations)
        
print('Processing %i significant events for last month (USGS)...' % nEvents)
for iEvent in range(0, nEvents):
    
    #   get the details for this event...  
    latEvent = data['features'][iEvent]['geometry']['coordinates'][1]
    lonEvent = data['features'][iEvent]['geometry']['coordinates'][0]
    depEvent = data['features'][iEvent]['geometry']['coordinates'][2]
    magEvent = data['features'][iEvent]['properties']['mag']
    namEvent = data['features'][iEvent]['properties']['place']
    dtEvent = UTCDateTime(data['features'][iEvent]['properties']['time']/1000.)

    #   check arbitrary threshold 
    #   and only quakes since last daily run 
    if magEvent >= magThreshold and dtEvent > (UTCDateTime.now() - 60*60*24) :
        
        print('Event magnitude >= %1.1f' % magThreshold)
        print(namEvent, magEvent, depEvent, dtEvent)

        #   set-up figure
        fig, axs = plt.subplots(iPlotCount, 1, 
                                figsize=(10, 3*iPlotCount)) 
#                                sharey=True)
        plt.suptitle("M" + str(magEvent) + ' - ' + namEvent + ' - ' + dtEvent.ctime())
        
        #   loop through each station in the list 
        iStation = 0 
        for namStation in namStations:       
            
            #   work out epicentral distance from station, degrees and km 
            metadata = client.get_stations(network='AM', station=namStation, location='00', channel='*HZ', level='resp')
            idSEED = netStation + '.' + namStation + '.' + locStation + '.' + chaStation 
            latStation = metadata[0].get_coordinates(idSEED, UTCDateTime())['latitude']
            lonStation = metadata[0].get_coordinates(idSEED, UTCDateTime())['longitude']
            eleStation = metadata[0].get_coordinates(idSEED, UTCDateTime())['elevation']
            distDeg = locations2degrees(latEvent, lonEvent, latStation, lonStation)
            distKm, _, _ = gps2dist_azimuth(latStation, lonStation, latEvent, lonEvent)   
            print("Epicentral distance (degrees):", distDeg) 
            
            #   plot seismograms - raw, plus 2 filtered 
            t1 = dtEvent
            t2 = dtEvent + T_END
            
            # Download and filfter data
            st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                                      starttime=t1, endtime=t2, attach_response=True)
            st.merge(method=0, fill_value='latest')
            st.detrend(type="demean")
            st.remove_response()
            st.filter("bandpass", freqmin=Flow1, freqmax=Fhigh1, corners=4)
            
            t1Crop = dtEvent + (distDeg * 6.) * .85 
            st.trim(t1Crop, t2)
            
            axs[iStation].plot(st[0].times(reftime=dtEvent), st[0].data*1000, linewidth=1,
                    color="darkred")
            ymin, ymax = axs[iStation].get_ylim()
            axs[iStation].grid(True)
            axs[iStation].set_xlabel("Time after earthquake (s)")
            axs[iStation].set_title("{:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
                st[0].stats.network, st[0].stats.station, st[0].stats.location,
                st[0].stats.channel, Flow1, Fhigh1))
            axs[iStation].set_ylabel("Ground velocity (mm/s)")

            # Now plot the theoretical arrival times
            for phase in PHASES:
                phase = [phase]
                tt = model.get_travel_times(source_depth_in_km=depEvent,
                                            distance_in_degree=distDeg,
                                            phase_list=phase)
                if len(tt) > 0:
                    axs[iStation].vlines(tt[0].time, ymin, ymax, color="blue",
                          linewidth=1.2, zorder=3, linestyle="-", alpha=0.5)
                    axs[iStation].text(tt[0].time*1.005, ymax, phase[0], 
                                       fontsize=10, color="blue",
                        horizontalalignment="left", verticalalignment="top")

            iStation += 1 
            
        fnPlotThis = 'plotEvent_RShakes_USGS_JSON' + '_' + namEvent + '.png'
        plt.tight_layout() 
        plt.savefig(fnPlotThis, dpi=150)
