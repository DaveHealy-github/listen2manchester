#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 11:02:25 2021

script to read the USGS GeoJSON feed and plot RShake data for each quake

@author: davidhealy
"""

import datetime 
import urllib
import json 
import numpy as np 
import matplotlib.pyplot as plt

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.taup import plot_travel_times
from obspy.geodetics.base import locations2degrees
from obspy.geodetics import gps2dist_azimuth

#   some initialisation... 
PHASES = ["P", "pP", "PP", "PcP", "PKiKP", "sP", "S"]   # list of phases to compute theoretical times for
MODEL = 'iasp91'  # Velocity model to predict travel-times through
model = TauPyModel(model=MODEL)

NMINUTES = 30. 
T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*NMINUTES    # length in seconds of data to plot after origin time
Flow1 = 0.7
Flow2 = 0.1 
Fhigh1 = 2. 
Fhigh2 = 0.8 
magThreshold = 5.5

#   now get some station meta data, esp. lat & long 
#   now get some station meta data, esp. lat & long 
namStationList = ['R4DD7', 'R3FEA', 'R6C8A', 'RAEED', 'R36CB', 'RE28A', 'R6055']
labelStationList = ['Firs', 'Chorlton', 'MMU', 'UoM', 'Cheadle Hulme', 'Menai', 'Amlwch']
#namStation = 'R72C7'    #   in my office at Aberdeen 
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

#   get the data for the last 30 days from USGS; GeoJSON summary feed 
url = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/significant_month.geojson"
uh = urllib.request.urlopen(url)
datas = uh.read()
data = json.loads(datas)

#   check the last time we ran this script; 
#       only process 'new' quakes since the last run... 
with open("USGS_JSON_log.txt", mode='r') as file:
    sdtLog = file.read().split('\n')
    dtLog = datetime.datetime.strptime(sdtLog[0], '%Y-%m-%d %H:%M:%S.%f')
print(dtLog)
     
#   loop through the data, stripping out lat, long, depth and name 
nEvents = data['metadata']['count']
print('Processing %i significant events for last month (USGS)...' % nEvents)

for namStation in namStationList:

    print('Raspberry Shake: ', namStation)
    
    metadata = client.get_stations(network='AM', station=namStation, location='00', channel='*HZ', level='resp')
    idSEED = netStation + '.' + namStation + '.' + locStation + '.' + chaStation 
    latStation = metadata[0].get_coordinates(idSEED, UTCDateTime())['latitude']
    lonStation = metadata[0].get_coordinates(idSEED, UTCDateTime())['longitude']
    eleStation = metadata[0].get_coordinates(idSEED, UTCDateTime())['elevation']

    for iEvent in range(0, nEvents):
        
        #   get the details for this event...  
        latEvent = data['features'][iEvent]['geometry']['coordinates'][1]
        lonEvent = data['features'][iEvent]['geometry']['coordinates'][0]
        depEvent = data['features'][iEvent]['geometry']['coordinates'][2]
        magEvent = data['features'][iEvent]['properties']['mag']
        namEvent = data['features'][iEvent]['properties']['place']
        dtEvent = UTCDateTime(data['features'][iEvent]['properties']['time']/1000.)
    
        #   Firs station wasn't online until June 10 2022 
        if namStation == 'R4DD7' and dtEvent < '2022-06-10T00:00:00.0':
            continue 
        
        #   don't process events already plotted 
        if UTCDateTime(dtLog) > dtEvent:
            continue 
        
        #   arbitrary threshold for seeing better waveforms on RShakes...
        if magEvent >= magThreshold:
            
            print('*** Event magnitude >= %1.1f' % magThreshold)
            print(namEvent, magEvent, depEvent, dtEvent)
            
            #   work out epicentral distance from station, degrees and km 
            distDeg = locations2degrees(latEvent, lonEvent, latStation, lonStation)
            distm, _, _ = gps2dist_azimuth(latStation, lonStation, latEvent, lonEvent)
            distKm = distm / 1000.  
            print('Epicentral distance (degrees): %3.2f\n' % distDeg) 
            
            fig, ax = plt.subplots(figsize=(5, 5))
            ax = plot_travel_times(source_depth=depEvent, 
                                   phase_list=PHASES,
                                   ax=ax, fig=fig, verbose=True, show=False)
            ax = fig.gca()
            ax.plot([distDeg, distDeg], [0., NMINUTES], ':r')
            ax.grid(True)
            ax.autoscale(enable=True, axis='both', tight=True)
#            plt.show()

            fn = namStation + '/' + namStation + ' - travel times - M' + str(magEvent) + ' - ' + namEvent + ".png"
            fig.savefig(fn, dpi=300)
            
            fig, ax = plt.subplots(figsize=(5, 5), subplot_kw={'projection': 'polar'})
            arrivals = model.get_ray_paths(source_depth_in_km=depEvent, 
                                            distance_in_degree=distDeg, 
                                            phase_list=PHASES)
            arrivals.plot_rays(ax=ax, fig=fig, legend=True, show=False)            
            ax = fig.gca()
            ax.grid(True)
#            plt.show()            

            fn = namStation + '/' + namStation + ' - ray paths - M' + str(magEvent) + ' - ' + namEvent + ".png" 
            fig.savefig(fn, dpi=300)
            
            #   plot seismograms - raw, plus 2 filtered 
            t1 = dtEvent - T_START
            t2 = dtEvent + T_END
            # Download and filfter data
            st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                                      starttime=t1, endtime=t2, attach_response=True)
            st.merge(method=0, fill_value='latest')
            st.detrend(type="demean")
            st.remove_response()
            st.filter("bandpass", freqmin=Flow1, freqmax=Fhigh1, corners=4)
            
            # first-pass through theoretical arrival times, save arrivals 
            ttSave = [] 
            for phase in PHASES:
                phase = [phase]
                tt = model.get_travel_times(source_depth_in_km=depEvent,
                                            distance_in_degree=distDeg,
                                            phase_list=phase)
                if len(tt) > 0:
                    ttSave.append(tt[0].time)
                    
            # Set-up figure
            fig, ax = plt.subplots(figsize=(12, 4))
            sMainTitle = ("M %2.1f - %s - %s\nEpicentral distance = %3.1f (degrees), %6.1f (km)" % (magEvent, namEvent, dtEvent.ctime(), distDeg, distKm))
            fig.suptitle(sMainTitle)
    
            #   use the phase arrivals to window the waveform data...
            tPhaseFirst = np.min(ttSave)*.8 
            tPhaseLast = np.max(ttSave)*1.1
            mask = ((st[0].times(reftime=dtEvent) >= tPhaseFirst) & 
                    (st[0].times(reftime=dtEvent) <= tPhaseLast))
    
            #   now plot the waveform data
            ax.plot(st[0].times(reftime=dtEvent)[mask], 
                    st[0].data[mask]*1000, 
                    linewidth=1,
                    color="darkred")
    
            #   now get y limits for this plotted data window...         
            ymin, ymax = ax.get_ylim()
            ylim = np.max([abs(ymin), abs(ymax)])
            
            #   plot the theoretical arrival times
            for phase in PHASES:
                phase = [phase]
                tt = model.get_travel_times(source_depth_in_km=depEvent,
                                            distance_in_degree=distDeg,
                                            phase_list=phase)
                if len(tt) > 0:
                    ax.vlines(tt[0].time, -ylim, ylim, color="blue",
                          linewidth=1.2, zorder=3, linestyle="--", alpha=0.5)
                    ax.text(tt[0].time*1.005, ylim*.95, phase[0], fontsize=12,
                        horizontalalignment="left", verticalalignment="top")
                    
            ax.grid(True)        
            ax.set_xlabel("Time after earthquake (s)")
            ax.set_title("{:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
                st[0].stats.network, st[0].stats.station, st[0].stats.location,
                st[0].stats.channel, Flow1, Fhigh1))
            ax.set_ylabel("Ground velocity (mm/s)")
            ax.autoscale(enable=True, axis='both', tight=True)
            
            plt.tight_layout()
            fn = namStation + '/' + namStation + ' - M' + str(magEvent) + ' - ' + namEvent + ".png" 
            plt.savefig(fn, dpi=300)
            
with open("USGS_JSON_log.txt", mode='w') as file:
    file.write('%s' % datetime.datetime.now())
