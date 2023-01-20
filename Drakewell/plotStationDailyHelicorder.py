#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 12:46:19 2022

@author: davidhealy
"""
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

#   now get some station meta data, esp. lat & long 
namStations = ['R72C7',     #   in my office at Aberdeen 
               'R6055',     #   GeoMon Amlwch 
               'RE28A',     #   Menai Bridge
               'R3FEA']     #   Chorlton 
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

START_TIME = "2022-01-12T00:00:00"
dtStart = UTCDateTime(START_TIME)
T_END = 60*60*24

#   loop through each station in the list 
iStation = 0 
for namStation in namStations:       
    
    #   plot seismograms - raw, plus 2 filtered 
    t1 = dtStart
    t2 = dtStart + T_END
    
    # Download and filfter data
    st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                              starttime=t1, endtime=t2, attach_response=True)
    st.merge(method=0, fill_value='latest')
#    st.detrend(type="demean")
    st.decimate(1)
#    st.remove_response()
#    st.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)

    st.plot(type="dayplot", interval=60, linewidth=0.3, 
            right_vertical_labels=False,
#            vertical_scaling_range=1e6, 
            one_tick_per_line=True,
            color=['k', 'r', 'b', 'g'], 
            show_y_UTC_label=False)     