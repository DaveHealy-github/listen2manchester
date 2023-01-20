#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 14:17:26 2022

@author: davidhealy
"""
import matplotlib.pyplot as plt

from obspy import read
from obspy.io.xseed import Parser
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*60*24     # length in seconds of data to plot after origin time
Fhigh = 1.0 

#   now get some station meta data, esp. lat & long 
namStationList = ['R4DD7', 'R3FEA', 'R6C8A', 'RAEED', 'R36CB', 'LBWR']
labelStationList = ['Firs', 'Chorlton', 'MMU', 'UoM', 'Cheadle Hulme', 'Ladybower']
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client1 = Client('RASPISHAKE')

EQ_TIME = "2022-06-26T00:00:00"
dtEvent = UTCDateTime(EQ_TIME)
t1start = dtEvent + T_START
t1end = dtEvent + T_END

i = 0 
for namStation in namStationList:

    if "LBWR" in namStation:
        client1=Client('http://eida.bgs.ac.uk')
        netStation = 'GB'
        locStation = '00' 
        chaStation = 'HHZ'
        
    # Download and filter data for station1
    resp = client1.get_stations(dtEvent, network=netStation, station=namStation, location=locStation, channel=chaStation, level="response")
    print(resp)
    
    st1 = client1.get_waveforms(netStation, namStation, locStation, chaStation,
                          starttime=t1start, endtime=t1end, attach_response=True)
    
    ppsd = PPSD(st1[0].stats, metadata=resp, 
                    ppsd_length=3600, overlap=0.5)
    ppsd.add(st1)    
    ppsd.plot(cmap=pqlx, 
              xaxis_frequency=True, 
              period_lim=[0.5, 50.], 
              show_mean=True, 
              show_coverage=False,
              show_histogram=False,
              grid=True)