#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:58:30 2021

@author: davidhealy
"""
#import matplotlib.pyplot as plt

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*60*24    # length in seconds of data to plot after origin time
Flow = 0.5
Fhigh = 40.   

#   now get some station meta data, esp. lat & long 
#namStation1 = 'R3FEA'    #   Oswald Road, Chorlton  
namStation1 = 'R4DD7'    #   Firs, Fallowfield  
netStation1 = 'AM'
locStation1 = '00'
chaStation1 = 'EHZ'      #   vertical component of geophone 
client1 = Client('RASPISHAKE')
idSEED1 = netStation1 + '.' + namStation1 + '.' + locStation1 + '.' + chaStation1 

#   now get some station meta data, esp. lat & long 
namStation2 = 'R6C8A'    #   MMU Oxford Road, Manchester  
netStation2 = 'AM'
locStation2 = '00'
chaStation2 = 'EHZ'      #   vertical component of geophone 
client2 = Client('RASPISHAKE')
idSEED2 = netStation2 + '.' + namStation2 + '.' + locStation2 + '.' + chaStation2 

#   now get some station meta data, esp. lat & long 
namStation3 = 'RAEED'    #   UoM Oxford Road, Manchester  
netStation3 = 'AM'
locStation3 = '00'
chaStation3 = 'EHZ'      #   vertical component of geophone 
client3 = Client('RASPISHAKE')
idSEED3 = netStation3 + '.' + namStation3 + '.' + locStation3 + '.' + chaStation3 

#   event datetime, e.g. from USGS or BGS
EQ_TIME1 = "2022-06-15T00:00:00"
dtEvent1 = UTCDateTime(EQ_TIME1)
t1start = dtEvent1 - T_START
t1end = dtEvent1 + T_END

# Download and filter data for station1
st1 = client1.get_waveforms(netStation1, namStation1, locStation1, chaStation1,
                          starttime=t1start, endtime=t1end, attach_response=True)
st1.merge(method=0, fill_value='latest')
st1.detrend(type="demean")
st1.remove_response()

st1f = st1.copy()
st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
#st1f.filter("highpass", freq=Flow, corners=2)
st1.trim(t1start, t1end)
st1f.trim(t1start, t1end)

st2 = client1.get_waveforms(netStation2, namStation2, locStation2, chaStation2,
                          starttime=t1start, endtime=t1end, attach_response=True)
st2.merge(method=0, fill_value='latest')
st2.detrend(type="demean")
st2.remove_response()

st2f = st2.copy()
st2f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
#st2f.filter("highpass", freq=Flow, corners=2)
st2.trim(t1start, t1end)
st2f.trim(t1start, t1end)

st3 = client1.get_waveforms(netStation3, namStation3, locStation3, chaStation3,
                          starttime=t1start, endtime=t1end, attach_response=True)
st3.merge(method=0, fill_value='latest')
st3.detrend(type="demean")
st3.remove_response()

st3f = st3.copy()
st3f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
#st3f.filter("highpass", freq=Flow, corners=2)
st3.trim(t1start, t1end)
st3f.trim(t1start, t1end)

VSCALE = 1e-4 

# Now plot the waveform data
st1f.plot(type="dayplot", interval=60, right_vertical_labels=False,
         vertical_scaling_range=VSCALE, one_tick_per_line=True,
         color=['k', 'r', 'b', 'g'], show_y_UTC_label=False,
         outfile="R4DD7 helicorder.png")

st2f.plot(type="dayplot", interval=60, right_vertical_labels=False,
         vertical_scaling_range=VSCALE, one_tick_per_line=True,
         color=['k', 'r', 'b', 'g'], show_y_UTC_label=False,
         outfile="R6C8A helicorder.png")

st3f.plot(type="dayplot", interval=60, right_vertical_labels=False,
         vertical_scaling_range=VSCALE, one_tick_per_line=True,
         color=['k', 'r', 'b', 'g'], show_y_UTC_label=False,
         outfile="RAEED helicorder.png")
