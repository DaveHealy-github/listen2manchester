#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 09:37:58 2022

@author: davidhealy
"""
import numpy as np
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.signal import PPSD

nDays = 5 
T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*60*24*nDays    # length in seconds of data to plot after origin time
Fhigh = 1.   

#   now get some station meta data, esp. lat & long 
namStation = 'R36CB'
labelStation = 'Cheadle Hulme'
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

#   define the time window
#EQ_TIME = "2022-08-22T00:00:00" # Summer
EQ_TIME = "2022-04-25T00:00:00" # Spring
dtEvent = UTCDateTime(EQ_TIME)
t1start = dtEvent + T_START
t1end = dtEvent + T_END - 1 
print(t1start, t1end)

#   loop through the stations 
print("Getting streams...")
#   get the response for this station
resp = client.get_stations(t1start, network=netStation, station=namStation, location=locStation,
                  channel=chaStation, level="response")
# Download and filter data for station1
st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                      starttime=t1start, endtime=t1end, attach_response=True)
stf = st.copy()
stf.filter("highpass", freq=Fhigh, corners=2)

#   build a new PPSD object for this station, this interval 
ppsd = PPSD(stf[0].stats, resp,
            ppsd_length=3600, overlap=0.5,
            period_smoothing_width_octaves=0.025,
            period_step_octaves=0.0125,
            period_limits=(0.025, 50),
            db_bins=(-200, 20, 0.25))
ppsd.add(stf)

fig = ppsd.plot_temporal([1./30., 1./10., 1./5.], legend=False, grid=True, show=False)
fig.set_size_inches(8, 4)
ax = fig.axes[0]
line10Hz = ax.lines[1]
x10Hz = line10Hz.get_xdata()
y10Hz = line10Hz.get_ydata()
ax.set_xlabel("Date & time")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
plt.legend(['30 Hz', '10 Hz', '5 Hz'])
#plt.ylim(-125, -75)
plt.tight_layout()
outfile='SelectedPSD_weekdays_paper.png'
plt.savefig(outfile, dpi=300)
plt.show()

#   loop approach 
dBdata = np.zeros([nDays, 47])
ix = 0  
for d in range(0, nDays):
    for t in range(0, 47):
        dBdata[d,t] = y10Hz[ix]
        ix += 1 

dBdatamedian = np.zeros([47,1])
for t in range(0, 47):
    dBdatamedian[t] = np.median(dBdata[0:,t])
    
xHalfhours = list(range(1,48))
fig, ax = plt.subplots(figsize=(4.5,4)) 
ax.boxplot(dBdata, sym='')
ax.plot(xHalfhours, dBdatamedian, '-r', label='10 Hz')
ax.set_xlabel("Hour of day")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
plt.xticks(ticks=[0,8,16,24,32,40], labels=['0', '4', '8', '12', '16', '20'])
plt.legend()
plt.ylim(-104, -81)
plt.tight_layout()
outfile='SelectedPSD_weekdays_halfhourly_paper_' + namStation + '.png'
plt.savefig(outfile, dpi=300)

del fig, resp, st, stf, ppsd, x10Hz, y10Hz, dBdata, dBdatamedian

#   define the 1st weekend time window
nDays = 2 
#EQ_TIME = "2022-08-20T00:00:00" # Summer
EQ_TIME = "2022-04-23T00:00:00" # Spring 
T_END = 60*60*24*nDays    # length in seconds of data to plot after origin time
dtEvent = UTCDateTime(EQ_TIME)
t1start = dtEvent + T_START
t1end = dtEvent + T_END - 1 
print(t1start, t1end)

#   loop through the stations 
print("Getting weekend 1 streams...")
#   get the response for this station
resp = client.get_stations(t1start, network=netStation, station=namStation, location=locStation,
                  channel=chaStation, level="response")
# Download and filter data for station1
st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                      starttime=t1start, endtime=t1end, attach_response=True)
stf = st.copy()
stf.filter("highpass", freq=Fhigh, corners=2)

#   build a new PPSD object for this station, this interval 
ppsd = PPSD(stf[0].stats, resp,
            ppsd_length=3600, overlap=0.5,
            period_smoothing_width_octaves=0.025,
            period_step_octaves=0.0125,
            period_limits=(0.025, 50),
            db_bins=(-200, 20, 0.25))
ppsd.add(stf)

#   define the 2nd weekend time window
nDays = 3 
#EQ_TIME = "2022-08-27T00:00:00" # Summer
EQ_TIME = "2022-04-30T00:00:00" # Spring 
T_END = 60*60*24*nDays    # length in seconds of data to plot after origin time
dtEvent = UTCDateTime(EQ_TIME)
t1start = dtEvent + T_START
t1end = dtEvent + T_END - 1 
print(t1start, t1end)

#   loop through the stations 
print("Getting weekend 2 streams...")
# Download and filter data for station1
st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                      starttime=t1start, endtime=t1end, attach_response=True)
stf = st.copy()
stf.filter("highpass", freq=Fhigh, corners=2)

ppsd.add(stf)

fig = ppsd.plot_temporal([1./30., 1./10., 1./5.], legend=False, grid=True, show=False)
fig.set_size_inches(8, 4)
ax = fig.axes[0]
print(len(ax.lines))
line10HzPart1 = ax.lines[2]
line10HzPart2 = ax.lines[3]
x10Hz = np.concatenate((line10HzPart1.get_xdata(), line10HzPart2.get_xdata()))
y10Hz = np.concatenate((line10HzPart1.get_ydata(), line10HzPart2.get_ydata()))
ax.set_xlabel("Date & time")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
plt.legend(['30 Hz', '10 Hz', '5 Hz'])
#plt.ylim(-125, -75)
plt.tight_layout()
outfile='SelectedPSD_weekends_paper.png'
plt.savefig(outfile, dpi=300)
plt.show()

#   loop approach 
nDays = 5 
dBdata = np.zeros([nDays, 47])
ix = 0  
for d in range(0, nDays):
    for t in range(0, 47):
        dBdata[d,t] = y10Hz[ix]
        ix += 1 

dBdatamedian = np.zeros([47,1])
for t in range(0, 47):
    dBdatamedian[t] = np.median(dBdata[0:,t])
    
xHalfhours = list(range(1,48))
fig, ax = plt.subplots(figsize=(4.5,4)) 
ax.boxplot(dBdata, sym='')
ax.plot(xHalfhours, dBdatamedian, '-r', label='10 Hz')
ax.set_xlabel("Hour of day")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
plt.xticks(ticks=[0,8,16,24,32,40], labels=['0', '4', '8', '12', '16', '20'])
plt.legend()
plt.ylim(-104, -81)
plt.tight_layout()
outfile='SelectedPSD_weekends_halfhourly_paper_' + namStation + '.png'
plt.savefig(outfile, dpi=300)
