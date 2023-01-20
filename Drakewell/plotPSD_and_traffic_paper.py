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
#namStation = 'R36CB'
#labelStation = 'Cheadle Hulme'
namStation = 'R3FEA'
labelStation = 'Chorlton'
# namStation = 'R4DD7'
# labelStation = 'Fallowfield'
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

#   define the time window
EQ_TIME = "2022-04-25T00:00:00" 
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
fig, axs = plt.subplots(1, 4, figsize=(18,4)) 
axs[0].boxplot(dBdata, whis=(5,95))
axs[0].plot(xHalfhours, dBdatamedian, '-r', label='10 Hz')
axs[0].set_xlabel("Hour of day")
axs[0].set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
axs[0].grid(True)
axs[0].legend()
#axs[0].xticks(ticks=[0,8,16,24,32,40], labels=['0', '4', '8', '12', '16', '20'])
axs[0].set_ylim(-105, -70)

#   now the traffic data from Drakewell camera for same 5 days 
#   read in data file, in half-hour increments 
#   columns for light vehciles, heavy vehicles, total vehicles 
#   reshape data into an array based on 24 hours, 47 half-hours 
#   find median value at each half-hour over the 5 days 

from datetime import datetime
#   read in file of traffic counts  
#fnTC = 'Drakewell 1426 - 30 minute data - 25Apr22-29Apr22.csv' # R36CB Cheadle Hulme 
fnTC = 'Drakewell 1312 - 30 minute data - 25Apr22-29Apr22.csv' # R3FEA Chorlton, Oswald Rd 
#fnTC = 'Drakewell 1017 - 30 minute data - 05Sep22-09Sep22.csv' # R4DD7 Fallowfield 
iLine = 0 
#iColHeavy = 5
#iColMax = 7 
iColHeavy = 8
iColMax = 14 
#   open the file 
with open(fnTC, 'r') as reader:
    dtList = [] 
    fTotalCountList = []
    fHeavyCountList = []
    for line in reader:
        iLine += 1 
        if iLine >= 5:
            sLineTokens = line.split(',')
            #print(sLineTokens)
            fCount = 0. 
            fHeavyCount = 0. 
            for i, s in enumerate(sLineTokens):
                if i == 0:
                    s = s.replace('"', '')
                    s = s.strip()
                    dDateTime = datetime.strptime(s, 
                                '%d/%m/%Y %H:%M')
                    dtList.append(dDateTime)
                elif i <= iColMax: 
                    fCount += float(s)
                    if i >= iColHeavy and i <= iColMax:
                        fHeavyCount += float(s) 
                else:
                    continue 
            fTotalCountList.append(fCount)
            fHeavyCountList.append(fHeavyCount)
            
print("Read in %i rows of data." % iLine)
fLightCountList = []
for i in range(len(fTotalCountList)):
    fLightCountList.append(fTotalCountList[i] - fHeavyCountList[i]) 

#   loop approach 
lightTrafficdata = np.zeros([nDays, 47])
heavyTrafficdata = np.zeros([nDays, 47])
ix = 0  
for d in range(0, nDays):
    for t in range(0, 47):
        lightTrafficdata[d,t] = fLightCountList[ix]
        heavyTrafficdata[d,t] = fHeavyCountList[ix]
        ix += 1 

lightTrafficdatamedian = np.zeros([47,1])
heavyTrafficdatamedian = np.zeros([47,1])
for t in range(0, 47):
    lightTrafficdatamedian[t] = np.median(lightTrafficdata[0:,t])
    heavyTrafficdatamedian[t] = np.median(heavyTrafficdata[0:,t])

axs[1].boxplot(heavyTrafficdata, whis=(5,95))
axs[1].plot(xHalfhours, heavyTrafficdatamedian, '-g', label='Heavy vehicles') 
axs[1].grid(True)
axs[1].set_ylabel('Number of vehicles')
axs[1].legend()
axs[1].set_xlim(0, 47)
axs[1].set_ylim(60, 120)
#axs[1].xticks(ticks=[0,8,16,24,32,40], labels=['0', '4', '8', '12', '16', '20'])

axs[2].boxplot(lightTrafficdata, whis=(5,95))
axs[2].plot(xHalfhours, lightTrafficdatamedian, '-b', label='Light vehicles') 
axs[2].grid(True)
axs[2].set_ylabel('Number of vehicles')
axs[2].legend()
axs[2].set_xlim(0, 47)
axs[2].set_ylim(0, 700)
#axs[2].xticks(ticks=[0,8,16,24,32,40], labels=['0', '4', '8', '12', '16', '20'])

axs[3].boxplot((heavyTrafficdata+lightTrafficdata), whis=(5,95))
axs[3].plot(xHalfhours, heavyTrafficdatamedian+lightTrafficdatamedian, '-k', label='All vehicles') 
axs[3].grid(True)
axs[3].set_ylabel('Number of vehicles')
axs[3].legend()
axs[3].set_xlim(0, 47)
axs[3].set_ylim(0, 800)
#axs[3].xticks(ticks=[0,8,16,24,32,40], labels=['0', '4', '8', '12', '16', '20'])

plt.setp(axs, xticks=[0,8,16,24,32,40], xticklabels=['0', '4', '8', '12', '16', '20'])
plt.tight_layout()
fn = namStation + '_PSD_and_Traffic.png'
plt.savefig(fn, dpi=300)
plt.show()
