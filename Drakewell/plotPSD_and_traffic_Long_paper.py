#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 09:37:58 2022

@author: davidhealy
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd 
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.signal import PPSD
from scipy import signal # further FFT functionality
from matplotlib.dates import DateFormatter

nDays = 5 
T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*60*24*nDays    # length in seconds of data to plot after origin time
Fhigh = 1.   

#   now get some station meta data, esp. lat & long 
namStation = 'R36CB'
labelStation = 'Cheadle Hulme'
#namStation = 'R3FEA'
#labelStation = 'Chorlton'
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
st.merge(method=0, fill_value='latest')
stf = st.copy()
stf.filter("highpass", freq=Fhigh, corners=2)

window_size = 512
recording_rate = 100
frequencies, times, amplitudes = signal.spectrogram(st[0].data, 
                                                    fs=recording_rate, 
                                                    window='hamming', 
                                                    nperseg=window_size, 
                                                    noverlap=window_size - 100, 
                                                    detrend=False, 
                                                    scaling='density')
decibels = 20. * np.log10(amplitudes)
index2 = pd.DatetimeIndex([(dtEvent + ns).datetime for ns in times])
days = mdates.DayLocator(interval=1)

fig, axs = plt.subplots(figsize=(8,4))
pcm = axs.pcolormesh(index2, frequencies, decibels, 
                        cmap="plasma", 
                        shading='auto',
                        vmin=-25, vmax=125)
axs.set_ylabel("Frequency (Hz)")
#axs.set_xlabel("Time from start (s)")
axs.set_ylim(1.0, 50.)
#plt.colorbar(pcm, ax=axs, orientation='horizontal')
axs.xaxis.set_major_locator(days)
axs.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
plt.xticks()
axs.set_xlim(t1start, t1end)
plt.tight_layout() 
plt.savefig('plotCheadle_Spect_5days_paper.png', dpi=300)

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
ax.legend(['30 Hz', '10 Hz', '5 Hz']) 
line10Hz = ax.lines[1]
x10Hz = line10Hz.get_xdata()
y10Hz = line10Hz.get_ydata()
plt.autoscale(enable=True, axis='x', tight=True)
ax.xaxis.set_major_locator(days)
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
ax.set_xlim(t1start, t1end)
fig.autofmt_xdate(rotation=0)
plt.tight_layout() 
plt.savefig('plotCheadle_SelFreq_5days_paper.png', dpi=300)

#   now the traffic data from Drakewell camera for same 5 days 
#   read in data file, in half-hour increments 
#   columns for light vehciles, heavy vehicles, total vehicles 
#   reshape data into an array based on 24 hours, 47 half-hours 
#   find median value at each half-hour over the 5 days 

from datetime import datetime
#   read in file of traffic counts  
fnTC = 'Drakewell 1426 - 30 minute data - 25Apr22-29Apr22.csv' # R36CB Cheadle Hulme 
#fnTC = 'Drakewell 1312 - 30 minute data - 25Apr22-29Apr22.csv' # R3FEA Chorlton, Oswald Rd 
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

fig, axs = plt.subplots(figsize=(8,4))
axs.plot(dtList, fLightCountList, label='Light vehicles', color='C3') 
axs2 = axs.twinx()
axs2.plot(dtList, fHeavyCountList, label='Heavy vehicles', color='C4') 
#axs.plot(dtList, fTotalCountList, label='Total vehicles', color='C5') 
axs.grid(True)
axs.set_ylabel('Number of vehicles')
#axs.set_title('Traffic counts, TfGM Drakewell camera 1426')
fig.legend() 
plt.autoscale(enable=True, axis='x', tight=True)
axs.xaxis.set_major_locator(days)
axs.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
axs.set_xlim(t1start, t1end)
plt.tight_layout() 
plt.savefig('plotCheadle_Traffic_5days_paper.png', dpi=300)

