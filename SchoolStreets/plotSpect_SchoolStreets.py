#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 09:37:58 2022

@author: davidhealy
"""
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from scipy import signal # further FFT functionality

#   now get some station meta data, esp. lat & long 
#namStation = 'R1770'
#labelStation = 'Brokburn School, Chorlton'
#namStation = 'RA4D0'
#labelStation = 'Claude Road, Chorlton'
namStation = 'R0174'
labelStation = 'Ivy Green, Chorlton'
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

# now get the stream for the School Streets day
#sStart = "2023-02-02T00:00:00" # Week 1 
#sStart = "2023-02-09T07:50:00" # Week 2 
sStart = "2023-02-16T07:50:00" # Week 3 
dtStart = UTCDateTime(sStart)
T_START = 0
t1start = dtStart + T_START
T_END = 60*60*2
t1end = dtStart + T_END
print(t1start, t1end)

print("Getting stream for School Streets (spectrogram)...")
#   get the response for this station
resp = client.get_stations(t1start, network=netStation, station=namStation, location=locStation,
                  channel=chaStation, level="response")
# Download and filter data for station1
st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                      starttime=t1start, endtime=t1end, attach_response=True)
st.merge(method=0, fill_value='latest')

#   look at the spectrogram for the School Streets day - any issues?
window_size = 512
recording_rate = 100
frequencies, times, amplitudes = signal.spectrogram(st[0].data, fs=recording_rate, 
                                                    window='hamming', nperseg=window_size, 
                                                    noverlap=None, 
                                                    detrend=False, scaling="density")
index2 = pd.DatetimeIndex([(dtStart + ns).datetime for ns in times])
decibels = 20 * np.log10(amplitudes)

fig, ax = plt.subplots(figsize=(10, 3))
pcm = ax.pcolormesh(index2, frequencies, decibels, 
                        cmap="plasma", 
                        shading='auto', 
                        vmin=40, vmax=120)
ax.set_ylabel("Frequency (Hz)")
ax.set_ylim(1., 50.)
#plt.colorbar(pcm, ax=ax, orientation='horizontal', label='dB', location='top')
fig.autofmt_xdate()
plt.tight_layout() 
outfile='Spectrogram_SchoolStreets_Morning' + namStation + '_wk3.png'
plt.savefig(outfile, dpi=300)

# now get the stream for the School Streets day
#sStart = "2023-02-02T00:00:00" # Week 1 
#sStart = "2023-02-09T14:30:00" # Week 2 
sStart = "2023-02-16T14:30:00" # Week 3 
dtStart = UTCDateTime(sStart)
T_START = 0
t1start = dtStart + T_START
T_END = 60*60*2
t1end = dtStart + T_END
print(t1start, t1end)

print("Getting stream for School Streets (spectrogram)...")
#   get the response for this station
resp = client.get_stations(t1start, network=netStation, station=namStation, location=locStation,
                  channel=chaStation, level="response")
# Download and filter data for station1
st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                      starttime=t1start, endtime=t1end, attach_response=True)
st.merge(method=0, fill_value='latest')

#   look at the spectrogram for the School Streets day - any issues?
window_size = 512
recording_rate = 100
frequencies, times, amplitudes = signal.spectrogram(st[0].data, fs=recording_rate, 
                                                    window='hamming', nperseg=window_size, 
                                                    noverlap=None, 
                                                    detrend=False, scaling="density")
index2 = pd.DatetimeIndex([(dtStart + ns).datetime for ns in times])
decibels = 20 * np.log10(amplitudes)

fig, ax = plt.subplots(figsize=(10, 3))
pcm = ax.pcolormesh(index2, frequencies, decibels, 
                        cmap="plasma", 
                        shading='auto', 
                        vmin=40, vmax=120)
ax.set_ylabel("Frequency (Hz)")
ax.set_ylim(1., 50.)
#plt.colorbar(pcm, ax=ax, orientation='horizontal', label='dB', location='top')
fig.autofmt_xdate()
plt.tight_layout() 
outfile='Spectrogram_SchoolStreets_Afternoon' + namStation + '_wk3.png'
plt.savefig(outfile, dpi=300)
