#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 11:43:12 2021

@author: davidhealy
"""
import numpy as np
import matplotlib.pyplot as plt

from scipy import signal # further FFT functionality
from obspy.clients.fdsn import Client
from obspy import UTCDateTime 
from scipy.signal import find_peaks

T_START = 0     #   length in seconds of data to plot before origin time
T_END = 60*60    #   length in seconds of data to plot after origin time
Flow = 10. 
Fhigh = 20. 

#   now get some station meta data, esp. lat & long 
namStation = 'RA4D0'    #   Chorlton, Claude Road    
#namStation = 'R0174'    #   Chorlton, Ivy Green    
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

#   start of 'event' 
EQ_TIME = "2022-10-19T14:00:00"
dtEvent = UTCDateTime(EQ_TIME)
sEventName1 = 'Raw data'
sEventName2 = 'Filtered data'

#   get the data
#   plot seismograms  
t1 = dtEvent + T_START
t2 = dtEvent + T_END

# Download and filter data for this station
st1 = client.get_waveforms(netStation, namStation, locStation, chaStation,
                          starttime=t1, endtime=t2, attach_response=True)
st1.merge(method=0, fill_value='latest')

st1v = st1.copy()
st1v.detrend(type="demean")
st1v.remove_response()
st1v.trim(t1, t2)

st1f = st1v.copy()
#st1f.filter("highpass", freq=1., corners=2)
st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
st1f.trim(t1, t2)

window_size = 512
recording_rate = 100
frequencies, times, amplitudes = signal.spectrogram(st1[0].data, fs=recording_rate, 
                                                    window='hamming', nperseg=window_size, 
                                                    noverlap=window_size - 100, 
                                                    detrend=False, scaling="density")
decibels = 20 * np.log10(amplitudes)

# Now plot the waveform data
fig, axs = plt.subplots(2, 1, figsize=(36,9))

# axs[0].plot(st1[0].times(reftime=dtEvent+1), st1[0].data, linewidth=1)
# axs[0].grid(True)
# axs[0].set_ylabel("Counts")
# axs[0].set_xlim(T_START+1, T_END)
# axs[0].set_ylim(-1.e4, +4.e4)
# plt.xticks([0, 300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700])

# axs[0].plot(st1f[0].times(reftime=dtEvent+1), (st1f[0].data*1000)**2., linewidth=1)
# axs[0].grid(True)
# axs[0].set_ylabel("Velocity squared (mm/s)^2")
# axs[0].set_xlim(T_START+1., T_END)
# axs[0].set_ylim(0., 0.001)

pcm = axs[0].pcolormesh(times, frequencies, decibels, 
                        cmap="plasma", 
                        shading='auto',
                        vmin=40, vmax=110)
axs[0].set_ylabel("Frequency (Hz)")
axs[0].set_xlim(T_START+1., T_END)
axs[0].set_ylim(1.0, 20.)
#plt.colorbar(pcm, ax=axs[1], orientation='vertical')

peaks = find_peaks((st1f[0].data*1000)**2., 
                   distance=400, height=0.00015)

axs[1].plot(st1f[0].times(reftime=dtEvent+1), (st1f[0].data*1000)**2.)
[axs[1].axvline(st1f[0].times()[p], c='C3', linewidth=0.3) for p in peaks[0]]
#axs[1].set_xlim(0, len(st1f[0].times()))
axs[1].set_xlabel("Time from start (mins)")
axs[1].set_ylabel("Velocity squared (mm/s)^2")
axs[1].set_xlim(T_START+1., T_END)
axs[1].set_ylim(0., 0.001)
axs[1].set_title('Number of peaks: %3i' % len(peaks[0]))

plt.setp(axs, xticks=[0, 300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700], 
              xticklabels=['0', '5', '10', '15', '20', '25', '30', '35', '40', '45'])
plt.tight_layout() 
fn = namStation + '_SchoolStreets_19Oct22_PM.png'
plt.savefig(fn, dpi=300)

print(len(peaks[0])) 
      
