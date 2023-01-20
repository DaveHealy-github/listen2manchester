#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 11:43:12 2021

@author: davidhealy
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from numpy.fft import * # a library for carrying out Fourier Transforms
from scipy import signal # further FFT functionality

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

T_START = 0     #   length in seconds of data to plot before origin time
T_END = 60*10    #   length in seconds of data to plot after origin time
Flow = 0.7     #   low cut-off for bandpass filter 
Fhigh = 2.      #   high cut-off for bandpass filter 
Flow2 = 4.      #   low cut-off for bandpass filter 
Fhigh2 = 20.      #   high cut-off for bandpass filter 

#   now get some station meta data, esp. lat & long 
namStation = 'R9E71'    #   RHS Bridgewater   
#namStation = 'R72C7'    #   Aberdeen Uni   
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

#   start of 'event' 
EQ_TIME = "2021-11-16T01:40:00"
dtEvent = UTCDateTime(EQ_TIME)
sEventName1 = 'Raw data'
sEventName2 = 'Filtered data'

#   get the data
#   plot seismograms  
t1 = dtEvent - T_START
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
st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)
st1f.trim(t1, t2)

st1f2 = st1v.copy()
st1f2.filter("bandpass", freqmin=Flow2, freqmax=Fhigh2, corners=4)
st1f2.trim(t1, t2)

#   code from Richard Strange, Oxford 
y_values = np.fft.fft(st1[0].data)
no_of_datapoints = len(y_values)
time_interval = 0.01 #    sample interval for RS seismometers
yf_values = 2.0/no_of_datapoints * np.abs(y_values[:no_of_datapoints//2])
x_values = fftfreq(no_of_datapoints, d=time_interval)
xf_values = fftfreq(no_of_datapoints, d=time_interval)[:no_of_datapoints//2]

window_size = 512
recording_rate = 100
frequencies, times, amplitudes = signal.spectrogram(st1[0].data, fs=recording_rate, 
                                                    window='hamming', nperseg=window_size, 
                                                    noverlap=window_size - 100, 
                                                    detrend=False, scaling="density")
decibels = 20 * np.log10(amplitudes)

# Now plot the waveform data
fig, axs = plt.subplots(5, 1, figsize=(12,18))

axs[0].plot(st1[0].times(reftime=dtEvent+1), st1[0].data, linewidth=1)
axs[0].grid(True)
axs[0].set_xlabel("Time from start (s)")
axs[0].set_title("{:} - RaspberryShake {:}.{:}.{:}.{:} - raw".format(
    sEventName1, st1[0].stats.network, st1[0].stats.station, st1[0].stats.location,
    st1[0].stats.channel))
axs[0].set_ylabel("Counts")
axs[0].set_xlim(T_START+1, T_END)

axs[1].plot(st1f[0].times(reftime=dtEvent+1), st1f[0].data*1000, linewidth=1)
axs[1].grid(True)
axs[1].set_xlabel("Time from start (s)")
axs[1].set_title("{:} - RaspberryShake {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    sEventName2, st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
    st1f[0].stats.channel, Flow, Fhigh))
axs[1].set_ylabel("Ground velocity (mm/s)")
axs[1].set_xlim(T_START+1, T_END)

axs[2].plot(st1f2[0].times(reftime=dtEvent+1), st1f2[0].data*1000, linewidth=1)
axs[2].grid(True)
axs[2].set_xlabel("Time from start (s)")
axs[2].set_title("{:} - RaspberryShake {:}.{:}.{:}.{:} - bandpass filter: {:}-{:} Hz".format(
    sEventName2, st1f2[0].stats.network, st1f2[0].stats.station, st1f2[0].stats.location,
    st1f2[0].stats.channel, Flow2, Fhigh2))
axs[2].set_ylabel("Ground velocity (mm/s)")
axs[2].set_xlim(T_START+1, T_END)

pcm = axs[3].pcolormesh(times, frequencies, decibels, 
                        cmap="plasma", 
                        shading='auto',
                        vmin=-100, vmax=200)
axs[3].set_ylabel("Frequency (Hz)")
axs[3].set_xlabel("Time from start (s)")
axs[3].set_title("Spectrogram")
#plt.colorbar(pcm, ax=axs[2], orientation='horizontal')

axs[4].plot(xf_values, yf_values, lw=1)
axs[4].grid(True)
axs[4].set_xlabel("Frequency (Hz)")
axs[4].set_ylabel("Amplitude")
axs[4].set_title("Fast Fourier Transform")
axs[4].set_xlim(0, 50)
axs[4].set_ylim(0, 10)

plt.tight_layout() 
plt.savefig('plotSpectrogram6HourSlicesR9E71.png', dpi=300)
