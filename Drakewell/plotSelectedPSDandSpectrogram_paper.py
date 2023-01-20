#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 09:37:58 2022

@author: davidhealy
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.signal import PPSD
from scipy import signal # further FFT functionality

nDays = 10
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
EQ_TIME = "2022-04-23T00:00:00" # Spring
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

window_size = 256
recording_rate = 100
frequencies, times, amplitudes = signal.spectrogram(st[0].data, fs=recording_rate, 
                                                    window='hamming', nperseg=window_size, 
                                                    noverlap=window_size - 100, 
                                                    detrend=False, scaling="density")
index2 = pd.DatetimeIndex([(dtEvent + ns).datetime for ns in times])
decibels = 20 * np.log10(amplitudes)

fig, ax = plt.subplots(figsize=(10, 4))
pcm = ax.pcolormesh(index2, frequencies, decibels, 
                        cmap="plasma", 
                        shading='auto', 
                        vmin=20, vmax=90)
ax.set_ylabel("Frequency (Hz)")
ax.set_ylim(1., 50.)
plt.colorbar(pcm, ax=ax, orientation='horizontal', label='dB', location='top')
#fig.autofmt_xdate()
plt.tight_layout() 
outfile='SelectedPSDandSpectrogram_Spect_paper.png'
plt.savefig(outfile, dpi=300)

#   build a new PPSD object for this station, this interval 
ppsd = PPSD(stf[0].stats, resp,
            ppsd_length=3600, overlap=0.5,
            period_smoothing_width_octaves=0.025,
            period_step_octaves=0.0125,
            period_limits=(0.025, 50),
            db_bins=(-200, 20, 0.25))
ppsd.add(stf)

fig = ppsd.plot_temporal([1./30., 1./10., 1./5.], legend=False, grid=True, show=False)
fig.set_size_inches(10, 4)
ax = fig.axes[0]
line10Hz = ax.lines[1]
x10Hz = line10Hz.get_xdata()
y10Hz = line10Hz.get_ydata()
ax.set_xlabel("Date & time")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
plt.legend(['30 Hz', '10 Hz', '5 Hz'])
plt.xlim(t1start, t1end)
plt.ylim(-125, -80)
plt.tight_layout()
outfile='SelectedPSDandSpectrogram_paper.png'
plt.savefig(outfile, dpi=300)
plt.show()
