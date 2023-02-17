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
from obspy.signal import PPSD
from scipy import signal # further FFT functionality

nDays = 5 
T_START = 0                 # length in seconds of data to plot before origin time
T_END = 60*60*24*nDays      # length in seconds of data to plot after origin time
Fhigh = 1.                  # high-pass filter low end cutoff, i.e. only frequencies above this 

#   now get some station meta data, esp. lat & long 
namStation = 'RA4D0'
labelStation = 'Claude Road, Chorlton'
#namStation = 'R0174'
#labelStation = 'Ivy Green, Chorlton'
sTitle = namStation + ' - ' + labelStation
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

#   define the time window
sStart = "2023-01-23T00:00:00" # Spring
dtStart = UTCDateTime(sStart)
t1start = dtStart + T_START
t1end = dtStart + T_END + 60*30
print(t1start, t1end)

#   loop through the stations 
print("Getting streams for background...")
#   get the response for this station
resp = client.get_stations(t1start, network=netStation, station=namStation, location=locStation,
                  channel=chaStation, level="response")
# Download and filter data for station1
st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                      starttime=t1start, endtime=t1end, attach_response=True)
stf = st.copy()
stf.filter("highpass", freq=Fhigh, corners=2)

#   build a new PPSD object for this station, this interval 
binSizeMinutes = 1  
#   length of 600 = 10 mins (in seconds), overlap of 0.5 = 5 minute bins
ppsd = PPSD(stf[0].stats, resp,
            ppsd_length=binSizeMinutes*60*2, overlap=0.5,
            period_smoothing_width_octaves=0.025,
            period_step_octaves=0.0125,
            period_limits=(0.025, 50),
            db_bins=(-200, 20, 0.25))
ppsd.add(stf)

selectedPeriods = [1./30., 1./20., 1./10.]   # inverse of frequencies at 30, 20, 10 Hz
fig = ppsd.plot_temporal(selectedPeriods, legend=False, grid=True, show=False)
fig.set_size_inches(8, 4)
ax = fig.axes[0]
line20Hz = ax.lines[1]
line10Hz = ax.lines[2]
x10Hz = line10Hz.get_xdata()
y10Hz = line10Hz.get_ydata()
x20Hz = line20Hz.get_xdata()
y20Hz = line20Hz.get_ydata()
ax.set_xlabel("Date & time")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
ax.set_title(r"$\bf{" + sTitle + "}$ \n5 days from " + sStart[:10], loc='left')
plt.legend(['30 Hz', '20 Hz', '10 Hz'])
plt.xlim(t1start, t1end)
plt.ylim(-120, -50)
#plt.tight_layout()
outfile='SchoolStreets_PSDfreq_' + namStation + '_5days.png'
plt.savefig(outfile, dpi=300)
plt.show()

binPerHour = 60/binSizeMinutes 
nBins = int(binPerHour * 24 - 1)

#   loop approach 
dBdata10Hz = np.zeros([nDays, nBins])
dBdata20Hz = np.zeros([nDays, nBins])
ix = 0  
for d in range(0, nDays):
    for t in range(0, nBins):
        dBdata10Hz[d,t] = y10Hz[ix]
        dBdata20Hz[d,t] = y20Hz[ix]
        ix += 1 

dBdatamedian10Hz = np.zeros([nBins])
dBdatamedian20Hz = np.zeros([nBins])
dBdatastd10Hz = np.zeros([nBins])
dBdatastd20Hz = np.zeros([nBins])
for t in range(0, nBins):
    dBdatamedian10Hz[t] = np.median(dBdata10Hz[0:,t])
    dBdatamedian20Hz[t] = np.median(dBdata20Hz[0:,t])
    dBdatastd10Hz[t] = np.std(dBdata10Hz[0:,t])
    dBdatastd20Hz[t] = np.std(dBdata20Hz[0:,t])
    
# now get the stream for the School Streets day
sStart = "2023-02-02T08:30:00" # Week 1 
#sStart = "2023-02-09T08:30:00" # Week 2 
#sStart = "2023-02-16T08:30:00" # Week 3 
dtStart = UTCDateTime(sStart)
t1start = dtStart + T_START
T_END = 60*40
t1end = dtStart + T_END
print(t1start, t1end)

print("Getting stream for School Streets spectrogram AM...")
#   get the response for this station
resp = client.get_stations(t1start, network=netStation, station=namStation, location=locStation,
                  channel=chaStation, level="response")
# Download and filter data for station1
st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                      starttime=t1start, endtime=t1end, attach_response=True)

#   look at the spectrogram for the School Streets day - any issues?
window_size = 512
recording_rate = 100
frequencies, times, amplitudes = signal.spectrogram(st[0].data, fs=recording_rate, 
                                                    window='hamming', nperseg=window_size, 
                                                    noverlap=window_size - 100, 
                                                    detrend=False, scaling="density")
index2 = pd.DatetimeIndex([(dtStart + ns).datetime for ns in times])
decibels = 20 * np.log10(amplitudes)

fig, ax = plt.subplots(figsize=(7, 4))
pcm = ax.pcolormesh(index2, frequencies, decibels, 
                        cmap="plasma", 
                        shading='auto', 
                        vmin=40, vmax=120)
ax.set_ylabel("Frequency (Hz)")
ax.set_ylim(1., 50.)
ax.set_title(r"$\bf{" + sTitle + "}$ \nMorning session - " + sStart[:10], loc='left')
#plt.colorbar(pcm, ax=ax, orientation='horizontal', label='dB', location='top')
fig.autofmt_xdate()
plt.tight_layout() 
outfile='SchoolStreets_Spectrogram_AM_' + namStation + '_' + sStart[:10] + '.png'
plt.savefig(outfile, dpi=300)

del st 

# now get the stream for the School Streets day
sStart = "2023-02-02T15:10:00" # Week 1 
#sStart = "2023-02-09T15:10:00" # Week 2 
#sStart = "2023-02-16T15:10:00" # Week 3 
dtStart = UTCDateTime(sStart)
t1start = dtStart + T_START
T_END = 60*40
t1end = dtStart + T_END
print(t1start, t1end)

print("Getting stream for School Streets spectrogram PM...")
#   get the response for this station
resp = client.get_stations(t1start, network=netStation, station=namStation, location=locStation,
                  channel=chaStation, level="response")
# Download and filter data for station1
st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                      starttime=t1start, endtime=t1end, attach_response=True)

#   look at the spectrogram for the School Streets day - any issues?
window_size = 512
recording_rate = 100
frequencies, times, amplitudes = signal.spectrogram(st[0].data, fs=recording_rate, 
                                                    window='hamming', nperseg=window_size, 
                                                    noverlap=window_size - 100, 
                                                    detrend=False, scaling="density")
index2 = pd.DatetimeIndex([(dtStart + ns).datetime for ns in times])
decibels = 20 * np.log10(amplitudes)

fig, ax = plt.subplots(figsize=(7, 4))
pcm = ax.pcolormesh(index2, frequencies, decibels, 
                        cmap="plasma", 
                        shading='auto', 
                        vmin=40, vmax=120)
ax.set_ylabel("Frequency (Hz)")
ax.set_ylim(1., 50.)
ax.set_title(r"$\bf{" + sTitle + "}$ \nAfternoon session - " + sStart[:10], loc='left')
#plt.colorbar(pcm, ax=ax, orientation='horizontal', label='dB', location='top')
fig.autofmt_xdate()
plt.tight_layout() 
outfile='SchoolStreets_Spectrogram_PM_' + namStation + '_' + sStart[:10] + '.png'
plt.savefig(outfile, dpi=300)

del st 

sStart = "2023-02-02T00:00:00" # Week 1 
#sStart = "2023-02-09T00:00:00" # Week 2 
#sStart = "2023-02-16T00:00:00" # Week 3 
dtStart = UTCDateTime(sStart)
t1start = dtStart + T_START
T_END = 60*60*24
t1end = dtStart + T_END
print(t1start, t1end)

print("Getting stream for School Streets...")
#   get the response for this station
resp = client.get_stations(t1start, network=netStation, station=namStation, location=locStation,
                  channel=chaStation, level="response")
# Download and filter data for station1
st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                      starttime=t1start, endtime=t1end, attach_response=True)
stf = st.copy()
stf.filter("highpass", freq=Fhigh, corners=2)

#   build a new PPSD object for this station, this interval 
binSizeMinutes = 1  
#   1 minute bins (2 minute length, 50% overlap)
ppsd = PPSD(stf[0].stats, resp,
            ppsd_length=binSizeMinutes*60*2, overlap=0.5,
            period_smoothing_width_octaves=0.025,
            period_step_octaves=0.0125,
            period_limits=(0.025, 50),
            db_bins=(-200, 20, 0.25))
ppsd.add(stf)

binPerHour = 60/binSizeMinutes 
nBins = int(binPerHour * 24 - 1)

x20Hz = np.zeros(nBins)
y20Hz = np.zeros(nBins)

selectedPeriods = [1./30., 1./20., 1./10.]   # inverse of frequencies at 30, 20, 10 Hz
fig = ppsd.plot_temporal(selectedPeriods, legend=False, grid=True, show=False)
fig.set_size_inches(7, 4)
ax = fig.axes[0]
line20Hz = ax.lines[1]
x20Hz[0:len(line20Hz.get_xdata())] = line20Hz.get_xdata()
y20Hz[0:len(line20Hz.get_ydata())] = line20Hz.get_ydata()
ax.set_xlabel("Date & time")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
ax.set_title(r"$\bf{" + sTitle + "}$ \nWhole day - " + sStart[:10], loc='left')
plt.legend(['30 Hz', '20 Hz', '10 Hz'])
plt.xlim(t1start, t1end)
plt.ylim(-120, -60)
plt.tight_layout()
outfile='SchoolStreets_PSDfreq_' + namStation + '_' + sStart[:10] + '_24hours.png'
plt.savefig(outfile, dpi=300)
plt.show()

#   loop approach 
dBdataSS20Hz = np.zeros([nBins])
ix = 0  
for t in range(0, nBins):
    dBdataSS20Hz[t] = y20Hz[ix]
    ix += 1 

xBins = list(range(1,nBins+1))

fig, ax = plt.subplots(figsize=(7,4)) 
plt.fill_between(xBins, 
                 dBdatamedian20Hz - dBdatastd20Hz, 
                 dBdatamedian20Hz + dBdatastd20Hz,
                 facecolor='lightblue', label='5-day envelope')
ax.plot(xBins, dBdataSS20Hz, '-r', label='Closure day (20 Hz)')
ax.set_xlabel("Hour of day")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
ax.set_title(r"$\bf{" + sTitle + "}$ \nWhole day - " + sStart[:10], loc='left')
plt.xticks(ticks=[0,4*binPerHour,8*binPerHour,12*binPerHour,16*binPerHour,20*binPerHour,24*binPerHour], labels=['0', '4', '8', '12', '16', '20', '24'])
plt.legend()
plt.xlim(0, max(xBins))
plt.ylim(-115, -60)
plt.tight_layout()
outfile='SchoolStreets_PSD20Hz_' + namStation + '_' + sStart[:10] + '_24hours.png'
plt.savefig(outfile, dpi=300)

fig, ax = plt.subplots(figsize=(7,4)) 
plt.fill_between(xBins, 
                 dBdatamedian20Hz - dBdatastd20Hz, 
                 dBdatamedian20Hz + dBdatastd20Hz,
                 facecolor='lightblue', label='5-day envelope (20 Hz)')
ax.plot(xBins, dBdataSS20Hz, '-r', label='Closure day (20 Hz)')
ax.axvline((8.5*binPerHour), 0, 1, linestyle='-.')
ax.axvline((9.17*binPerHour), 0, 1, linestyle='-.')
ax.set_xlabel("Time of day")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
ax.set_title(r"$\bf{" + sTitle + "}$ \nMorning session - " + sStart[:10], loc='left')
plt.xticks(ticks=[8*binPerHour,9*binPerHour,9.83*binPerHour], labels=['8:00 am', '9 am', '9:50 am'])
plt.legend()
plt.xlim(7.83*binPerHour, 9.83*binPerHour)
plt.ylim(-100, -60)
plt.tight_layout()
outfile='SchoolStreets_PSD20Hz_' + namStation + '_' + sStart[:10] + '_MorningClosure.png'
plt.savefig(outfile, dpi=300)

fig, ax = plt.subplots(figsize=(7,4)) 
plt.fill_between(xBins, 
                 dBdatamedian20Hz - dBdatastd20Hz, 
                 dBdatamedian20Hz + dBdatastd20Hz,
                 facecolor='lightblue', label='5-day envelope (20 Hz)')
ax.plot(xBins, dBdataSS20Hz, '-r', label='Closure day (20 Hz)')
ax.axvline((15.17*binPerHour), 0, 1, linestyle='-.')
ax.axvline((15.83*binPerHour), 0, 1, linestyle='-.')
ax.set_xlabel("Time of day")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.grid(True)
ax.set_title(r"$\bf{" + sTitle + "}$ \nAfternoon session - " + sStart[:10], loc='left')
plt.xticks(ticks=[14.5*binPerHour, 15*binPerHour, 16*binPerHour, 16.5*binPerHour], labels=['2:30 pm', '3 pm', '4 pm', '4:30 pm'])
plt.legend()
plt.xlim(14.5*binPerHour, 16.5*binPerHour)
plt.ylim(-100, -60)
plt.tight_layout()
outfile='SchoolStreets_PSD20Hz_' + namStation + '_' + sStart[:10] + '_AfternoonClosure.png'
plt.savefig(outfile, dpi=300)
