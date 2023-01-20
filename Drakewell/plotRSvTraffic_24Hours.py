#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 11:43:12 2021

@author: davidhealy
"""
import os
import glob 
import obspy.signal 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib.dates import DateFormatter
from scipy import signal # further FFT functionality
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read 
from obspy.signal import PPSD
#from hampel import hampel

T_START = 0     #   length in seconds of data to plot before origin time
T_END = 60*60*24    #   length in seconds of data to plot after origin time

#   now get some station meta data, esp. lat & long 
#namStation = 'R6C8A'    #   MMU Oxford Road   
#namStation = 'R36CB'    #   Cheadle Hulme   
namStation = 'R3FEA'    #   Chorlton   
#namStation = 'R72C7'    #   Aberdeen Uni   
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

#   start of 'event' 
EQ_TIME = "2022-04-28T00:00:00"
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
st1f.filter("highpass", freq=1.0, corners=2)
st1f.trim(t1, t2)

window_size = 256
recording_rate = 100
frequencies, times, amplitudes = signal.spectrogram(st1[0].data, fs=recording_rate, 
                                                    window='hamming', nperseg=window_size, 
                                                    noverlap=window_size - 100, 
                                                    detrend=False, scaling="density")
decibels = 20 * np.log10(amplitudes)

# Now plot the waveform data
fig, axs = plt.subplots(3, 1, figsize=(12,12))

axs[0].plot(st1[0].times(reftime=dtEvent+1), st1[0].data, linewidth=1)
axs[0].grid(True)
axs[0].set_xlabel("Time from start (s)")
axs[0].set_title("{:} - RaspberryShake {:}.{:}.{:}.{:} - raw".format(
    sEventName1, st1[0].stats.network, st1[0].stats.station, st1[0].stats.location,
    st1[0].stats.channel))
axs[0].set_ylabel("Counts")
axs[0].set_xlim(T_START+1, T_END)
axs[0].set_ylim(-0.01e6, 0.04e6)

axs[1].plot(st1f[0].times(reftime=dtEvent+1), st1f[0].data*1000, linewidth=1)
axs[1].grid(True)
axs[1].set_xlabel("Time from start (s)")
axs[1].set_title("{:} - RaspberryShake {:}.{:}.{:}.{:} - highpass filtered 1 Hz".format(
    sEventName2, st1f[0].stats.network, st1f[0].stats.station, st1f[0].stats.location,
    st1f[0].stats.channel))
axs[1].set_ylabel("Velocity (mm/s)")
axs[1].set_xlim(T_START+1, T_END)
axs[1].set_ylim(-0.05, 0.05)

pcm = axs[2].pcolormesh(times, frequencies, decibels, 
                        cmap="plasma", 
                        shading='auto',
                        vmin=-25, vmax=125)
axs[2].set_ylabel("Frequency (Hz)")
axs[2].set_xlabel("Time from start (s)")
axs[2].set_title("Spectrogram")
axs[2].set_ylim(1.0, 50.)
#plt.colorbar(pcm, ax=axs[1], orientation='vertical')

plt.tight_layout() 
fn = namStation + '_counts_spectrogram.png'
plt.show() 
#plt.savefig(fn, dpi=300)

#   tidy up old work files
# for f in glob.glob("*.npz"):
#     os.remove(f)
# for f in glob.glob("*.mseed"):
#     os.remove(f)

#############################################################################
#   get a specified time window of RShake data from a station 
#############################################################################

#   define the time window of interest 
start = UTCDateTime("2022-04-28T00:00:00")
nDays = 1. 
twSpan = nDays*24.*60.*60.-1
end = start + twSpan
datelist = pd.date_range(start.datetime, end.datetime, freq="D")
Flow = 5. 
Fhigh = 20. 

#   define the RShake station of interest 
network = "AM"
#station = "R9E71"   #  RHS Bridgewater
#name = "RHS Bridgewater, Worsley" 
#station = "R6C8A"   #  MMU Oxford Road
#name = "MMU Oxford Road"
#station = "RAEED"   #  UoM Oxford Road
#name = "UoM Oxford Road"
#station = "R36CB"   #  St James Cheadle Hulme
#name = "St James' CHS, Cheadle Hulme"
station = "R3FEA"   #  Oswald Road Chorlton 
name = "Oswald Road Primary, Chorlton"
location = "00"
channel = "EHZ"
client = Client('RASPISHAKE')
sFullName = network + "-" + station + "-" + location + "-" + channel 

#   get the waveforms from the station(s), as Obspy streams 
#   save each day to a separate (.mseed) file 
for day in datelist:
    
    fn = sFullName + "_" + day.strftime("%Y-%m-%d.mseed")
    print(fn)
    
    if day != datelist[-1] and os.path.isfile(fn):
        continue
    else:
        st = client.get_waveforms(network, station, location, channel,
                             UTCDateTime(day), UTCDateTime(day)+86400, attach_response=True)
        st.merge(method=0, fill_value='latest')
        print(st)
        st.write(fn)
        
resp = client.get_stations(UTCDateTime(day), network=network, station=station, location=location,
                      channel=channel, level="response")
print(resp)

#   combine streams in probabilistic power spectral density objects
#   and save each day's data in a separate .npz file    
for day in datelist:
    
    fn_in = sFullName + "_" + day.strftime("%Y-%m-%d.mseed")
    fn_out = sFullName + "_" + day.strftime("%Y-%m-%d.npz")
    
    if day != datelist[-1] and os.path.isfile(fn_out):
        continue
    else:    
        st = read(fn_in)
        st.attach_response(resp)
        print(st)
        ppsd = PPSD(st[0].stats, metadata=resp,
                    ppsd_length=600, overlap=0.5,
                    period_smoothing_width_octaves=0.025,
                    period_step_octaves=0.0125,
                    period_limits=(0.008, 50),
                    db_bins=(-200, 20, 0.25))
        ppsd.add(st)
        ppsd.save_npz(fn_out[:-4])
        del st, ppsd

#   Reload daily PSDs from the disk and create a single PPSD object
for day in datelist:
    fn = sFullName + "_" + day.strftime("%Y-%m-%d.npz")
    if day == datelist[0]:
        ppsd = PPSD.load_npz(fn)
    else:
        ppsd.add_npz(fn)
        
# Define frequency bands of interest:
freqs = [(0.1,1.0),(1.0,20.0),(4.0,14.0),(5.0,20.0),(5.0,40.0)]

def rms(s, f):
    # Parseval: the RMS in time domain is the sqrt of the integral of the power spectrum
    return np.sqrt(np.trapz(s, f))

per = ppsd.period_bin_centers
displacement_RMS = []
for psd in ppsd.psd_values:
    RMS = {}
    for fmin, fmax in freqs:
        ix = np.where((per>=1.0/fmax) & (per<=1.0/fmin))

        # acceleration power spectrum in Hz
        spec = psd.copy()[ix][::-1]
        f = 1.0/per.copy()[ix][::-1]

        # remove NaNs from the list
        valid = np.where(np.isfinite(spec))[0]
        spec = spec[valid]
        f = f[valid]

        w2f = (2.0 * np.pi * f)
        
        # The acceleration amplitude spectrum (dB to Power! = divide by 10 and not 20!)
        amp = 10.0**(spec/10.) 
        
        # velocity spectrum (divide by omega**2)
        vamp = amp / w2f**2
         
        # displacement spectrum (divide by omega**2)
        damp =  vamp / w2f**2

        RMS["%.1f-%.1f"%(fmin, fmax)] = rms(damp, f)

    displacement_RMS.append(RMS)

index = pd.DatetimeIndex([d.datetime for d in ppsd.times_processed])
displacement_RMS = pd.DataFrame(displacement_RMS, index=index)
band = "5.0-40.0"
d = displacement_RMS[band]
#pf = 12 
fig, axs = plt.subplots(1,1,figsize=(12,4))
#axs.plot(d.index[::pf], d[::pf]*1000, '-r')
axs.plot(d.index, d*1000, '-r')
print(len(d))
axs.set_ylabel("Displacement, mm")
axs.set_title('Averaged (RMS) displacement for %s - Filter: [%s] Hz' 
          % ("%s.%s.%s.%s" % (network, station, location, channel),
             band))
axs.set_ylim(0,0.00016)
axs.grid(True)
plt.autoscale(enable=True, axis='x', tight=True)
fig.autofmt_xdate()
fn = station + '_rmsD.png'
plt.savefig(fn, dpi=300)
plt.show()

from datetime import datetime
#   read in file of traffic counts  
#fnTC = 'Drakewell 1312 - 1 hour data - 27Apr22-10May22.csv' # R3FEA/Chorlton 
fnTC = 'Drakewell 1312 - 5 minute data - 28Apr22.csv' # R3FEA/Chorlton 
#fnTC = 'Drakewell 1414 - 1 hour data - 27Apr22-10May22.csv' # RAEED/UoM
#fnTC = 'Drakewell 1426 - 1 hour data - 27Apr22-10May22.csv' # R36CB/Cheadle Hulme 
#fnTC = 'Drakewell 1426 - 1 hour data - 28Apr22.csv' # R36CB/Cheadle Hulme 
#fnTC = 'Drakewell 1426 - 5 minute data - 28Apr22.csv' # R36CB/Cheadle Hulme 
#fnTC = 'Drakewell 1421 - 1 hour data - 27Apr22-10May22.csv' # R6C8A/MMU
#fnTC = 'Drakewell1305.csv' # R9E71/RHS
#fnTC = 'Drakewell 1414 - 1 hour data - 03May22-16May22.csv' # RAEED/UoM
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

fig, axs = plt.subplots(1,1,figsize=(12,4))
axs.plot(dtList, fHeavyCountList, label='Heavy vehicles') 
#axs.plot(dtList, fLightCountList, label='Light vehicles') 
#axs.plot(dtList, fTotalCountList, label='Total vehicles') 
axs.grid(True)
axs.set_ylabel('Number of vehicles')
axs.set_title('Traffic counts, TfGM Drakewell camera 1426')
axs.legend() 
plt.autoscale(enable=True, axis='x', tight=True)
fig.autofmt_xdate()
fn = station + '_Heavy_Drakewell1426_Traffic.png'
plt.savefig(fn, dpi=300)
plt.show()

fig, axs = plt.subplots(1,1,figsize=(12,4))
#axs.plot(dtList, fHeavyCountList, label='Heavy vehicles') 
axs.plot(dtList, fLightCountList, label='Light vehicles') 
#axs.plot(dtList, fTotalCountList, label='Total vehicles') 
axs.grid(True)
axs.set_ylabel('Number of vehicles')
axs.set_title('Traffic counts, TfGM Drakewell camera 1426')
axs.legend() 
plt.autoscale(enable=True, axis='x', tight=True)
fig.autofmt_xdate()
fn = station + '_Light_Drakewell1426_Traffic.png'
plt.savefig(fn, dpi=300)
plt.show()

fig, axs = plt.subplots(1,1,figsize=(12,4))
#axs.plot(dtList, fHeavyCountList, label='Heavy vehicles') 
#axs.plot(dtList, fLightCountList, label='Light vehicles') 
axs.plot(dtList, fTotalCountList, label='Total vehicles') 
axs.grid(True)
axs.set_ylabel('Number of vehicles')
axs.set_title('Traffic counts, TfGM Drakewell camera 1426')
axs.legend() 
plt.autoscale(enable=True, axis='x', tight=True)
fig.autofmt_xdate()
fn = station + '_All_Drakewell1426_Traffic.png'
plt.savefig(fn, dpi=300)
plt.show()
