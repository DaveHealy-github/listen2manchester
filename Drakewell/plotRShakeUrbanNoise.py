#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 09:46:18 2022

@author: davidhealy
    - cloned from seismosocialdistancing.ipynb by Thomas Lecoq 
"""

#   tools we need 
import os
import glob 
import obspy.signal 

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

from matplotlib.dates import DateFormatter
from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
from obspy.signal import PPSD
#from hampel import hampel

#   tidy up old work files
for f in glob.glob("*.npz"):
    os.remove(f)
for f in glob.glob("*.mseed"):
    os.remove(f)

#############################################################################
#   get a specified time window of RShake data from a station 
#############################################################################

#   define the time window of interest 
start = UTCDateTime("2022-04-27T00:00:00")
nDays = 14. 
twSpan = nDays*24.*60.*60.-1
end = start + twSpan
datelist = pd.date_range(start.datetime, end.datetime, freq="D")
Flow = 8. 
Fhigh = 20. 

#   define the RShake station of interest 
network = "AM"
#station = "R9E71"   #  RHS Bridgewater
#name = "RHS Bridgewater, Worsley" 
station = "R6C8A"   #  MMU Oxford Road
name = "MMU Oxford Road"
#station = "RAEED"   #  UoM Oxford Road
#name = "UoM Oxford Road"
#station = "R36CB"   #  St James Cheadle Hulme
#name = "St James' CHS, Cheadle Hulme"
#station = "R3FEA"   #  Oswald Road Chorlton 
#name = "Oswald Road Primary, Chorlton"
location = "00"
channel = "EHZ"
client = Client('RASPISHAKE')
sFullName = network + "-" + station + "-" + location + "-" + channel 

# # Download and filter data for station1
# st1 = client.get_waveforms(network, station, location, channel,
#                       starttime=start, endtime=end, attach_response=True)
# #st1.decimate(1)
# st1.merge(method=0, fill_value='latest')
# st1.detrend(type="demean")
# st1.remove_response()

# st1f = st1.copy()
# st1f.filter("bandpass", freqmin=Flow, freqmax=Fhigh, corners=4)

# # Envelope of filtered data
# st1_env = obspy.signal.filter.envelope(st1f[0].data*1000)

# # #   manual 'decimation' of streams 
# # iSample = 60000
# # nDec = int(len(st1f[0].data) / iSample)+1
# # iThis = 0  
# # st1_dec = np.zeros([nDec,1])
# # for i in range(len(st1f[0].data)):
# #     if i == 0:
# #         st1_dec[iThis] = st1f[0].data[i]
# #         iThis += 1
# #     elif i % iSample == 0:
# #         st1_dec[iThis] = st1f[0].data[i]
# #         iThis += 1
# #     else:
# #         continue 

# #   plot the velocity data as a function of time 
# fig = plt.figure(figsize=(6,2))
# plt.plot(st1f[0].times(reftime=start), st1f[0].data*1000, linewidth=1,
#         color="darkred")
# plt.grid(True)
# plt.xlabel("Time from start (s)")
# plt.title("{:}.{:}.{:}.{:} {:} - bandpass filter: {:}-{:} Hz".format(
#     st1[0].stats.network, st1[0].stats.station, st1[0].stats.location,
#     st1[0].stats.channel, name, Flow, Fhigh))
# plt.ylabel("Ground velocity (mm/s)")
# plt.xlim(0, twSpan)

# del st1, st1f 

#   plot a spectrogram (power spectral density) from 1-50 Hz

#   plot RMS displacement, after Lecocq et al. 

#   plot envelope of velocity, after Diaz et al. 
#   plot the velocity data as a function of time 
# fig = plt.figure(figsize=(20,5))
# #plt.plot(st1f[0].times(reftime=start), np.sqrt((st1f[0].data*1000)**2), linewidth=1,
# plt.plot(st1f[0].times(reftime=start), st1_env, linewidth=1,
#         color="darkblue")
# plt.grid(True)
# plt.xlabel("Time from start (s)")
# plt.title("{:}.{:}.{:}.{:} {:} - bandpass filter: {:}-{:} Hz".format(
#     st1[0].stats.network, st1[0].stats.station, st1[0].stats.location,
#     st1[0].stats.channel, name, Flow, Fhigh))
# plt.ylabel("Ground velocity envelope (mm/s)")
# plt.xlim(0, twSpan)
# plt.ylim(0, 0.8)

#   get the waveforms from the station(s), as Obspy streams 
#   save each day to a separate (.mseed) file 
for day in datelist:
    
    fn = sFullName + "_" + day.strftime("%Y-%m-%d.mseed")
    print(fn)
    
    if day != datelist[-1] and os.path.isfile(fn):
        continue
    else:
        st = client.get_waveforms(network, station, location, channel,
                             UTCDateTime(day)-1801, UTCDateTime(day)+86400+1801, attach_response=True)
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
freqs = [(0.1,1.0),(1.0,20.0),(4.0,14.0),(8.0,20.0),(10.,40.)]

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
band = "8.0-20.0"
d = displacement_RMS[band]
# dser = pd.Series(d)
# ts_imputation = hampel(dser, window_size=5, n=1, imputation=True) 
# print(ts_imputation)
#print(len(d))

fig, axs = plt.subplots(1,1,figsize=(12,4))
axs.plot(d.index, d*1000, '-r')
axs.set_ylabel("Displacement, mm")
axs.set_title('Seismic noise (rms displacement) for %s - Filter: [%s] Hz' 
          % ("%s.%s.%s.%s" % (network, station, location, channel),
             band))
#axs[0].set_xlim(d.index.min(), d.index.max())
#axs.set_ylim(0,0.00006)
axs.grid(True)
#axs[0].gca().set_axisbelow(True)
plt.autoscale(enable=True, axis='x', tight=True)
fig.autofmt_xdate()
fn = station + '_rmsD.png'
plt.savefig(fn, dpi=600)

#fig, axs = plt.subplots(1,1,figsize=(10,10))
fig = ppsd.plot_spectrogram(cmap='rainbow', clim=[-100,-70], show=False)
fig.axes
ax = fig.axes[0]
plt.gca().invert_yaxis()
ax.set_ylabel('Frequency, Hz')
plt.yticks([0.2,0.1,0.05,0.033,0.025], ['5','10','20','30','40'])
ax.set_ylim(0.2,0.025)
ax.set_title('Seismic noise for %s' 
          % ("%s.%s.%s.%s" % (network, station, location, channel)))
plt.gcf().set_size_inches(13, 5.)
plt.gcf().autofmt_xdate()

fn = station + '_PSD.png'
plt.savefig(fn, dpi=600)
plt.show()

from datetime import datetime
#   read in file of fault dips, if required  
#fnTC = 'Drakewell 1312 - 1 hour data - 27Apr22-10May22.csv'
fnTC = 'Drakewell 1414 - 1 hour data - 27Apr22-10May22.csv'
#fnTC = 'Drakewell 1426 - 1 hour data - 27Apr22-10May22.csv'
#fnTC = 'Drakewell 1421 - 1 hour data - 27Apr22-10May22.csv'
#fnTC = 'Drakewell1305.csv'
#fnTC = 'Drakewell 1414 - 1 hour data - 03May22-16May22.csv'
iLine = 0 
iColHeavy = 7
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
axs.plot(dtList, fLightCountList, label='Light vehicles') 
axs.plot(dtList, fTotalCountList, label='Total vehicles') 
# for dbi in db:
#     plt.axvspan(dbi, dbi+datetime.timedelta(days=1),
#                 facecolor='lightgreen', edgecolor="none",
#                 alpha=0.2, zorder=-10)
axs.grid(True)
axs.set_ylabel('Number of vehicles')
axs.set_title('Traffic counts, TfGM Drakewell camera 1414')
axs.legend() 
plt.autoscale(enable=True, axis='x', tight=True)
fig.autofmt_xdate()

fn = 'Drakewell1414_Traffic.png'
plt.savefig(fn, dpi=600)
plt.show()
