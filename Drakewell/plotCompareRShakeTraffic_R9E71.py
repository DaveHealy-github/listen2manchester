#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 11:18:40 2021

script to plot ambient noise (e.g., 4-14 Hz) from Raspberry Shakes 
over a specified time window

@author: davidhealy
    - cloned from seismosocialdistancing.ipynb by Thomas Lecoq 
"""

import os
import glob 

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
from obspy.signal import PPSD
from hampel import hampel

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

#   tidy up old work files
# for f in glob.glob("*.npz"):
#     os.remove(f)
# for f in glob.glob("*.mseed"):
#     os.remove(f)

#   define the time window of interest 
start = UTCDateTime("2021-11-14")
end = UTCDateTime("2021-11-27") # means "now"
datelist = pd.date_range(start.datetime, end.datetime, freq="D")

#   define the stations of interest 
network = "AM"
station = "R9E71"   #  RHS Bridgewater
location = "00"
channel = "EHZ"
sFullName = network + "-" + station + "-" + location + "-" + channel 

#   get the waveforms from the station(s), as Obspy streams 
#   save each day to a separate (.mseed) file 
c = Client(base_url='https://fdsnws.raspberryshakedata.com/')
for day in datelist:
    
    fn = sFullName + "_" + day.strftime("%Y-%m-%d.mseed")
    print(fn)
    
    if day != datelist[-1] and os.path.isfile(fn):
        continue
    else:
        st = c.get_waveforms(network, station, location, channel,
                             UTCDateTime(day)-1801, UTCDateTime(day)+86400+1801, attach_response=True)
        st.merge(method=0, fill_value='latest')
#        st.decimate(10)
        print(st)
        st.write(fn)
        
resp = c.get_stations(UTCDateTime(day), network=network, station=station, location=location,
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
                    ppsd_length=1800, overlap=0.5,
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
freqs = [(0.1,1.0),(1.0,20.0),(4.0,14.0),(4.0,20.0)]

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
band = "4.0-20.0"
d = displacement_RMS[band]
dser = pd.Series(d)
ts_imputation = hampel(dser, window_size=5, n=1, imputation=True) 
print(ts_imputation)
print(len(d))

fig, axs = plt.subplots(1,1,figsize=(15,5))
axs.plot(d.index, d, '-r')
axs.set_ylabel("Displacement, nm")
axs.set_title('Seismic noise for %s - Filter: [%s] Hz' 
          % ("%s.%s.%s.%s" % (network, station, location, channel),
             band))
#axs[0].set_xlim(d.index.min(), d.index.max())
axs.set_ylim(0,1.e-7)
axs.grid(True)
#axs[0].gca().set_axisbelow(True)
fig.autofmt_xdate()


fig, axs = plt.subplots(1,1,figsize=(15,5))
axs.plot(ts_imputation, '-b')
# Get normal business days and set their background color to green
# db = pd.bdate_range(start.datetime, end.datetime)
# for dbi in db:
#     plt.axvspan(dbi, dbi+datetime.timedelta(days=1),
#                 facecolor='lightgreen', edgecolor="none",
#                 alpha=0.2, zorder=-10)
scale = 1e9
#ticks = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x*scale))
#axs[0].gca().yaxis.set_major_formatter(ticks)
axs.set_ylabel("Displacement, nm")
axs.set_title('Seismic noise for %s - Filter: [%s] Hz' 
          % ("%s.%s.%s.%s" % (network, station, location, channel),
             band))
#axs[0].set_xlim(d.index.min(), d.index.max())
axs.set_ylim(0,0.6e-7)
axs.grid(True)
#axs[0].gca().set_axisbelow(True)
fig.autofmt_xdate()
plt.savefig('plotCompareRShakeTraffic_R9E71_SeismicRMSd.png', dpi=300)
plt.show()

from datetime import datetime
#   read in file of fault dips, if required  
fnTC = 'Drakewell1305.csv'
iLine = 0 
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
                                '%d-%m-%Y %H:%M:%S')
                    dtList.append(dDateTime)
                elif i <= 12: 
                    fCount += float(s)
                    if i >= 7 and i <= 12:
                        fHeavyCount += float(s) 
                else:
                    continue 
            fTotalCountList.append(fCount)
            fHeavyCountList.append(fHeavyCount)
            
print("Read in %i rows of data." % iLine)
fLightCountList = []
for i in range(len(fTotalCountList)):
    fLightCountList.append(fTotalCountList[i] - fHeavyCountList[i]) 

fig, axs = plt.subplots(1,1,figsize=(15,5))
axs.plot(dtList, fHeavyCountList, label='Heavy vehicles') 
axs.plot(dtList, fLightCountList, label='Light vehicles') 
axs.plot(dtList, fTotalCountList, label='Total vehicles') 
# for dbi in db:
#     plt.axvspan(dbi, dbi+datetime.timedelta(days=1),
#                 facecolor='lightgreen', edgecolor="none",
#                 alpha=0.2, zorder=-10)
axs.grid(True)
axs.set_xlabel('Date & time')
axs.set_ylabel('Number of vehicles')
axs.set_title('Traffic counts, TfGM Drakewell camera 1305')
axs.legend() 
fig.autofmt_xdate()

plt.savefig('plotCompareRShakeTraffic_R9E71_Traffic_Decimated.png', dpi=300)
plt.show()
