#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 10:04:13 2022

@author: davidhealy

plot PSD for separate Raspberry Shakes for a defined time interval
plot mean of PSD for each station/time  
"""
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.signal import PPSD, spectral_estimation

T_START = 0     # length in seconds of data to plot before origin time
T_END = 60*60*24     # length in seconds of data to plot after origin time
Fhigh = 1.   

#   now get some station meta data, esp. lat & long 
namStationList = ['R4DD7', 'R3FEA', 'R6C8A', 'R36CB', 'LBWR']
labelStationList = ['Firs', 'Chorlton', 'MMU', 'Cheadle Hulme', 'Ladybower']
netStation = 'AM'
locStation = '00'
chaStation = 'EHZ'      #   vertical component of geophone 
client = Client('RASPISHAKE')

#   define the time window
EQ_TIME = "2022-06-26T00:00:00"
dtEvent = UTCDateTime(EQ_TIME)
t1start = dtEvent + T_START
t1end = dtEvent + T_END - 1 
print(t1start, t1end)

#   loop through the stations 
print("Getting streams for each station...")
for namStation in namStationList:

    if "LBWR" in namStation:
        client=Client('http://eida.bgs.ac.uk')
        netStation = 'GB'
        locStation = '00' 
        chaStation = 'HHZ'
        
    #   get the response for this station
    resp = client.get_stations(t1start, network=netStation, station=namStation, location=locStation,
                      channel=chaStation, level="response")
    # Download and filter data for station1
    st = client.get_waveforms(netStation, namStation, locStation, chaStation,
                          starttime=t1start, endtime=t1end, attach_response=True)
    # st.merge(method=0, fill_value='latest')
    # st.detrend(type="demean")
    # st.remove_response()
    stf = st.copy()
    stf.filter("highpass", freq=Fhigh, corners=2)
    # stf.trim(t1start, t1end)
    
    #   build a new PPSD object for this station, this interval 
    ppsd = PPSD(stf[0].stats, resp,
                ppsd_length=1800, overlap=0.5,
                period_smoothing_width_octaves=0.025,
                period_step_octaves=0.0125,
                period_limits=(0.025, 50),
                db_bins=(-200, 20, 0.25))
    ppsd.add(stf)

    #   save this PPSD to a numpy binary file, *.npz
    fnout = "RS_PSD_" + namStation
    ppsd.save_npz(fnout)
    del st, stf, ppsd, resp 
        
#   loop through the separate *.npz files, merging into one PPSD object 
print("Merging streams for all stations...")
ppsds = {}
for namStation in namStationList:

    fnin = "RS_PSD_" + namStation + ".npz"
    mseedid = namStation 
    ppsds[mseedid] = PPSD.load_npz(fnin) 
    ppsds[mseedid].add_npz(fnin)
        
print(ppsds.items())

#   plot the PSDs as separate lines, as per Holmgren & Werner, 2021 Figure 2b
fig, ax = plt.subplots(figsize=(5,4)) 
iSta = 0
for mseedid, ppsd in ppsds.items():
    [period, meanpsd] = ppsd.get_mean()
    ax.plot(1./period, meanpsd, label=namStationList[iSta])
    iSta = iSta + 1 
ax.plot(1./spectral_estimation.get_nlnm()[0], spectral_estimation.get_nlnm()[1], color='gray')
ax.plot(1./spectral_estimation.get_nhnm()[0], spectral_estimation.get_nhnm()[1], color='gray')
ax.axis([0.5, 50., -220, -40])
ax.set_xscale('log')
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Amplitude ($m^2$/$s^4$/Hz)(dB)")
ax.text(0.95,-108,"NHMN", ha='right', color='gray')
ax.text(0.95,-175,"NLMN", ha='right', color='gray')
ax.grid(True)
ax.legend(loc='lower right')
plt.xticks([1., 2., 5., 10., 20., 50.], 
           ['1', '2', '5', '10', '20', '50'])
plt.tight_layout()
outfile='RS_PSD_mean.png'
plt.savefig(outfile, dpi=300)
