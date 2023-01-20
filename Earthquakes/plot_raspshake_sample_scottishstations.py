#!/usr/bin/env python
"""
Script to download and plot RaspberryShake station data
Also computes and plots theoretical phase arrival times and raypaths
Stephen Hicks
Imperial College London
Feb 2020
Modified by:
    Dave Healy
    Feb 2020 
"""

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
import matplotlib.pyplot as plt

# Start of parameters to define
NETWORK = "AM"   # AM = RaspberryShake network
#STATIONS = ["R72C7", "R92E0", "R98E6"]  # Station code of station to get data for
STATIONS = ["R92E0", "R98E6"]  # Station code of station to get data for
CHANNEL = "EHZ"  # channel to grab data for (e.g. EHZ, SHZ, EHE, EHN)
T_START = 0 # Length in seconds of data to plot before origin time
T_END = 1000  # Length in seconds of data to plot after origin time
PHASES = ["P", "S"]   # list of phases to compute theoretical times for

#   M6.9 Kuril Islands, 2020-02-13T10:33:44
#LABEL = "M6.9 Kuril Islands on Scottish Raspberry Shakes"  # Title to plot on figure
#EQ_TIME = "2020-02-13T10:33:44"  # origin time of earthquake
#EVT_LAT = 45.67  # Latitude of event
#EVT_LON = -148.89  # Longitude of event 
#EVT_Z = 150  # Depth of event, km

#   M7.7 Jamaica, 2020-01-28T19:10:24
#LABEL = "M7.7 Jamaica on Scottish Raspberry Shakes"  # Title to plot on figure
#EQ_TIME = "2020-01-28T19:10:24"  # origin time of earthquake
#EVT_LAT = 19.421  # Latitude of event
#EVT_LON = 78.763  # Longitude of event 
#EVT_Z = 14.8  # Depth of event, km

#   M5.6 Turkey, 2020-01-24T17:55:14
LABEL = "M5.6 Turkey on Scottish Raspberry Shakes"  # Title to plot on figure
EQ_TIME = "2020-01-22T19:22:16"  # origin time of earthquake
EVT_LAT = 39.076  # Latitude of event
EVT_LON = -27.843  # Longitude of event 
EVT_Z = 8.8  # Depth of event, km

#   station coordinates 
#STA_LAT = [57.16, 56.94, 57.72]  # Latitude of station  
#STA_LON = [-2.09, -2.26, -5.69]  # Longitude of station
STA_LAT = [56.94, 57.72]  # Latitude of station  
STA_LON = [-2.26, -5.69]  # Longitude of station

F1 = 0.4  # High-pass filter corner
F2 = 0.7  # Low-pass filter corner 
MODEL = 'iasp91'  # Velocity model to predict travel-times through

# End of parameters to define

# Set-up figure
fig = plt.figure(figsize=(8, 10))
plt.suptitle(LABEL)


# Define fdsn client to get data from
client = Client('http://fdsnws.raspberryshakedata.com')

# Define start and end time
orig_time = UTCDateTime(EQ_TIME)
t1 = orig_time - T_START
t2 = orig_time + T_END

# loop through station list 
s = 0 
for statThis in STATIONS:
    s = s + 1 
#    statThis = [statThis] 
    print(statThis)
    # Download and filter data
    st = client.get_waveforms(NETWORK, statThis, "00", CHANNEL,
                          starttime=t1, endtime=t2, attach_response=True)
    st.merge()
    st.detrend(type="demean")
    st.remove_response()
    st.filter("bandpass", freqmin=F1, freqmax=F2, corners=4)
    st.trim(t1, t2)

    # Now plot the waveform data
    ax = plt.subplot(4,1,s)
    
    ax.plot(st[0].times(reftime=orig_time), st[0].data*1000, linewidth=1,
            color="darkred")
    
    # Set-up taup travel-time model
    dist = locations2degrees(EVT_LAT, EVT_LON, STA_LAT[s-1], STA_LON[s-1])
    model = TauPyModel(model=MODEL)    

    ymin, ymax = ax.get_ylim()
    
    for phase in PHASES:
        phase = [phase]
        tt = model.get_travel_times(source_depth_in_km=EVT_Z,
                                    distance_in_degree=dist,
                                    phase_list=phase)
        ax.vlines(tt[0].time, ymin, ymax, color="blue",
                  linewidth=1.2, zorder=3, linestyle="--", alpha=0.5)
        ax.text(tt[0].time*1.02, ymax, phase[0], fontsize=12,
                horizontalalignment="left", verticalalignment="top")
    
    ax.set_xlabel("Time from start (s)")
    ax.set_title("{:}.{:}.{:}.{:}\nBandpass filter: {:}-{:} Hz".format(
            st[0].stats.network, st[0].stats.station, st[0].stats.location,
            st[0].stats.channel, F1, F2))
    ax.set_ylabel("Ground velocity (mm/s)")
    plt.grid(True) 
    plt.ylim(-max(abs(ymin),abs(ymax)), max(abs(ymin),abs(ymax)))
    plt.xlim(250, 750)
    
# Save and plot the figure
plt.tight_layout(pad=3.5)
plt.savefig("traces_scottish_raspbs.png")
plt.show()