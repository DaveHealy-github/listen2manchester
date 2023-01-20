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
import matplotlib.pyplot as plt
#from pylab import *

# Start of parameters to define
NETWORK = "AM"   # AM = RaspberryShake network
STATION = "R72C7"  # Station code of station to get data for
CHANNEL = "EHZ"  # channel to grab data for (e.g. EHZ, SHZ, EHE, EHN)
# Idaho, 31-3-2020, M6.5 
EQ_TIME = "2021-10-02T06:29:18"  # origin time of earthquake
T_START = 0 # Length in seconds of data to plot before origin time
T_END = 60*1500  # Length in seconds of data to plot after origin time
PHASES = ["P", "S"]   # list of phases to compute theoretical times for
EVT_LAT = -21.104  # Latitude of event
EVT_LON = -174.895  # Longitude of event 
EVT_Z = 535.8  # Depth of event
STA_LAT = 51.33  # Latitude of station  
STA_LON = -0.49  # Longitude of station
F1 = 0.5  # High-pass filter corner
F2 = 1.5 # Low-pass filter corner 
LABEL = "M7.3 Vanuatu, 02 Oct 2021"  # Title to plot on figure
MODEL = 'iasp91'  # Velocity model to predict travel-times through
# End of parameters to define

# Define fdsn client to get data from
client = Client('http://fdsnws.raspberryshakedata.com')

# Define start and end time
orig_time = UTCDateTime(EQ_TIME)
t1 = orig_time - T_START
t2 = orig_time + T_END
# Download and filfter data
st = client.get_waveforms(NETWORK, STATION, "00", CHANNEL,
                          starttime=t1, endtime=t2, attach_response=True)
st.merge(method=0, fill_value='latest')
st.detrend(type="demean")
#st.remove_response()
st.filter("bandpass", freqmin=F1, freqmax=F2, corners=4)
st.trim(t1, t2)

# Set-up figure
fig = plt.figure(figsize=(16, 4))
ax = plt.gca()
#plt.suptitle(LABEL)

# Now plot the waveform data
ax.plot(st[0].times(reftime=orig_time), st[0].data*1000, linewidth=1,
        color="darkred")
ymin, ymax = ax.get_ylim()
ax.set_xlabel("Time from start (seconds)")
ax.set_title("{:} - {:}.{:}.{:}.{:}\nBandpass filter: {:}-{:} Hz".format(
    LABEL, st[0].stats.network, st[0].stats.station, st[0].stats.location,
    st[0].stats.channel, F1, F2))
ax.set_ylabel("Ground velocity (mm/s)")
plt.grid(True) 
plt.xlim(1000, 1500)
plt.ylim(-max(abs(ymin),abs(ymax)), max(abs(ymin),abs(ymax)))

# Save and plot the figure
plt.savefig("traces.png")
plt.show()
