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
STATION = "R72C7"  # Station code of station to get data for
CHANNEL = "EHZ"  # channel to grab data for (e.g. EHZ, SHZ, EHE, EHN)
EQ_TIME = "2020-02-14T13:00:00"  # origin time of earthquake
T_START = 0 # Length in seconds of data to plot before origin time
T_END = 35*60  # Length in seconds of data to plot after origin time
PHASES = ["P", "S"]   # list of phases to compute theoretical times for
EVT_LAT = 17.916  # Latitude of event
EVT_LON = -66.813  # Longitude of event 
EVT_Z = 10  # Depth of event
STA_LAT = 51.33  # Latitude of station  
STA_LON = -0.49  # Longitude of station
F1 = 0.3  # High-pass filter corner
F2 = 0.5  # Low-pass filter corner 
LABEL = "Raspberry Shake sample"  # Title to plot on figure
MODEL = 'iasp91'  # Velocity model to predict travel-times through
# End of parameters to define

# Define fdsn client to get data from
client = Client('http://fdsnws.raspberryshakedata.com')

# Define start and end time
orig_time = UTCDateTime(EQ_TIME)
t1 = orig_time - T_START
t2 = orig_time + T_END
# Download and filfter data
print(STATION)
st = client.get_waveforms(NETWORK, STATION, "00", CHANNEL,
                          starttime=t1, endtime=t2, attach_response=True)
st.merge()
st.detrend(type="demean")
st.remove_response()
#st.filter("bandpass", freqmin=F1, freqmax=F2, corners=4)
st.trim(t1, t2)

# Set-up figure
fig = plt.figure(figsize=(8, 4))
#plt.suptitle(LABEL)
ax = plt.gca()

# Now plot the waveform data
ax.plot(st[0].times(reftime=orig_time), st[0].data*1000, linewidth=1,
        color="darkred")
ymin, ymax = ax.get_ylim()
ax.set_xlabel("Time from start (s)")
ax.set_title("{:}.{:}.{:}.{:}\nBandpass filter: {:}-{:} Hz".format(
    st[0].stats.network, st[0].stats.station, st[0].stats.location,
    st[0].stats.channel, F1, F2))
ax.set_ylabel("Ground velocity (mm/s)")
plt.grid(True) 
plt.ylim(-max(abs(ymin),abs(ymax)), max(abs(ymin),abs(ymax)))
# Save and plot the figure
plt.savefig("traces.png")
plt.show()