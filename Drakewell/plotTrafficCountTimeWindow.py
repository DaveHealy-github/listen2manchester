#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:18:34 2021

@author: davidhealy
"""
import matplotlib.pyplot as plt
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

#   plot the data 
fig = plt.figure(figsize=(12,4))
#plt.plot(dtList, fTotalCountList, label='Total, all vehicles') 
plt.plot(dtList, fHeavyCountList, label='Heavy vehicles') 
plt.plot(dtList, fLightCountList, label='Light vehicles') 
plt.grid(True)
plt.xlabel('Date & time')
plt.ylabel('Number of vehicles')
plt.title('Traffic counts from TfGM Drakewell cameras')
plt.legend() 
plt.savefig('plotTrafficCountTimeWindow.png', dpi=300)
