#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 10:01:06 2021

@author: davidhealy
"""
import matplotlib.pyplot as plt
from datetime import datetime

#   read in file of fault dips, if required  
fnAQ = '2021-11-12-211118092446Edited.csv'
iLine = 0 
#   open the file 
with open(fnAQ, 'r') as reader:
    dtList = [] 
    fPM10List = []
    fNOList = [] 
    fNO2List = [] 
    fNOXasNO2List = [] 
    fO3List = [] 
    fPM25List = [] 
    for line in reader:
        iLine += 1 
        if iLine >= 7:
            sLineTokens = line.split(',')
            print(sLineTokens)
            for i, s in enumerate(sLineTokens):
                if i == 1 and len(s) > 0:
                    if "24:" in s[0:3]:
                        s = s.replace("24:", "00:")
                        sDate = sLineTokens[i-1]
                        day = int(sDate[0:2])
                        sLineTokens[i-1] = sDate.replace(str(day), str(day+1))
                    if "\n" in s[8:10]:
                        s = s.replace("\n", "")
                        sDateTime = sLineTokens[i-1] + ' ' + s
                        dDateTime = datetime.strptime(sDateTime, "%d/%m/%Y %H:%M:%S")
                        dtList.append(dDateTime)
                        fPM10List.append(float('nan'))
                        fNOList.append(float('nan'))
                        fNO2List.append(float('nan'))
                        fNOXasNO2List.append(float('nan'))
                        fO3List.append(float('nan'))
                        fPM25List.append(float('nan'))
                        continue
                    else:
                        sDateTime = sLineTokens[i-1] + ' ' + s
                        dDateTime = datetime.strptime(sDateTime, "%d/%m/%Y %H:%M:%S")
                        dtList.append(dDateTime)
                elif i == 2: 
                    if len(s) > 0:
                        fPM10List.append(float(s))
                    else:
                        fPM10List.append(float('nan'))
                elif i == 4:
                    if len(s) > 0:
                        fNOList.append(float(s))
                    else:
                        fNOList.append(float('nan'))
                elif i == 6:
                    if len(s) > 0:
                        fNO2List.append(float(s))
                    else:
                        fNO2List.append(float('nan'))
                elif i == 8:
                    if len(s) > 0:
                        fNOXasNO2List.append(float(s))
                    else:
                        fNOXasNO2List.append(float('nan'))
                elif i == 10: 
                    if len(s) > 0:
                        fO3List.append(float(s))
                    else:
                        fO3List.append(float('nan'))
                elif i == 12:
                    if len(s) > 0:
                        fPM25List.append(float(s))
                    else:
                        fPM25List.append(float('nan'))
                else:
                    continue 
    
print(iLine)

fig = plt.figure(figsize=(12,4))
plt.plot(dtList, fPM10List, label='PM$_{10}$') 
plt.plot(dtList, fPM25List, label='PM$_{2.5}$')
plt.plot(dtList, fNOList, label='NO') 
plt.plot(dtList, fNO2List, label='NO$_2$') 
#plt.plot(dtList, fNOXasNO2List, label='NO$_x$ as NO$_2$')
#plt.plot(dtList, fO3List, label='O$_3$')
plt.grid(True)
plt.xlabel('Date & time')
plt.ylabel('Concentration, $\mu$g $m^{-3}$')
plt.title('Air quality from DEFRA Salford M60')
plt.legend() 
plt.savefig('plotAirQualityTimeWindow.png', dpi=300)
