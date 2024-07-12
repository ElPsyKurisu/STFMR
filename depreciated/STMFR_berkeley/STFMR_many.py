# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 09:39:43 2020

@author: Isaac
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


if len(sys.argv) > 1: filenames = sys.argv[1:]
else: filenames = ["Pt-CoFe-BPBO_021621\9GHz.lvm"]


#else: print("Please list one or more STFMR files for me to read")



data = []
len_files = len(filenames)
frequencies = []
fields = []
Vmixes =[]
labels = []
ii = 0
for file in filenames:
    with open(file, 'r') as file:
        dataframe = pd.read_csv(file, delimiter = '\t', skiprows = 5, 
                               skipfooter = 1,usecols = [0,2,3], engine = 'python')#for skipfooter
    
    points = dataframe.to_numpy()
    datum = np.transpose(points)
    
    #datum[0] = freq, datum[1] = field, datum[2] = Vmix
 
    frequencies.append(datum[0][0])
    fields.append(datum[1])
    Vmixes.append(datum[2])
    
  
    
    if filenames[ii][0:10] == "Pt-CoFe-Si": labels.append("") #reference")
    else: 
        print ii
        labels.append("%.1f GHz"% frequencies[ii])
    #labels.append("%.1f GHz"% frequencies[ii])
    
    ii += 1
fig, ax = plt.subplots()
for i in range(len_files):
    ax.plot(fields[i],Vmixes[i], '.', label = labels[i])
ax.ticklabel_format(axis='y', style = 'scientific', scilimits = (2,-2))
ax.set_xlabel("Field (Oe)")
ax.set_ylabel("V$_{mix}$ (V)")
ax.legend()

#ax.set_title()

plt.show()