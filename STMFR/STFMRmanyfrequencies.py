# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 09:39:43 2020

@author: Isaac
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

if len(sys.argv) > 1: 
    filenames = sys.argv[1:-3]
    NM_thick_s = sys.argv[-3] #nm         string versions
    FM_thick_s = sys.argv[-2] #nm
    Ms_s = sys.argv[-1]       #emu/cc
    print("NM thickness = %s nm, FM thickness = %s nm, and Ms = %s emu/cc" %(NM_thick_s, FM_thick_s, Ms_s))
    NM_thick = np.float(NM_thick_s)
    FM_thick = np.float(FM_thick_s)
    Ms = np.float(Ms_s)                 #float versions
    
    
else: 
    filenames = ["Pt-CoFe-BPBO_021621\dev2_031121\\6GHz_r.lvm",
                 "Pt-CoFe-BPBO_021621\dev2_031121\\7GHz_r.lvm",
                 "Pt-CoFe-BPBO_021621\dev2_031121\\8GHz_r.lvm",
                 "Pt-CoFe-BPBO_021621\dev2_031121\\9GHz_r.lvm",
                 "Pt-CoFe-BPBO_021621\dev2_031121\\10GHz_r.lvm",
                 "Pt-CoFe-BPBO_021621\dev2_031121\\11GHz_r.lvm"]
    NM_thick = 50. #nm
    FM_thick = 4.13 #nm
    Ms = 800 #emu/cc
    print("NM thickness = %s nm, FM thickness = %s nm, and Ms = %s emu/cc" %(NM_thick, FM_thick, Ms))
    

#......ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo......
#Define important stuff
num_files = len(filenames)

H0_list = []
S_list = []
A_list = []
W_list = []
Offset_list = []

frequencies = [] 
poling_fields = []
fields = []
Vmixes = []
labels = []

prelim_fit_ranges = []
peak_indexes = []
fit_ranges = []

SHA_list = []

def stfmr_fit(x,S,A,W,H0,Offset):
    return S*(W**2/(W**2+(x-H0)**2))+A*((x-H0)*W/(W**2+(x-H0)**2)) + Offset
#symmetric (S) + antisymetric (A) lorentzian with center H0, width W, with Offset

def symmetric(x,S,W,H0,Offset):
    return S*(W**2/(W**2+(x-H0)**2)) + Offset

def antisymmetric(x,A,W,H0,Offset):
    return A*((x-H0)*W/(W**2+(x-H0)**2)) + Offset

g_e = 2.0 #lande g factor for electron
bohr_mag = 9.274e-21 #erg/G, CGS units
hbar = 1.055e-27 # erg*s, CGS units
gyromagnetic_ratio = g_e * bohr_mag / hbar # = 1.76e7 1/Gs
#then the prefactor gyromagnetic_ratio/2pi = 2.8e-3 GHz/G
mu_naught = 4*np.pi #CGS (emu) units
e = 1.602e-20 #abC, or Abcoulomb. CGS (emu) units.


def Kittel_fit(Bres,Meff): #working in CGS (emu) units where mu_naught = 4pi
    return 1.0e-9*(gyromagnetic_ratio/(2*np.pi))*np.sqrt((Bres)*(Bres + mu_naught*Meff)) #prefactor of 1e-9 to convert to GHz

def Calculate_SHA(S,A,Ms,t,d,Meff,Bres):
    if Bres > 0: return (S/A)*(e*mu_naught*Ms*t*d/hbar)*np.sqrt(1+mu_naught*Meff/Bres) #make sure everything is in CGS units, so t and d are in cm...
    if Bres < 0: return -1*(S/A)*(e*mu_naught*Ms*t*d/hbar)*np.sqrt(1+mu_naught*Meff/(-1*Bres))
    #I differentiate so we can compare the two sides better
        

#......ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo......
#Collect and label data from files
units_bool = True
for file in filenames:
    with open(file, 'r') as file:
        dataframe = pd.read_csv(file, delimiter = '\t', skiprows = 5, 
                               skipfooter = 1,usecols = [0,1,2,3], engine = 'python')#for skipfooter
    
    points = dataframe.to_numpy()
    datum = np.transpose(points)
    
    #datum[0] = freq, datum[1] = poling field, datum[2] = field, datum[3] = Vmix
    frequencies.append(datum[0][0]) #in GHz, careful
    poling_fields.append(datum[1][0])
    fields.append(datum[2])
    Vmixes.append(datum[3])
    
    
    if units_bool: #only give the GHz label to the first item in the legend
        labels.append("%i GHz"% np.float(datum[0][0]))
        units_bool = False
    else: labels.append("%i"% np.float(datum[0][0]))

#......ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo......
#Find peaks (simple min or max) and make a little range around the peak
#fig, ax = plt.subplots() #in case we want to plot the plus minus 10 points...

for ii in range(num_files):
    Peak_Vmix_index = np.argmin(Vmixes[ii]) #try argmin to fit the other side
    
    #ax.plot(fields[ii],Vmixes[ii], '.', label = labels[ii])
    B0_Sign = -1.0 #Fit to the positive (1) field resonance peaks or negative (-1) field peaks?
    #Hardcoding File peak positions if necessary, for some reason the if filenames thing wasn't working.
    if ii == 0: Peak_Vmix_index = np.argmin([np.abs(field - B0_Sign*500) for field in fields[ii]])  
    if ii == 1: Peak_Vmix_index = np.argmin([np.abs(field - B0_Sign*650) for field in fields[ii]])  
    if ii == 2: Peak_Vmix_index = np.argmin([np.abs(field - B0_Sign*800) for field in fields[ii]])  
    if ii == 3: Peak_Vmix_index = np.argmin([np.abs(field - B0_Sign*1000) for field in fields[ii]])  
    if ii == 4: Peak_Vmix_index = np.argmin([np.abs(field - B0_Sign*1200) for field in fields[ii]])  
    if ii == 5: Peak_Vmix_index = np.argmin([np.abs(field - B0_Sign*1500) for field in fields[ii]])  
    
    rangeAround_Peak_index = np.arange(Peak_Vmix_index-10,Peak_Vmix_index+10)  #make a little range around the peak for fitting
    #ax.plot(fields[ii][Peak_Vmix_index],Vmixes[ii][Peak_Vmix_index], 'ks', markersize = 4) #plot it to make sure we understand what we are looking at
    #ax.plot(fields[ii][rangeAround_Peak_index], Vmixes[ii][rangeAround_Peak_index], 'k.', markersize = 3)
    
    prelim_fit_ranges.append(rangeAround_Peak_index)
    peak_indexes.append(Peak_Vmix_index)
    

    
#ax.ticklabel_format(axis='y', style = 'scientific', scilimits = (2,-2))
#ax.set_xlabel("Field (Oe)")
#ax.set_ylabel("V$_{mix}$ (V)")
#ax.legend()

#......ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo......
# Fit to the small ranges around the peaks to get some fit values, then do it again to within twice the width


for ii in range(num_files):
    print("Fitting File %s" %labels[ii])
    prelim_popt, prelim_pcov = curve_fit(stfmr_fit, fields[ii][prelim_fit_ranges[ii]], Vmixes[ii][prelim_fit_ranges[ii]], 
                           p0 = [Vmixes[ii][peak_indexes[ii]], Vmixes[ii][peak_indexes[ii]], #fit parameters S, A
                                 np.abs(fields[ii][prelim_fit_ranges[ii][-1]]-fields[ii][prelim_fit_ranges[ii][0]]), #W
                                 fields[ii][peak_indexes[ii]], 0]) #H0, Offset

    extent = np.int(1.2*np.abs(prelim_popt[2])/np.abs(fields[ii][1]-fields[ii][0]))  #setting size of the final fit, roundingto nearest integer
    B0_index = np.argmin([np.abs(field - prelim_popt[3]) for field in fields[ii]])   #Finding true center of the fit
    fit_range_max = B0_index + extent
    fit_range_min = B0_index - extent
    if fit_range_max > len(fields[ii]): fit_range_max = len(fields[ii]) #making sure we only fit to data we have...
    if fit_range_min < 0: fit_range_min = 0
    
   
 
    
    fit_ranges.append(np.arange(fit_range_min, fit_range_max))

    popt, pcov = curve_fit(stfmr_fit, fields[ii][fit_ranges[ii]], Vmixes[ii][fit_ranges[ii]], p0 = prelim_popt)    
    
    S_list.append(popt[0])
    A_list.append(popt[1])
    W_list.append(popt[2])
    H0_list.append(popt[3])
    Offset_list.append(popt[4])



fig1, ax1 = plt.subplots(figsize = (4,5), dpi = 200) 


for ii in range(num_files):   
    H_lin = np.linspace(fields[ii][fit_ranges[ii][0]], fields[ii][fit_ranges[ii][-1]], 1000) #fields to show stfmr_fit function
    
    ax1.plot(fields[ii],Vmixes[ii], '.', label = labels[ii])
    
    ax1.plot(H_lin,stfmr_fit(H_lin, S_list[ii], A_list[ii], W_list[ii], H0_list[ii], Offset_list[ii]),
             'k-') #plot it to make sure the fits worked reasonably well
    
    #Plotting the components is optional, and can get cluttered
    #ax1.plot(H_lin,symmetric(H_lin, S_list[ii], W_list[ii], H0_list[ii], Offset_list[ii]),color = 'gray', linestyle = '-')
    #ax1.plot(H_lin,antisymmetric(H_lin, A_list[ii], W_list[ii], H0_list[ii], Offset_list[ii]),color = 'gray', linestyle = '-')
    
    prelim_fit_ranges.append(rangeAround_Peak_index)
    peak_indexes.append(Peak_Vmix_index)
    
ax1.ticklabel_format(axis='y', style = 'scientific', scilimits = (2,-2))
ax1.set_xlabel("Field (Oe)")
ax1.set_ylabel("V$_{mix}$ (V)")
ax1.legend() 

    
#......ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo......
#Fit the Kittel equation to find anisotropy field Hanis and effective magnetization Meff


Kittel_popt, Kittel_pcov = curve_fit(Kittel_fit, np.abs(H0_list), frequencies, p0 = [Ms/2]) #Meff
Meff = Kittel_popt[0]


fig2, ax2 = plt.subplots(figsize = (3,3), dpi = 200)

ax2.plot(np.abs(H0_list),frequencies,'sk', label = 'data')

H0s = np.linspace(0,np.max(np.abs(H0_list)),1000)
ax2.plot(H0s,Kittel_fit(np.abs(H0s),Meff), '-b',label = 'Kittel fit')
ax2.set_xlabel(r'$H_0$ (Oe)')
ax2.set_ylabel("Frequency (GHz)")
ax2.legend()

print("M_eff = %.2e" % (Meff))

#......ooooooOOOOOOoooooo............ooooooOOOOOOoooooo............ooooooOOOOOOoooooo......
#Home stretch: calculate and plot the spin hall angle as a function of frequency!!

fig3, ax3 = plt.subplots(figsize = (3,3), dpi = 200)
for ii in range(num_files):
    SHA_list.append(Calculate_SHA(S_list[ii], A_list[ii], Ms, NM_thick*1e-7, FM_thick*1e-7, Meff, H0_list[ii]))
ax3.plot(frequencies,SHA_list,'ks')
ax3.set_xlabel("Frequency (GHz)")
ax3.set_ylabel(r"$\theta_{SH}$")
#ax3.set_yscale("log")

print("SHA: %.2f +/- %.2f" %(np.mean(SHA_list), np.std(SHA_list)))


plt.show()



