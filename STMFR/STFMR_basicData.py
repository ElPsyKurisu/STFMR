# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:20:50 2021

@author: Isaac
"""

"""
Created by Geo Fratian - UC Berkeley 8/6/2021 5:11 PM

STFMR analysis Program
"""

"""
Program Overview:
Using txt files from stfmr and a google doc to grab sample info
In case no google doc is created can alternatively hardcode sample information - 
namely Ms, Magnet Thickness, and Conductor Thickness.
Note - DO NOT USE RUN BUTTON, INSTEAD USE OPEN IN INTEGRATED TERMINAL
python STFMR_basicData.py folderName\devname
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import scipy.signal as sci
from scipy.optimize import curve_fit
#from gdstfmr import* #get rid of this line if there is no internet

labels = []
argLabel_bool = False
if len(sys.argv) < 2: directories = [os.getcwd() + "\\1116b\\dev12"]
else: 
    arguments = sys.argv[1:]
    if 'labels' in arguments:
        print("Using labels")
        labels = arguments[arguments.index('labels')+1:]
        print(labels)
        argLabel_bool = True
        nonlabel_arguments = arguments[:arguments.index('labels')]
        directories = nonlabel_arguments
    else: directories = arguments


num_dir = len(directories)
fields = [] #nested list of all the field data for each frequency
vmixes = [] #nested list of all the vmix data for each frequency
files = []
numbers_of_files = []
for i in range(num_dir):
    directory = directories[i]
    temp_file_list = []
    for foo in os.listdir(directory):
        if foo.endswith(".lvm"): temp_file_list.append(foo)
    files.append(temp_file_list)
    if not argLabel_bool: 
        for f in files[i]: 
            labels.append(f[:-4])
    numbers_of_files.append(len(files[i]))
    for j in range(numbers_of_files[i]): #open each file in freq_dir
        with open(directory + '\\' + files[i][j], 'r') as file:
            lines = file.readlines()     
        #data starts from line 6
        len_of_data = len(lines) - 6 
        #setup individual files to be added to fields/vmixes
        field=[]
        vmix = []
        print("reading file %s" % files[i][j])
        for q in range(len_of_data): #build fields and vmixes
            temp = np.fromstring(lines[q+6], dtype=float, sep='\t')
            field.append(temp[0])
            vmix.append(temp[1])
        fields.append(field)
        vmixes.append(vmix)

#DATA COLLECTION FROM STMFR FILES IS COMPLETE

#TIME TO PLOT DATA/FITS
#FIG1 PLOTS THE STMFR FIT 
fig1, ax1 = plt.subplots() #(figsize = (3,3), dpi = 200) 
for ii in range(len(fields)):
    #H_lin = np.linspace(fields[ii][fit_ranges[ii][0]], fields[ii][fit_ranges[ii][-1]], 1000) #fields to show stfmr_fit function
    
    ax1.plot(fields[ii],vmixes[ii], '-', label = labels[ii], )
    
    #ax1.plot(H_lin,stfmr_fit(H_lin, S_list[ii], A_list[ii], W_list[ii], H0_list[ii], Offset_list[ii]),
    #         'k-') #plot it to make sure the fits worked reasonably well
    
    #Plotting the components is optional, and can get cluttered
    #ax1.plot(H_lin,symmetric(H_lin, S_list[ii], W_list[ii], H0_list[ii], Offset_list[ii]),color = 'gray', linestyle = '-')
    #ax1.plot(H_lin,antisymmetric(H_lin, A_list[ii], W_list[ii], H0_list[ii], Offset_list[ii]),color = 'gray', linestyle = '-')    
ax1.ticklabel_format(axis='y', style = 'scientific', scilimits = (2,-2))
ax1.set_xlabel("Field (Oe)")
ax1.set_ylabel("V$_{mix}$ (V)")
ax1.legend() 
plt.show()

"""
#GRAB DATA FROM GOOGLE SHEET FOR USE LATER
sample_df = data #pandas df with sample info
Mag_thick_lst = []
Con_thick_lst = []
Ms_lst = []
for b in range(num_dir):
    sample_name = directories[b].split('\\')[0]
    sample_info = sample_df.loc[sample_name]
    Mag_thick, Con_thick, Ms = sample_info[0], sample_info[1], sample_info[2]
    if Ms == '':
        print('Please Input Ms Into Google Doc for %s' %(sample_name))
        sys.exit()
    Mag_thick_lst.append(float(Mag_thick)), Con_thick_lst.append(float(Con_thick)), Ms_lst.append(float(Ms))
#DATA COLLECTION FROM GOOGLE SHEET IS COMPLETE

#DATA ANALYSIS BEGINS HERE
#Initialize all variables
H0_list = [] #center list
S_list = [] #symmetric list
A_list = [] #Anti-symmetric list
W_list = [] #width list
Offset_list = []
frequencies = []
labels = []
for mm in range(number_of_files): #build labels
    labels.append(freq_dir[mm][:-4])
frequencies_temp = []
counter = 0 #to count each char in filename
for gg in range(number_of_files):
    uncleanfreq = freq_dir[gg] #grab individual filename
    for char in uncleanfreq:
        if char == 'G':
            frequencies_temp = uncleanfreq[:counter]
            break
        else:
            counter += 1
    frequencies.append(float(frequencies_temp))
    counter = 0

#DEFINE FUNCTIONS FOR USE IN DATA ANALYSIS
def stfmr_fit(x,S,A,W,H0,Offset):
    return S*(W**2/(W**2+(x-H0)**2))+A*((x-H0)*W/(W**2+(x-H0)**2)) + Offset
#symmetric (S) + antisymetric (A) lorentzian with center H0, width W, with Offset

def symmetric(x,S,W,H0,Offset):
    return S*(W**2/(W**2+(x-H0)**2)) + Offset

def antisymmetric(x,A,W,H0,Offset):
    return A*((x-H0)*W/(W**2+(x-H0)**2)) + Offset

#CONSTANTS FOR USE IN FITTING FUNCTIONS
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

#NOW WE WANT TO TRY AND FIND THE PEAK POSITIONS
#First smooth data with savgol filter
vmixes_smooth = []
for ll in range(number_of_files):
    temp = sci.savgol_filter(vmixes[ll], 51, 10)
    vmixes_smooth.append(temp)

#use find_peaks from scipy to get index of peaks
peaks_lst = []
troughs_lst = []
for rr in range(number_of_files):
    peaks, _ = sci.find_peaks(vmixes_smooth[rr], width=8)
    troughs, _ = sci.find_peaks(-1*vmixes_smooth[rr], width=8)
    peaks_lst.append(peaks)
    troughs_lst.append(troughs)

#Now we want to find the correct peaks (max) and troughs (min)
correct_peaks = []
correct_troughs = []
for mm in range(number_of_files):
    vmix_smooth = vmixes_smooth[mm]
    correct_peaks.append(np.argmax(vmix_smooth[peaks_lst[mm]]))
    correct_troughs.append(np.argmin(vmix_smooth[troughs_lst[mm]]))

#NOW WE NEED TO MAKE A RANGE AROUND EACH CORRECT PEAK/TROUGH
peak_ranges = []
trough_ranges=[]
peak_indexes = []
trough_indexes = []
for mm in range(number_of_files):
    peakvmixindex = peaks_lst[mm][correct_peaks[mm]]
    troughvmixindex = troughs_lst[mm][correct_troughs[mm]]
    range_around_peak_vmixindex = np.arange(peakvmixindex-15, peakvmixindex+15)
    range_around_trough_vmixindex = np.arange(troughvmixindex-15, troughvmixindex+15)
    trough_ranges.append(range_around_trough_vmixindex)
    peak_ranges.append(range_around_peak_vmixindex)
    peak_indexes.append(peakvmixindex)
    trough_indexes.append(troughvmixindex) 

#BOOL TO ALLOW FOR EITHER FIT AROUND TROUGHS OR PEAK
troughs_bool = False

#TIME TO FIT THE DATA USING THE FUNCTIONS WE DEFINED EARLIER
fit_ranges = []
for ii in range(number_of_files):
    print("Fitting File %s" %labels[ii])
    #put stuff to plot here to simplify the curve fit
    if troughs_bool:
        temp = np.array(trough_ranges[ii])
    else:
        temp = np.array(peak_ranges[ii])
    field = fields[ii]
    field_arr = np.array(field)
    fields_to_fit = np.array(field_arr[temp])
    vmix = vmixes[ii]#changed to vmix not vmixes_smooth
    vmix_arr = np.array(vmix)
    vmix_to_fit = np.array(vmix_arr[temp])

    #DEPENDING ON THE BOOL, FITS STMFR TO TROUGH VERSUS PEAK FIT

    if troughs_bool:
        prelim_popt, prelim_pcov = curve_fit(stfmr_fit, fields_to_fit, vmix_to_fit, 
                           p0 = [vmixes[ii][trough_indexes[ii]], vmixes[ii][trough_indexes[ii]], #fit parameters S, A
                                 np.abs(fields[ii][trough_ranges[ii][-1]]-fields[ii][trough_ranges[ii][0]]), #W
                                 fields[ii][trough_indexes[ii]], 0], maxfev=5000) #H0, Offset

    else:
        prelim_popt, prelim_pcov = curve_fit(stfmr_fit, fields_to_fit, vmix_to_fit, 
                            p0 = [vmixes[ii][peak_indexes[ii]], vmixes[ii][peak_indexes[ii]], #fit parameters S, A
                                    np.abs(fields[ii][peak_ranges[ii][-1]]-fields[ii][peak_ranges[ii][0]]), #W
                                    fields[ii][peak_indexes[ii]], 0], maxfev=5000) #H0, Offset

    extent = int(1.2*np.abs(prelim_popt[2])/np.abs(fields[ii][1]-fields[ii][0]))  #setting size of the final fit, roundingto nearest integer
    B0_index = np.argmin([np.abs(field - prelim_popt[3]) for field in fields[ii]])   #Finding true center of the fit
    fit_range_max = B0_index + extent
    fit_range_min = B0_index - extent
    if fit_range_max > len(fields[ii]): fit_range_max = len(fields[ii]) #making sure we only fit to data we have...
    if fit_range_min < 0: fit_range_min = 0
    
    fit_ranges.append(np.arange(fit_range_min, fit_range_max))

    #now to change everything to array again for easier plotting
    fit_ranges_to_plot = fit_ranges[ii]
    fit_ranges_arr = np.array(fit_ranges_to_plot)
    fields_to_fit_2 = np.array(field_arr[fit_ranges_arr])
    vmix_to_fit_2 = np.array(vmix_arr[fit_ranges_arr])
    popt, pcov = curve_fit(stfmr_fit, fields_to_fit_2, vmix_to_fit_2, p0 = prelim_popt)    
    
    S_list.append(popt[0]) #build lists here from fit
    A_list.append(popt[1])
    W_list.append(popt[2])
    H0_list.append(popt[3])
    Offset_list.append(popt[4])
"""

"""
#NEXT THING TO PLOT IS THE KITTEL_FIT
Kittel_popt, Kittel_pcov = curve_fit(Kittel_fit, np.abs(H0_list), frequencies, p0 = [Ms_lst[0]/2]) #Meff
Meff = Kittel_popt[0]
fig2, ax2 = plt.subplots(figsize = (3,3), dpi = 200)
ax2.plot(np.abs(H0_list),frequencies,'sk', label = 'data')
H0s = np.linspace(0,np.max(np.abs(H0_list)),1000)
ax2.plot(H0s,Kittel_fit(np.abs(H0s),Meff), '-b',label = 'Kittel fit')
ax2.set_xlabel(r'$H_0$ (Oe)')
ax2.set_ylabel("Frequency (GHz)")
ax2.legend()
print("M_eff = %.2e" % (Meff)) #Gives us M_eff

#NEXT THING TO PLOT IS THE SHA AS FUNCTION OF FREQUENCIES
#USE DATA FROM GOOGLE DOC MAG_THICK ETC
'''
NM_thick = 50 #nm
FM_thick = 4.13 '''
#NOTE CURRENTLY ONLY SUPPORTS 1 SAMPLE AT A TIME
SHA_list = []
fig3, ax3 = plt.subplots(figsize = (3,3), dpi = 200)
for ii in range(number_of_files):
    SHA_list.append(Calculate_SHA(S_list[ii], A_list[ii], Ms_lst[0], Con_thick_lst[0]*1e-7, Mag_thick_lst[0]*1e-7, Meff, H0_list[ii]))
ax3.plot(frequencies,SHA_list,'ks')
ax3.set_xlabel("Frequency (GHz)")
ax3.set_ylabel(r"$\theta_{SH}$")
#ax3.set_yscale("log")
print(SHA_list)
print("SHA: %.2f +/- %.2f" %(np.mean(SHA_list), np.std(SHA_list)))

#LAST IS OPTIONAL AND SHOWS HOW THE PEAKS AND SMOOTHING FUNCTIONS WORK
fig4, ax4 = plt.subplots()
for l in range(number_of_files):
    ax4.plot(fields[l], vmixes[l], '-', label = freq_dir[l][:-4])
    ax4.plot(fields[l], vmixes_smooth[l], label = freq_dir[l][:-4] + 'smooth')
    field = fields[l]
    vmix = vmixes_smooth[l]
    #PLOTS THE PEEAKS USED FOR THE STFMR FIT (RANGE AROUND PEAK)
    ax4.plot(field[peaks_lst[l][correct_peaks[l]]], vmix[peaks_lst[l][correct_peaks[l]]], 'g.')
    #PLOTS ALL PEAK POINTS AS BLACK
    for yy in range(len(peaks_lst[l])):
        ax4.plot(field[peaks_lst[l][yy]], vmix[peaks_lst[l][yy]], 'k.')
    #PLOTS ALL TROUGH POINTS AS BLUE
    for xx in range(len(troughs_lst[l])):
        ax4.plot(field[troughs_lst[l][xx]], vmix[troughs_lst[l][xx]], 'b.')
ax4.legend()
"""
#plt.show()




