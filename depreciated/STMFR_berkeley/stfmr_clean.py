import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import scipy.signal as sci
from scipy.optimize import curve_fit
from gdstfmr import*



if len(sys.argv) == 1: directories = ["BPBO-LSMO_072121b\\dev13"]
else: directories = sys.argv[1:] 

num_dir = len(directories)
freq_dir = []
fields = []
vmixes = []
for i in range(num_dir):
    directory = directories[i]
    files = os.listdir(directory)
    number_of_files = len(files)

    for a in range(number_of_files):
            freq_dir.append(files[a])

    for j in range(number_of_files):
        with open(directory + '\\' + freq_dir[j], 'r') as file:
            lines = file.readlines()
        
        print(lines[6]) #starting point for data

        len_of_data = len(lines) - 6

        print(len_of_data)

        field=[]
        vmix = []

        for q in range(len_of_data):
            temp = np.fromstring(lines[q+6], dtype=float, sep='\t')
            field.append(temp[0])
            vmix.append(temp[1])
            #print(temp)
        fields.append(field)
        vmixes.append(vmix)

print(len(fields), len(vmixes))
#initialize all variables
Ms = 800 #emu/cc #use ISAACS VALUE BECAUSE OUR SQUID DATA IS HIGH AF

H0_list = [] #center list
S_list = [] #symmetric list
A_list = [] #Anti-symmetric list
W_list = [] #width list
Offset_list = []

frequencies = []
labels = []
for mm in range(number_of_files):
    labels.append(freq_dir[mm][:-4])
print(labels)
frequencies_temp = []
tsty = []
gg = 0
yy = 0
while gg < number_of_files:
    uncleanfreq = freq_dir[gg]
    for jj in uncleanfreq:
        if jj == 'G':
            tsty = uncleanfreq[:yy]
            break
        else:
            yy += 1
            #frequencies_temp.append(jj)
    frequencies_temp.append(tsty)
    yy = 0
    gg+= 1
print(frequencies_temp)
for uu in range(number_of_files):
    frequencies.append(float(frequencies_temp[uu]))

print(frequencies, 'here')

#New Functions for STMFR analysis:
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

#Write code here to find peak positions bet
#first lets smooth the data with savgol filter
vmixes_smooth = []
for ll in range(number_of_files):
    yuh = sci.savgol_filter(vmixes[ll], 51, 10)
    vmixes_smooth.append(yuh)

Peak_Vmix_index = []
for ii in range(number_of_files):
    min1 = np.argmin(vmixes_smooth[ii])
    Peak_Vmix_index.append(min1)
#print(Peak_Vmix_index)
min_vmixes = []
for oo in range(number_of_files):
    bet = vmixes[oo][Peak_Vmix_index[oo]]
    min_vmixes.append(bet)
#print(min_vmixes)

#use find_peaks from scipy to get peak_vmix_index
peaks_lst = []
neg_peaks_lst = []
for rr in range(number_of_files):
    peaks, _ = sci.find_peaks(vmixes_smooth[rr], width=8)
    neg_peaks, _ = sci.find_peaks(-1*vmixes_smooth[rr], width=8)
    peaks_lst.append(peaks)
    neg_peaks_lst.append(neg_peaks)
#print(peaks_lst)
#print(neg_peaks_lst)

peaks_lst_true = []
neg_peaks_lst_true = []
for bb in range(number_of_files):
    temp = peaks_lst[bb]
    temp1 = neg_peaks_lst[bb]
    temp2 = list(temp)
    temp3 = list(temp1)
    peaks_lst_true.append(temp2)
    neg_peaks_lst_true.append(temp3)

#code to ensure that neg-peaks and peaks are same length
for qq in range(number_of_files):
    if len(neg_peaks_lst_true[qq]) > len(peaks_lst_true[qq]):
        neg_peaks_lst_true[qq] = neg_peaks_lst_true[qq][1:]
    elif len(neg_peaks_lst_true[qq]) < len(peaks_lst_true[qq]):
        peaks_lst_true[qq].pop()

#by this point we have all the peaks and troughs so we need to select the best ones
#we can try to use armin or argmax to get index of correct peaks to use
correct_peaks_back = []
correct_troughs_back = []
for ww in range(number_of_files):
    vmix = vmixes_smooth[ww]
    correct_peaks_back.append(np.argmax(vmix[peaks_lst_true[ww]]))
    correct_troughs_back.append(np.argmin(vmix[neg_peaks_lst_true[ww]]))

correct_troughs = correct_troughs_back
correct_peaks = correct_peaks_back

#please note lists were index from right to left for some reason on graph so i reversed them

#NOW WE NEED TO MAKE IT TO TAKE A RANGE AROUND THE VMIX INDEX
peak_ranges = []
trough_ranges=[]
peak_indexes = []
trough_indexes = []
for mm in range(number_of_files):
    peakvmixindex = peaks_lst_true[mm][correct_peaks[mm]]
    troughvmixindex = neg_peaks_lst_true[mm][correct_troughs[mm]]
    range_around_peak_vmixindex = np.arange(peakvmixindex-15, peakvmixindex+15)
    range_around_trough_vmixindex = np.arange(troughvmixindex-15, troughvmixindex+15)
    trough_ranges.append(range_around_trough_vmixindex)
    peak_ranges.append(range_around_peak_vmixindex)
    peak_indexes.append(peakvmixindex)
    trough_indexes.append(troughvmixindex) 

#BASICALLY JUST MADE A LIST OF THE RANGES AROUND THE PEAKS AND TROUGHS AND MADE A LIST OF THE CORRECT ONES FOR ALL FILES

peak_ranges_lst = []
troughs_ranges_lst = []
for mm in range(number_of_files):
    betsy = peak_ranges[mm]
    chloe = trough_ranges[mm]
    peak_ranges_lst.append(list(betsy))
    troughs_ranges_lst.append(list(chloe))
#AT THIS POINT WHAT THE WAT THE CODE SHOULD BE DOING IS GIVING INDEX RANGE AROUND SIGNAL FOR MAX OR TROUGH

#essentially time to change everything back to an array because im dumb and lists no bueno
#stuff to fit below goes here (im only going to do for peaks, ill recopy it for troughs later)

#IM GOING TO MAKE A BOOL TO ALLOW FOR EITHER TROUGHS OR PEAKS TO BE USED TO FIT DATA
troughs_bool = False


#copypasta of isaacs code and trying to make it work with my variables
fit_ranges = []
for ii in range(number_of_files):
    print("Fitting File %s" %labels[ii])
    #put stuff to plot here to simplify the curve fit
    if troughs_bool:
        woop = np.array(troughs_ranges_lst[ii])
    else:
        woop = np.array(peak_ranges_lst[ii])
    field = fields[ii]
    field_arr = np.array(field)
    fields_to_fit = np.array(field_arr[woop])
    vmix = vmixes[ii]#changed to vmix
    vmix_arr = np.array(vmix)
    vmix_to_fit = np.array(vmix_arr[woop])

    #DEPENDING ON THE BOOL, FITS STMFR TO TROUGH VERSUS PEAK FIT

    if troughs_bool:
        prelim_popt, prelim_pcov = curve_fit(stfmr_fit, fields_to_fit, vmix_to_fit, 
                           p0 = [vmixes[ii][trough_indexes[ii]], vmixes[ii][trough_indexes[ii]], #fit parameters S, A
                                 np.abs(fields[ii][troughs_ranges_lst[ii][-1]]-fields[ii][troughs_ranges_lst[ii][0]]), #W
                                 fields[ii][trough_indexes[ii]], 0], maxfev=5000) #H0, Offset

    else:
        prelim_popt, prelim_pcov = curve_fit(stfmr_fit, fields_to_fit, vmix_to_fit, 
                            p0 = [vmixes[ii][peak_indexes[ii]], vmixes[ii][peak_indexes[ii]], #fit parameters S, A
                                    np.abs(fields[ii][peak_ranges_lst[ii][-1]]-fields[ii][peak_ranges_lst[ii][0]]), #W
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
    
    S_list.append(popt[0])
    A_list.append(popt[1])
    W_list.append(popt[2])
    H0_list.append(popt[3])
    Offset_list.append(popt[4])

#taking break here, need to plot actual fits still from stmfrmany
fig1, ax1 = plt.subplots(figsize = (4,5), dpi = 200) 


for ii in range(number_of_files):   
    H_lin = np.linspace(fields[ii][fit_ranges[ii][0]], fields[ii][fit_ranges[ii][-1]], 1000) #fields to show stfmr_fit function
    
    ax1.plot(fields[ii],vmixes[ii], '.', label = labels[ii])
    
    ax1.plot(H_lin,stfmr_fit(H_lin, S_list[ii], A_list[ii], W_list[ii], H0_list[ii], Offset_list[ii]),
             'k-') #plot it to make sure the fits worked reasonably well
    
    #Plotting the components is optional, and can get cluttered
    #ax1.plot(H_lin,symmetric(H_lin, S_list[ii], W_list[ii], H0_list[ii], Offset_list[ii]),color = 'gray', linestyle = '-')
    #ax1.plot(H_lin,antisymmetric(H_lin, A_list[ii], W_list[ii], H0_list[ii], Offset_list[ii]),color = 'gray', linestyle = '-')
    
    #.append(rangeAround_Peak_index)
    #peak_indexes.append(Peak_Vmix_index)
    
ax1.ticklabel_format(axis='y', style = 'scientific', scilimits = (2,-2))
ax1.set_xlabel("Field (Oe)")
ax1.set_ylabel("V$_{mix}$ (V)")
ax1.legend() 



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

#Home stretch: calculate and plot the spin hall angle as a function of frequency!!
#need to make to take from google docs
NM_thick = 50 #nm
FM_thick = 4.13
SHA_list = []
fig3, ax3 = plt.subplots(figsize = (3,3), dpi = 200)
for ii in range(number_of_files):
    SHA_list.append(Calculate_SHA(S_list[ii], A_list[ii], Ms, NM_thick*1e-7, FM_thick*1e-7, Meff, H0_list[ii]))
ax3.plot(frequencies,SHA_list,'ks')
ax3.set_xlabel("Frequency (GHz)")
ax3.set_ylabel(r"$\theta_{SH}$")
#ax3.set_yscale("log")

print("SHA: %.2f +/- %.2f" %(np.mean(SHA_list), np.std(SHA_list)))


print(correct_peaks)
print(peaks_lst_true)
fig, ax = plt.subplots()
for l in range(number_of_files):
    ax.plot(fields[l], vmixes[l], '-', label = freq_dir[l][:-4])
    ax.plot(fields[l], vmixes_smooth[l], label = freq_dir[l][:-4] + 'smooth')
    field = fields[l]
    vmix = vmixes_smooth[l]
    #ax.plot(peaks_lst_true[l], vmix[peaks_lst_true[l]])
    for yy in range(len(peaks_lst_true[l])):
        ax.plot(field[peaks_lst_true[l][yy]], vmix[peaks_lst_true[l][yy]], 'k.')
    #for xx in range(len(neg_peaks_lst_true)):
        ax.plot(field[neg_peaks_lst_true[l][yy]], vmix[neg_peaks_lst_true[l][yy]], 'b.')
        ax.plot(field[peaks_lst_true[l][correct_peaks[l]]], vmix[peaks_lst_true[l][correct_peaks[l]]], 'g.')
ax.legend()
#correct peaks is a simple list that has the index with respect to peaks_lst_true in order of files from 0-3

print(data)
plt.show()
#print(data)



