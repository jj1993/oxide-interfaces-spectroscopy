"""
Script to analyse XLD and XMCD spectra
"""
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import math
import copy
import random
import h5py
import copy
import csv

XMCD4 = "csv/LMO4-OK-5.csv"
XMCD7 = "csv/LMO7-OK-7.csv"
XMCD7_RT = "csv/LMO7-OK-9.csv"
MAG4 = "csv/LMO4-OK-6.csv"
MAG7 =  "csv/LMO7-OK-10.csv"

LC4 = "csv/LMO4-OK-1.csv"
LC7 = "csv/LMO7-OK-1.csv"
LC7_RT = "csv/LMO7-OK-1.csv"

def load(File):
    energy = []
    spectrum = []
    with open(File, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            e, s = row
            energy.append(e)
            spectrum.append(s)

    return energy[1:], spectrum[1:]
    
def scaleIntensity(int_l, energy1, energy2):
    intensity = []
    for e1 in energy1:
        for n, e2 in enumerate(energy2):
            if e2 > e1:
                break
            e = e2
            
        # Build lineair fit between two nearest energies
        # I(e) = a*e + b
        try:
            a = (int_l[n] - int_l[n-1])/(e2 - e)
            b = int_l[n]-a*e2
            intensity.append(a*e1 + b)
        except: intensity.append(int_l[0])

    return intensity



if __name__ == '__main__':
    print "Gathering measurements data"
    
    engy_4, spect_4 = load(XMCD4)
    engy_7, spect_7 = load(XMCD7)
    engy_7rt, spect_7rt = load(XMCD7_RT)

    engy_4, spect_4lc = load(LC4)
    engy_7, spect_7lc = load(LC7)
    engy_7rt, spect_7rtlc = load(LC7_RT)

##    ref = np.mean([float(s) for n, s in enumerate(spect_4lc) if float(engy_4[n]) > 555])
##    spect_4 = [float(s)*.2/ref for s in spect_4]
##    ref = np.mean([float(s) for n, s in enumerate(spect_7lc) if float(engy_7[n]) > 555])
##    spect_7 = [float(s)*.2/ref for s in spect_7]
##    ref = np.mean([float(s) for n, s in enumerate(spect_7rtlc) if float(engy_7rt[n]) > 555])
##    spect_7rt = [float(s)*.2/ref for s in spect_7rt]

    plt.plot(engy_4, [0 for e in engy_4])
    plt.plot(engy_4, spect_4, 'g-',linewidth=2.0,  label="$\sigma^+-\sigma^-$ for $STO|LMO$ (4uc) at LT")
    plt.plot(engy_7, spect_7, 'b-', linewidth=2.0, label="$\sigma^+-\sigma^-$ for $STO|LMO$ (7uc) at LT")
    plt.plot(engy_7rt, spect_7rt, 'c-', linewidth=2.0, label="$\sigma^+-\sigma^-$ for $STO|LMO$ (7uc) at RT")
    
    legend = plt.legend(loc='upper right', shadow=True)
    title = plt.title("O K-edge XMCD data for different samples and temperatures")
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("Ratio $(\sigma^+-\sigma^-)/max(\sigma^-)$ ", fontsize=18)
    plt.show()

