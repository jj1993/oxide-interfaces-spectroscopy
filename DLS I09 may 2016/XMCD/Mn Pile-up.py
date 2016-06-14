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

LC4in = "csv/LMO4MNin-1.csv"
LC4ex = "csv/LMO4MNex-1.csv"
LC7 = "csv/LMO7MN-1.csv"
LC7_RT = "csv/LMO7MN-1.csv"
LC10 = "csv/LMO10MN-1.csv"

XMCD4in = "csv/LMO4MNin-3.csv"
XMCD4ex = "csv/LMO4MNex-5.csv"
XMCD7 = "csv/LMO7MN-7.csv"
XMCD7_RT = "csv/LMO7MN-9.csv"
XMCD10 = "csv/LMO10MN-5.csv"

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
    
    engy_4in, spect_4in = load(XMCD4in)
    engy_4ex, spect_4ex = load(XMCD4ex)
    engy_7, spect_7 = load(XMCD7)
    engy_7rt, spect_7rt = load(XMCD7_RT)
    engy_10, spect_10 = load(XMCD10)

    engy_4in, spect_4inlc = load(LC4in)
    engy_4ex, spect_4exlc = load(LC4ex)
    engy_7, spect_7lc = load(LC7)
    engy_7rt, spect_7rtlc = load(LC7_RT)
    engy_10, spect_10lc = load(LC10)

    plt.plot(engy_4in, [0 for e in engy_4in])
##    plt.plot(engy_4in, spect_4in, label="$\sigma^+-\sigma^-$ for $STO|LMO_4$ in-situ at LT")
    plt.plot(engy_4ex, spect_4ex, 'g-', linewidth=2.0, label="$\sigma^+-\sigma^-$ for $STO|LMO$ (4uc) ex-situ at LT")
    plt.plot(engy_7, spect_7, 'b-', linewidth=2.0, label="$\sigma^+-\sigma^-$ for $STO|LMO$ (7uc) at LT")
    plt.plot(engy_7rt, spect_7rt, 'c-', linewidth=2.0, label="$\sigma^+-\sigma^-$ for $STO|LMO$ (7uc) at RT")
    plt.plot(engy_10, spect_10, 'r-', linewidth=2.0, label="$\sigma^+-\sigma^-$ for $STO|LMO$ (10uc) at LT")
    
    legend = plt.legend(loc='bottom right', shadow=True)
    title = plt.title("Mn L-edge XMCD data for different samples and temperatures")
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("Ratio $(\sigma^+-\sigma^-)/max(\sigma^-)$ ", fontsize=18)
    plt.show()

