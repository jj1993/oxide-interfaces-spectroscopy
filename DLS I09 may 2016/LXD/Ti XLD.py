"""
Script to analyse XLD and XLD spectra
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

##XLD2in = "csv/LMO2TiLin-3.csv" #LH Dataset aborted during measurement...
XLD2ex = "csv/LMO2TiLex-3.csv"
XLD4in = "csv/LMO4TiLin-3.csv"
XLD7in = "csv/LMO7TiLin-3.csv"
XLD25ex = "csv/LMO25TiLex-3.csv"

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
    
##    engy_2in, spect_2in = load(XLD2in)
    engy_2ex, spect_2ex = load(XLD2ex)
    engy_4in, spect_4in = load(XLD4in)
    engy_7in, spect_7in = load(XLD7in)
    engy_25ex, spect_25ex = load(XLD25ex)

    plt.plot(engy_2ex, [0 for e in engy_2ex])
    plt.plot(engy_2ex, spect_2ex, '--', label="XLD for ex-situ $STO|LMO_2$ ")
##    plt.plot(engy_2in, spect_2in, label="XLD for in-situ $STO|LMO_2$")
    plt.plot(engy_4in, spect_4in, label="XLD for in-situ $STO|LMO_4$")
    plt.plot(engy_7in, spect_7in, label="XLD for in-situ $STO|LMO_{7}$")
    plt.plot(engy_25ex, spect_25ex, '--', label="XLD for ex-situ $STO|LMO_{25}$")
    
    legend = plt.legend(loc='bottom right', shadow=True)
    title = plt.title("Ti L-edge XLD data for different samples")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity (a.u.)")
    plt.show()

