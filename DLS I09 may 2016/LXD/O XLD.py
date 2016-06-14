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

XLD7in = "csv/LMO7OKin-3.csv"

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
    
    engy_7in, spect_7in = load(XLD7in)

    plt.plot(engy_7in, [0 for e in engy_7in])
    plt.plot(engy_7in, spect_7in, label="XLD for in-situ $STO|LMO_{7}$")
    
    legend = plt.legend(loc='bottom right', shadow=True)
    title = plt.title("O K-edge XLD data")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity (a.u.)")
    plt.show()

