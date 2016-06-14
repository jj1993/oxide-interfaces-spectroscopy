"""
Script to analyse LH and LH spectra
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

LH2in = "csv/LMO2OKin-1.csv"
LH2ex = "csv/LMO2OKex-1.csv"
LH4in = "csv/LMO4OKin-1.csv"
LH7in = "csv/LMO7OKin-1.csv"

XLD7inpol = "csv/LMO7OKin-3.csv"

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

    engy_2in, spect_2in = load(LH2in)
    engy_2ex, spect_2ex = load(LH2ex)
    engy_4in, spect_4in = load(LH4in)
    engy_7in, spect_7in = load(LH7in)

    engy_7in, spect_7inpol = load(XLD7inpol)

    ref = np.mean([float(s) for n, s in enumerate(spect_2in) if float(engy_2in[n]) > 555])
    print .2/ref
    spect_2in = [float(s)*.2/ref for s in spect_2in]
    ref = np.mean([float(s) for n, s in enumerate(spect_2ex) if float(engy_2ex[n]) > 555])
    print .2/ref
    spect_2ex = [float(s)*.2/ref for s in spect_2ex]
    ref = np.mean([float(s) for n, s in enumerate(spect_4in) if float(engy_4in[n]) > 555])
    print .2/ref
    spect_4in = [float(s)*.2/ref for s in spect_4in]
    ref = np.mean([float(s) for n, s in enumerate(spect_7in) if float(engy_7in[n]) > 555])
    print .2/ref
    spect_7in = [float(s)*.2/ref for s in spect_7in]
    spect_7inpol = [float(s)*.2/ref for s in spect_7inpol]

    plt.plot(engy_7in, [0 for e in engy_7in])
    plt.plot(engy_2in, spect_2in, 'k-', linewidth=2.0, label="LH for in-situ $STO|LMO$ (2uc)")
    plt.plot(engy_2ex, spect_2ex, 'm--',linewidth=2.0, label="LH for ex-situ $STO|LMO$ (2uc)")
    plt.plot(engy_4in, spect_4in, 'g-', linewidth=2.0, label="LH for in-situ $STO|LMO$ (4uc)")
    plt.plot(engy_7in, spect_7in, 'b-', linewidth=2.0, label="LH for in-situ $STO|LMO$ (7uc)")
    
    legend = plt.legend(loc='bottom right', shadow=True)
    title = plt.title("O K-edge LH data")
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("Intensity (a.u.)", fontsize=18)
    plt.show()

    plt.plot(engy_7in, [0 for e in engy_7in])
    plt.plot(engy_7in, spect_7inpol, 'b-',linewidth=2.0, label="XLD for in-situ $STO|LMO$ (7uc)")
    
    legend = plt.legend(loc='bottom right', shadow=True)
    title = plt.title("O K-edge XLD data")
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("Intensity (a.u.)", fontsize=18)
    plt.show()

