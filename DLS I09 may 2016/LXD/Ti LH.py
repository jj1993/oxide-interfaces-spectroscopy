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

LH2in = "csv/LMO2TiLin-1.csv"
LH2ex = "csv/LMO2TiLex-1.csv"
LH4in = "csv/LMO4TiLin-1.csv"
LH7in = "csv/LMO7TiLin-1.csv"
LH25ex = "csv/LMO25TiLex-1.csv"

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
    
    engy_2in, spect_2in = load(LH2in)
    engy_2ex, spect_2ex = load(LH2ex)
    engy_4in, spect_4in = load(LH4in)
    engy_7in, spect_7in = load(LH7in)
    engy_25ex, spect_25ex = load(LH25ex)

    engy_2ex, spect_2expol = load(XLD2ex)
    engy_4in, spect_4inpol = load(XLD4in)
    engy_7in, spect_7inpol = load(XLD7in)
    engy_25ex, spect_25expol = load(XLD25ex)

    ref = np.mean([float(s) for n, s in enumerate(spect_2in) if float(engy_2in[n]) > 480])
    spect_2in = [float(s)*.2/ref for s in spect_2in]
    ref = np.mean([float(s) for n, s in enumerate(spect_2ex) if float(engy_2ex[n]) > 480])
    print .2/ref
    spect_2ex = [float(s)*.2/ref for s in spect_2ex]
    spect_2expol = [float(s)*.2/ref for s in spect_2expol]
    ref = np.mean([float(s) for n, s in enumerate(spect_4in) if float(engy_4in[n]) > 480])
    print .2/ref
    spect_4in = [float(s)*.2/ref for s in spect_4in]
    spect_4inpol = [float(s)*.2/ref for s in spect_4inpol]
    ref = np.mean([float(s) for n, s in enumerate(spect_7in) if float(engy_7in[n]) > 480])
    print .2/ref
    spect_7in = [float(s)*.2/ref for s in spect_7in]
    spect_7inpol = [float(s)*.2/ref for s in spect_7inpol]
    ref = np.mean([float(s) for n, s in enumerate(spect_25ex) if float(engy_25ex[n]) > 480])
    print .2/ref
    spect_25ex = [float(s)*.2/ref for s in spect_25ex]
    spect_25expol = [float(s)*.2/ref for s in spect_25expol]

    plt.plot(engy_2ex, [0 for e in engy_2ex])
    plt.plot(engy_2ex, spect_2ex, 'm--', linewidth=2.0, label="LH for ex-situ $STO|LMO$ (2uc) ")
##    plt.plot(engy_2in, spect_2in, 'k', linewidth=2.0, label="LV for in-situ $STO|LMO_2$")
    plt.plot(engy_4in, spect_4in, 'g-', linewidth=2.0, label="LH for in-situ $STO|LMO$ (4uc)")
    plt.plot(engy_7in, spect_7in, 'b-', linewidth=2.0, label="LH for in-situ $STO|LMO$ (7uc)")
    plt.plot(engy_25ex, spect_25ex, 'r--',linewidth=2.0,  label="LH for ex-situ $STO|LMO$ (25uc)")
    
    legend = plt.legend(loc='bottom right', shadow=True)
    title = plt.title("Ti L-edge LH data for different samples")
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("Intensity (a.u.)", fontsize=18)
    plt.show()

    plt.plot(engy_2ex, [0 for e in engy_2ex])
    plt.plot(engy_2ex, spect_2expol, 'm--', linewidth=2.0, label="XLD for ex-situ $STO|LMO$ (2uc)")
##    plt.plot(engy_2in, spect_2in, label="XLD for in-situ $STO|LMO_2$")
    plt.plot(engy_4in, spect_4inpol, 'g-', linewidth=2.0, label="XLD for in-situ $STO|LMO$ (4uc)")
    plt.plot(engy_7in, spect_7inpol, 'b-', linewidth=2.0, label="XLD for in-situ $STO|LMO$ (7uc)")
    plt.plot(engy_25ex, spect_25expol, 'r--', linewidth=2.0, label="XLD for ex-situ $STO|LMO$ (25uc)")
    
    legend = plt.legend(loc='bottom right', shadow=True)
    title = plt.title("Ti L-edge XLD data for different samples")
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("Intensity (a.u.)", fontsize=18)
    plt.show()
