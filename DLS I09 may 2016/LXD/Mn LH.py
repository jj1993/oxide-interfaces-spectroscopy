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
from scipy import interpolate

smooth = 0.01

LH2in = "csv/LMO2MnLin-1.csv"
LV2in = "csv/LMO2MnLin-2.csv"
LH2ex = "csv/LMO2MnLex-1.csv"
LV2ex = "csv/LMO2MnLex-2.csv"
LH4in = "csv/LMO4MnLin-1.csv"
LV4in = "csv/LMO4MnLin-2.csv"
LH7in = "csv/LMO7MnLin-1.csv"
LV7in = "csv/LMO7MnLin-2.csv"
LH10ex = "csv/LMO10MnLex-1.csv"
LV10ex = "csv/LMO10MnLex-2.csv"

XLD2inpol = "csv/LMO2MnLin-3.csv"
XLD2expol = "csv/LMO2MnLex-3.csv"
XLD4inpol = "csv/LMO4MnLin-3.csv"
XLD7inpol = "csv/LMO7MnLin-3.csv"
XLD10expol = "csv/LMO10MnLex-3.csv"

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

def getScore(intensity, spect):
    s = 0
    for n, i in enumerate(intensity):
        diff = i - spect[n]
        s += diff**2
    return s

def getProportion(spectra, intensity, energy):

    l = [0.05*i for i in range(20)]
    bestScore = float("inf")
    print 'checking proportion'
    for e in l:
        for r in l:
            for t in l:
                if e + r + t == 1:
                    thisSpect = []
                    for n, s in enumerate(spectra[0]):
                        thisSpect.append(1*(s*e + spectra[1][n]*r + spectra[2][n]*t))
                    score = getScore(intensity, thisSpect)
                    if score < bestScore:
                        bestScore = score
                        proportion = (e, r, t)
                        spect = thisSpect

    print "Best proportion:",proportion
    
    return spect

if __name__ == '__main__':
    print "Gathering measurements data"
    
    # Manganese things
    f = open("../Mn-valences/energies.csv").read()
    L_energies = [float(e) for n, e in enumerate(f.split(","))]
    L_energies = [e+.3 for e in L_energies]
    
    f = open("../Mn-valences/Mn(II).csv").read()
    MNII = [float(s) for s in f.split(",")]
    for n, s in enumerate(MNII):
        if not s > 0:
            end = n-1
            break
    energies = [e for n, e in enumerate(L_energies) if n <= end]
    splineRep = interpolate.splrep(L_energies, MNII, s=smooth)
    MNII = interpolate.splev(energies, splineRep, der=0)
    for n, s in enumerate(MNII):
        if int(energies[n]) == 653:
            that = s
            ref = n
    
    f = open("../Mn-valences/Mn(III).csv").read()
    L_energies = [e-.7 for e in L_energies]
    MNIII = [float(i) for n, i in enumerate(f.split(","))]
    splineRep = interpolate.splrep(L_energies, MNIII, s=smooth)
    MNIII = interpolate.splev(energies, splineRep, der=0)
    this = MNIII[ref]
    MNIII = [s*that/this for s in MNIII]
    
    f = open("../Mn-valences/Mn(IV).csv").read()
    MNIV = [float(i) for n, i in enumerate(f.split(","))]
    L_energies = [e-.4 for e in L_energies]
    splineRep = interpolate.splrep(L_energies, MNIV, s=smooth)
    MNIV = interpolate.splev(energies, splineRep, der=0)
    this = MNIV[ref]
    MNIV = [s*that/this for s in MNIV]

    spectra = [MNII, MNIII, MNIV]

    # Data normalisation and to same energie scale
    engy_2in, spect_2in = load(LH2in)
    engy, LV = load(LV2in)
    engy_2in = [float(s) for s in engy_2in]
    spect_2in = [2*float(s)+float(LV[n]) for n, s in enumerate(spect_2in)]
    splineRep = interpolate.splrep(engy_2in, spect_2in, s=smooth)
    spect_2in = interpolate.splev(energies, splineRep, der=0)
    this = spect_2in[ref]
    spect_2in = [s*that/this for s in spect_2in]
    
    engy_2ex, spect_2ex = load(LH2ex)
    engy, LV = load(LV2ex)
    engy_2ex = [float(s) for s in engy_2ex]
    spect_2ex = [2*float(s)+float(LV[n]) for n, s in enumerate(spect_2ex)]
    splineRep = interpolate.splrep(engy_2ex, spect_2ex, s=smooth)
    spect_2ex = interpolate.splev(energies, splineRep, der=0)
    this = spect_2ex[ref]
    spect_2ex = [s*that/this for s in spect_2ex]
    
    engy_4in, spect_4in = load(LH4in)
    engy, LV = load(LV4in)
    engy_4in = [float(s) for s in engy_4in]
    spect_4in = [2*float(s)+float(LV[n]) for n, s in enumerate(spect_4in)]
    splineRep = interpolate.splrep(engy_4in, spect_4in, s=smooth)
    spect_4in = interpolate.splev(energies, splineRep, der=0)
    this = spect_4in[ref]
    spect_4in = [s*that/this for s in spect_4in]
    
    engy_7in, spect_7in = load(LH7in)
    engy, LV = load(LV7in)
    engy_7in = [float(s) for s in engy_7in]
    spect_7in = [2*float(s)+float(LV[n]) for n, s in enumerate(spect_7in)]
    splineRep = interpolate.splrep(engy_7in, spect_7in, s=smooth)
    spect_7in = interpolate.splev(energies, splineRep, der=0)
    this = spect_7in[ref]
    spect_7in = [s*that/this for s in spect_7in]
    
    engy_10ex, spect_10ex = load(LH10ex)
    engy, LV = load(LV10ex)
    engy_10ex = [float(s) for s in engy_10ex]
    spect_10ex = [2*float(s)+float(LV[n]) for n, s in enumerate(spect_10ex)]
    splineRep = interpolate.splrep(engy_10ex, spect_10ex, s=smooth)
    spect_10ex = interpolate.splev(energies, splineRep, der=0)
    this = spect_10ex[ref]
    spect_10ex = [s*that/this for s in spect_10ex]
    
    best = getProportion(spectra, spect_7in, energies)
    MNII = [s+30 for s in MNII]
    MNIII = [s+10 for s in MNIII]

##    plt.plot(energies, MNII, linewidth=2.0, label="$MnO$ (divalent)")
##    plt.plot(energies, MNIII, linewidth=2.0, label="$LaMnO_3$ (trivalent)")
##    plt.plot(energies, MNIV, linewidth=2.0, label="$SrMnO_3$ (tetravalent)")
    plt.plot(energies, spect_2in, linewidth=2.0, label="In-situ $STO|LMO$ (2uc)")
##    plt.plot(energies, best, linewidth=2.0, label="Best fit of di-, tri- and tetravalent Mn")
    plt.plot(energies, [0 for e in energies], 'k-')
    plt.plot(energies, spect_2ex, 'm--', linewidth=2.0, label="Ex-situ $STO|LMO$ (2uc)")
    plt.plot(energies, spect_4in, 'g-', linewidth=2.0, label="In-situ $STO|LMO$ (4uc)")
    plt.plot(energies, spect_7in, linewidth=2.0, label="In-situ $STO|LMO$ (7uc)")
    plt.plot(energies, spect_10ex, 'r--', linewidth=2.0, label="Ex-situ $STO|LMO$ (10uc)")


    legend = plt.legend(loc='bottom right', fontsize=16, shadow=True)
    plt.xlabel("Photon energy (eV)", fontsize=18)
    plt.ylabel("Drain current intensity (arb. units)", fontsize=16)
    plt.show()

##    
##    legend = plt.legend(loc='bottom right', shadow=True)
##    title = plt.title("Mn L-edge LH data for different samples")
##    plt.xlabel("Energy (eV)", fontsize=18)
##    plt.ylabel("Intensity (a.u.)", fontsize=18)
##    plt.show()

    

##    plt.plot(engy_2in, [0 for e in engy_2in])
##    plt.plot(engy_2in, spect_2inpol, 'k-', linewidth=2.0, label="XLD for in-situ $STO|LMO$ (2uc)")
##    plt.plot(engy_2ex, spect_2expol, 'm--', linewidth=2.0, label="XLD for ex-situ $STO|LMO$ (2uc)")
##    plt.plot(engy_4in, spect_4inpol, 'g-', linewidth=2.0, label="XLD for in-situ $STO|LMO$ (4uc)")
##    plt.plot(engy_7in, spect_7inpol, 'b-', linewidth=2.0, label="XLD for in-situ $STO|LMO$ (7uc)")
##    plt.plot(engy_10ex, spect_10expol, 'r--', linewidth=2.0, label="XLD for ex-situ $STO|LMO$ (10uc)")
##    
##    legend = plt.legend(loc='bottom right', shadow=True)
##    title = plt.title("Mn L-edge XLD data for different samples")
##    plt.xlabel("Energy (eV)", fontsize=18)
##    plt.ylabel("Intensity (a.u.)", fontsize=18)
##    plt.show()
##
