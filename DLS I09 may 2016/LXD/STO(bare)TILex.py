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
from scipy import interpolate

smooth = .03
FILE_LH = "../../Data/i09-75012.nxs"
FILE_LV = "../../Data/i09-75013.nxs"

import csv
def save(energies, spectrum, ext):
    name = 'STO(bare)TiLex'
    with open('csv/'+name+'-'+ext+'.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(('energy','intensity'))
        for n, s in enumerate(spectrum):
            writer.writerow((energies[n], s))
        f.close()

def getBackground(energies, spectrum):
    
    points = [1, 140, 222, 295, 330, 450, 600]
    B = 1
    answ = "n"

##    while answ == "n":
##        plt.plot(energies, spectrum)
##        for p in points:
##            plt.axvline(energies[p], color='k', linestyle='solid')
##        plt.show()
##            
##        answ = raw_input("Are the lines placed correctly? y/n: ")
##        while answ != "n" and answ != "y":
##            answ = raw_input("Use 'y' or 'n'. Are the lines correct?: ")
##
##        if answ == "y": break
##        print "Previous lines were places at ",points
##        for n, p in enumerate(points):
##            try: points[n] = int(raw_input("Line "+str(n+1)+" to: "))
##            except:
##                print "Invalid syntax"
##                points[n] = int(raw_input("Line "+str(n+1)+" to: "))
##        print points

    bg = range(points[0],points[1]) + [points[3]] + range(points[5],points[6])
    step = []
    for e in energies:
        # fermi-dirac distribution-like normalisation function
        n = .5/(math.exp((energies[points[2]]-e)/B)+1)
        m = .5/(math.exp((energies[points[4]]-e)/B)+1)
        step.append(n+m)

    x, y = [energies[n] for n in bg], [spectrum[n] for n in bg]
    a, b, c, d, e =  np.polyfit(x, y, 4)
    norm = [a*i**4 + b*i**3 + c*i**2 + d*i + e for i in energies]
    spectNorm = [s-norm[n] for n, s in enumerate(spectrum)]

    try:
        spectNorm = [s/maximum+step[n] for n, s in enumerate(spectNorm)]
        print "USED OLD MAXIMUM"
    except:
        global maximum
        maximum= 0.2*max(spectNorm)
        spectNorm = [s/maximum+step[n] for n, s in enumerate(spectNorm)]
        
    return norm, spectNorm

if __name__ == '__main__':
    print "Gathering measurements data"
    
    with h5py.File(FILE_LH, 'r') as f:
        energies = f["/entry1/instrument/pathgroup/pgmenergy"]
        
        spectrum = f["/entry1/instrument/smpmiamp39/smpmiamp39"]
        i_0 = f["/entry1/instrument/sm5iamp8/sm5iamp8"]

        LHEnergies = [e for e in energies]
        LHSpectrum = [s/i_0[n] for n, s in enumerate(spectrum)]

    with h5py.File(FILE_LV, 'r') as f:
        energies = f["/entry1/instrument/pathgroup/pgmenergy"]
        
        spectrum = f["/entry1/instrument/smpmiamp39/smpmiamp39"]
        i_0 = f["/entry1/instrument/sm5iamp8/sm5iamp8"]

        LVEnergies = [e for e in energies]
        LVSpectrum = [s/i_0[n] for n, s in enumerate(spectrum)]

    energies = np.linspace(LHEnergies[0], LHEnergies[-1], num=len(LHEnergies)*5, endpoint=True)

    print "============"
    print "First background"
    print "============"
    LHNorm, LHSpectrum_norm = getBackground(LHEnergies, LHSpectrum)
    splineRep = interpolate.splrep(LHEnergies[2:-2], LHSpectrum_norm[2:-2], s=smooth)
    LHSpectrum_norm = interpolate.splev(energies, splineRep, der=0)
    print "============"
    print "Second background"
    print "============"
    LVNorm, LVSpectrum_norm = getBackground(LVEnergies, LVSpectrum)
    splineRep = interpolate.splrep(LVEnergies, LVSpectrum_norm, s=smooth)
    LVSpectrum_norm = interpolate.splev(energies, splineRep, der=0)

    # Make XMCD
    XLD = [(LHSpectrum_norm[n]-s)*4 for n, s in enumerate(LVSpectrum_norm)]

    # Save data to CSV
    spectra = [LHSpectrum_norm, LVSpectrum_norm, XLD]
    for n, s in enumerate(spectra):
        save(energies, s, str(n+1))

    # Plot data
    plt.plot(energies, [0 for e in energies])
##    plt.plot(LVEnergies, LHSpectrum)
##    plt.plot(LVEnergies, LHNorm)
    plt.plot(energies, LHSpectrum_norm, label='LH')
    plt.plot(energies, LVSpectrum_norm, label='LV')
    plt.plot(energies, XLD, label='XLD, LH - LV')
    
    legend = plt.legend(loc='upper right', shadow=True)
    title = plt.title("Ti L XLD spectrum of ex-situ bare STO")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity (a.u.)")
    plt.show()   


