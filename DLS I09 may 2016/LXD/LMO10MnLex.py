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

smooth = .001
FILE_LH = "../../Data/i09-75371.nxs"
FILE_LV = "../../Data/i09-75372.nxs"

import csv
def save(energies, spectrum, ext):
    name = 'LMO10MnLex'
    with open('csv/'+name+'-'+ext+'.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(('energy','intensity'))
        for n, s in enumerate(spectrum):
            writer.writerow((energies[n], s))
        f.close()

def getBackground(energies, spectrum):

    last = len(energies)
    print last
    points = [1, 25, 65, 155, 185, 245, 278]
    B = 1
    answ = "n"

    # Substracting pre-edge
    minimum = min(spectrum[points[0]:points[1]])
    spectrum = [s-minimum for s in spectrum]

    # Simulating step-function
    step = []
    start = spectrum[points[1]]
    mid = spectrum[points[3]]
    end = spectrum[points[5]]
    diff1 = mid - start
    diff2 = end - mid
    for e in energies:
        # fermi-dirac distribution-like normalisation function
        n = diff1/(math.exp((energies[points[2]]-e)/B)+1)
        m = diff2/(math.exp((energies[points[4]]-e)/B)+1)
        step.append(n+m)
    fitspect = [s - step[n] for n, s in enumerate(spectrum)]

##    while answ == "n":
##        plt.plot(energies, spectrum)
##        plt.plot(energies, [s+start for s in step])
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

    # Fitting 4th order polynomial to bg without step
    bg = range(points[0],points[1]) + [points[3]] + range(points[5],points[6])
    x, y = [energies[n] for n in bg], [fitspect[n] for n in bg]
    a, b, c, d, e =  np.polyfit(x, y, 4)
    norm = [a*i**4 + b*i**3 + c*i**2 + d*i + e for i in energies]
    spectNorm = [s-norm[n]+step[n] for n, s in enumerate(fitspect)]
    
    return norm, spectNorm

if __name__ == '__main__':
    print "Gathering measurements data"
    
    with h5py.File(FILE_LH, 'r') as f:
        energies = f["/entry1/sm5amp8/pgmenergy"]
        
        spectrum = f["/entry1/smpmamp39/smpmamp39"]
        i_0 = f["/entry1/sm5amp8/sm5amp8"]

        LHEnergies = [e for e in energies[2:-2]]
        LHSpectrum = [s/i_0[n] for n, s in enumerate(spectrum[2:-2])]

    with h5py.File(FILE_LV, 'r') as f:
        energies = f["/entry1/sm5amp8/pgmenergy"]
        
        spectrum = f["/entry1/smpmamp39/smpmamp39"]
        i_0 = f["/entry1/sm5amp8/sm5amp8"]

        LVEnergies = [e for e in energies[2:-2]]
        LVSpectrum = [s/i_0[n] for n, s in enumerate(spectrum[2:-2])]

    energies = np.linspace(LHEnergies[0], LHEnergies[-1], num=len(LHEnergies)*5, endpoint=True)

    print "============"
    print "First background"
    print "============"
##    LHSpectrum_norm = LHSpectrum
    LHNorm, LHSpectrum_norm = getBackground(LHEnergies, LHSpectrum)
    splineRep = interpolate.splrep(LHEnergies, LHSpectrum_norm, s=smooth)
    LHSpectrum_norm = interpolate.splev(energies, splineRep, der=0)
    print "============"
    print "Second background"
    print "============"
##    LVSpectrum_norm = LVSpectrum
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
    plt.plot(LHEnergies, [0 for e in LHEnergies])
##    plt.plot(LVEnergies, LVSpectrum)
##    plt.plot(LVEnergies, LVNorm)
##    plt.plot(LHEnergies, LHSpectrum)
##    plt.plot(LHEnergies, LHNorm)
    plt.plot(energies, LHSpectrum_norm, label='LH')
    plt.plot(energies, LVSpectrum_norm, label='LV')
    plt.plot(energies, XLD, label='XLD, LH - LV')
    
    legend = plt.legend(loc='upper right', shadow=True)
    title = plt.title("Mn L XLD spectrum of ex-situ STO|LMO(10)")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity (a.u.)")
    plt.show()   


