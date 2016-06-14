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

smooth = .0003
FILE_LH = "../../Data/i09-75025.nxs"

import csv
def save(energies, spectrum, ext):
    name = 'LMO2OKex'
    with open('csv/'+name+'-'+ext+'.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(('energy','intensity'))
        for n, s in enumerate(spectrum):
            writer.writerow((energies[n], s))
        f.close()

def getBackground(energies, spectrum):
    
    last = len(energies)
    print 'LAST: ',last
    points = [2, 40, 65, 105, 222, 233]
    B = 1
    answ = "n"

    # Substracting pre-edge
    minimum = min(spectrum[points[0]:points[1]])
    spectrum = [s-minimum for s in spectrum]

    # Simulating step-function
    step = []
    start = spectrum[points[1]]
    mid = spectrum[points[3]]
    end = spectrum[points[4]]
    diff1 = mid - start
    diff2 = end - mid
    for e in energies:
        # fermi-dirac distribution-like normalisation function
        n = diff1/(math.exp((energies[points[2]]-e)/B)+1)
        step.append(n)
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
    bg = range(points[0],points[1]) + [points[3]] + range(points[4],points[5])
    x, y = [energies[n] for n in bg], [fitspect[n] for n in bg]
    a, b, c, d, e =  np.polyfit(x, y, 4)
    norm = [a*i**4 + b*i**3 + c*i**2 + d*i + e for i in energies]
    spectNorm = [s-norm[n]+step[n] for n, s in enumerate(fitspect)]
    
    return norm, spectNorm

if __name__ == '__main__':
    print "Gathering measurements data"
    
    with h5py.File(FILE_LH, 'r') as f:
        energies = f["/entry1/instrument/pathgroup/pgmenergy"]
        
        spectrum = f["/entry1/instrument/smpmiamp39/smpmiamp39"]
        i_0 = f["entry1/instrument/sm5iamp8/sm5iamp8"]

        LHEnergies = [e for e in energies[2:-2]]
        LHSpectrum = [s/i_0[n] for n, s in enumerate(spectrum[2:-2])]

    energies = np.linspace(LHEnergies[0], LHEnergies[-1], num=len(LHEnergies)*5, endpoint=True)

    print "============"
    print "First background"
    print "============"
    LHNorm, LHSpectrum_norm = getBackground(LHEnergies, LHSpectrum)
    splineRep = interpolate.splrep(LHEnergies, LHSpectrum_norm, s=smooth)
    LHSpectrum_norm = interpolate.splev(energies, splineRep, der=0)

    # Save data to CSV
    spectra = [LHSpectrum_norm]
    for n, s in enumerate(spectra):
        save(energies, s, str(n+1))

    # Plot data
    plt.plot(energies, [0 for e in energies])
##    plt.plot(LHEnergies, LHSpectrum)
##    plt.plot(LHEnergies, LHNorm)
    plt.plot(energies, LHSpectrum_norm, label='LH')
    
    legend = plt.legend(loc='upper right', shadow=True)
    title = plt.title("O K XLD spectrum of ex-situ STO|LMO(2)")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity (a.u.)")
    plt.show()   


