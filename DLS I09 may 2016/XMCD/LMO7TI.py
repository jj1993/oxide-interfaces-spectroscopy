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
FILE_LC = "../../Data/i09-75693.nxs"
FILE_LCflip = "../../Data/i09-75713.nxs"
FILE_RC = "../../Data/i09-75694.nxs"
FILE_RCflip = "../../Data/i09-75714.nxs"
FILE_LC_RT = "../../Data/i09-75680.nxs"
FILE_RC_RT = "../../Data/i09-75681.nxs"
delta_E = 0.05 #eV, absolute
delta_int = 0.01 #relative

import csv
def save(energies, spectrum, ext):
    name = 'LMO7TI'
    with open('csv/'+name+'-'+ext+'.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(('energy','intensity'))
        for n, s in enumerate(spectrum):
            writer.writerow((energies[n], s))
        f.close()

def getBackground(energies, spectrum):
    
    points = [1, 100, 170, 245, 280, 400, 445]
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
    try:
        for e in energies:
            # fermi-dirac distribution-like normalisation function
            n = diff1/(math.exp((energies[points[2]]-e)/B)+1)
            m = diff2/(math.exp((energies[points[4]]-e)/B)+1)
            step.append(n+m)
    except:
        diff1 = mid - start
        diff2 = end - mid
        global diff1
        global diff2
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

    reg = [points[1], points[3], points[5]]
    return norm, spectNorm, reg

def getError(intensity, energies, delta_E, delta_int):
    
    error = []
    for n, e in enumerate(energies):
        try:
            a = (intensity[n+1] - intensity[n-1])/(energies[n+1]-energies[n-1])
            i_err = intensity[n] * delta_int
            e_err = a * delta_E
            error.append(math.sqrt(i_err**2 + e_err**2))
        except: error.append(0)

    return error

def getMagnetism(int1, int2, energies, reg):
    # Getting error expectation for each data point    
    int1_error = getError(int1, energies, delta_E, delta_int)
    int2_error = getError(int2, energies, delta_E, delta_int)
    
    # Putting together the XMCD spectrum    
    s = []
    error = []
    for n, e in enumerate(energies):
        new = int1[n] - int2[n]
        err = math.sqrt(int1_error[n]**2 + int2_error[n]**2)
        s.append(new)
        error.append(err)
        
    # Read data
    A, B, A_err, B_err = 0, 0, 0, 0
    for n, e in enumerate(energies):
        if e > reg[0] and e < reg[1]:
            A += s[n]
            A_err += error[n]**2
        elif e > reg[1] and e < reg[2]:
            B += s[n]
            B_err += error[n]**2
    A_err = math.sqrt(A_err)
    B_err = math.sqrt(B_err)

    A = abs(A)
    B = abs(B)
    
    C = 10.0         # beV
    # Calculate m_s
    m_s = 1/C*(2*B - A)
    m_s_err = 1/C* math.sqrt((2*B_err)**2 + A_err**2)

    # Calculate m_o
    m_o = -2/3/C*(A + B)
    m_o_err = -2/3/C * math.sqrt(A_err**2 + B_err**2)
    
    print "m_s has a value of "+str(m_s)+" +- "+str(100*abs(m_s_err/m_s))+"%"
    print "m_o has a value of "+str(m_o)+" +- "+str(100*abs(m_o_err/m_o))+"%"
    print "Ratio: "+str(m_o/m_s)+" +- "+str(100*math.sqrt((m_s_err/m_s)**2 + (m_o_err/m_o)**2))+" %"

    return m_o, m_o_err, m_s, m_s_err

if __name__ == '__main__':
    print "Gathering measurements data"
    
    with h5py.File(FILE_LC, 'r') as f:
        energies = f["/entry1/smpmamp39/pgmenergy"]
        
        spectrum = f["/entry1/smpmamp39/smpmamp39"]
        i_0 = f["/entry1/instrument/sm5amp8/sm5amp8"]

        LCEnergies = [e for e in energies]
        LCSpectrum = [s/i_0[n] for n, s in enumerate(spectrum)]

    with h5py.File(FILE_LCflip, 'r') as f:
        energies = f["/entry1/smpmamp39/pgmenergy"]
        
        spectrum = f["/entry1/smpmamp39/smpmamp39"]
        i_0 = f["/entry1/instrument/sm5amp8/sm5amp8"]

        LCflipEnergies = [e for e in energies]
        LCflipSpectrum = [s/i_0[n] for n, s in enumerate(spectrum)]

    with h5py.File(FILE_RC, 'r') as f:
        energies = f["/entry1/smpmamp39/pgmenergy"]
        spectrum = f["/entry1/smpmamp39/smpmamp39"]
        
        i_0 = f["/entry1/instrument/sm5amp8/sm5amp8"]

        RCEnergies = [e for e in energies]
        RCSpectrum = [s/i_0[n] for n, s in enumerate(spectrum)]

    with h5py.File(FILE_RCflip, 'r') as f:
        energies = f["/entry1/smpmamp39/pgmenergy"]
        spectrum = f["/entry1/smpmamp39/smpmamp39"]
        
        i_0 = f["/entry1/instrument/sm5amp8/sm5amp8"]

        RCflipEnergies = [e for e in energies]
        RCflipSpectrum = [s/i_0[n] for n, s in enumerate(spectrum)]

    with h5py.File(FILE_RC_RT, 'r') as f:
        energies = f["/entry1/smpmamp39/pgmenergy"]
        spectrum = f["/entry1/smpmamp39/smpmamp39"]
        
        i_0 = f["/entry1/instrument/sm5amp8/sm5amp8"]

        RCEnergies_RT = [e for e in energies]
        RCSpectrum_RT = [s/i_0[n] for n, s in enumerate(spectrum)]

    with h5py.File(FILE_LC_RT, 'r') as f:
        energies = f["/entry1/smpmamp39/pgmenergy"]
        spectrum = f["/entry1/smpmamp39/smpmamp39"]
        
        i_0 = f["/entry1/instrument/sm5amp8/sm5amp8"]

        LCEnergies_RT = [e for e in energies]
        LCSpectrum_RT = [s/i_0[n] for n, s in enumerate(spectrum)]

    energies = np.linspace(LCEnergies[0], LCEnergies[-1], num=len(LCEnergies)*5, endpoint=True)

    print "============"
    print "First background"
    print "============"
    LCNorm, LCSpectrum_norm, reg = getBackground(LCEnergies, LCSpectrum)
    splineRep = interpolate.splrep(LCEnergies, LCSpectrum_norm, s=smooth)
    LCSpectrum_norm = interpolate.splev(energies, splineRep, der=0)
    print "============"
    print "Second background"
    print "============"
    LCflipNorm, LCflipSpectrum_norm, reg = getBackground(LCflipEnergies, LCflipSpectrum)
    splineRep = interpolate.splrep(LCflipEnergies, LCflipSpectrum_norm, s=smooth)
    LCflipSpectrum_norm = interpolate.splev(energies, splineRep, der=0)
    print "============"
    print "Third background"
    print "============"
    RCNorm, RCSpectrum_norm, reg = getBackground(RCEnergies, RCSpectrum)
    splineRep = interpolate.splrep(RCEnergies, RCSpectrum_norm, s=smooth)
    RCSpectrum_norm = interpolate.splev(energies, splineRep, der=0)
    print "============"
    print "Fourth background"
    print "============"
    RCflipNorm, RCflipSpectrum_norm, reg = getBackground(RCflipEnergies, RCflipSpectrum)
    splineRep = interpolate.splrep(RCflipEnergies, RCflipSpectrum_norm, s=smooth)
    RCflipSpectrum_norm = interpolate.splev(energies, splineRep, der=0)
    print "============"
    print "Fifth background"
    print "============"
    LCNorm, LCSpectrum_RT_norm, reg = getBackground(LCEnergies_RT, LCSpectrum_RT)
    splineRep = interpolate.splrep(LCEnergies_RT, LCSpectrum_RT_norm, s=smooth)
    LCSpectrum_RT_norm = interpolate.splev(energies, splineRep, der=0)
    print "============"
    print "Sixth background"
    print "============"
    RCNorm, RCSpectrum_RT_norm, reg = getBackground(RCEnergies_RT, RCSpectrum_RT)
    splineRep = interpolate.splrep(RCEnergies_RT, RCSpectrum_RT_norm, s=smooth)
    RCSpectrum_RT_norm = interpolate.splev(energies, splineRep, der=0)

    # Get magnetism
##    m_o, m_o_err, m_s, m_s_err = getMagnetism(RCSpectrum_norm, LCSpectrum_norm, energies, reg)

    # Make XMCD
    XMCD = [(LCSpectrum_norm[n]-s)*4 for n, s in enumerate(RCSpectrum_norm)]      
    XMCDflip = [(LCflipSpectrum_norm[n]-s)*4 for n, s in enumerate(RCflipSpectrum_norm)]
    XMCD_RT = [(LCSpectrum_RT_norm[n]-s)*4 for n, s in enumerate(RCSpectrum_RT_norm)]

    # Save data to CSV
    spectra = [LCSpectrum_norm, LCflipSpectrum_norm, RCSpectrum_norm, RCflipSpectrum_norm, LCSpectrum_RT_norm, RCSpectrum_RT_norm, XMCD, XMCDflip, XMCD_RT]
    for n, s in enumerate(spectra):
        save(energies, s, str(n+1))

    # Plot data
    plt.plot(energies, [0 for e in energies], linewidth=2.0)
##    plt.plot(RCEnergies, LCSpectrum)
##    plt.plot(RCEnergies, LCNorm)
    plt.plot(energies, LCSpectrum_norm, label='$\sigma^+H^+$', linewidth=2.0)
##    plt.plot(RCEnergies, LCflipSpectrum_norm, 'k--',label='$\sigma^+H^-$')
    plt.plot(energies, RCSpectrum_norm, label='$\sigma^-H^+$', linewidth=2.0)
##    plt.plot(RCEnergies, RCflipSpectrum_norm, 'k:', label='$\sigma^-H^-$')
    plt.plot(energies, XMCD,  linewidth=2.0, label='$\sigma^+-\sigma^-$ x4 in $H^+$')
    plt.plot(energies, XMCDflip, linewidth=2.0 , label='$\sigma^+-\sigma^-$ x4 in $H^-$')
##    plt.plot(energies, XMCD_RT,  linewidth=2.0, label='$\sigma^+-\sigma^-$ x4 in $H^+$, RT')
    
    legend = plt.legend(loc='upper right', shadow=True, fontsize=16)
##    title = plt.title("Ti L XMCD spectrum of ex-situ STO|LMO(7) at different temperatures")
    plt.xlabel("Photon energy (eV)", fontsize=18)
    plt.ylabel("Drain current intensity (arb. units)", fontsize=18)
##    plt.text(energies[reg[1]]-5,min(XMCD), "m_s = "+str(m_s)+" +- "+str(abs(100*m_s_err/m_s))+"%\n"+
##                                    "m_o = "+str(m_o)+" +- "+str(abs(100*m_o_err/m_o))+"%\n"+
##                                    "Ratio = "+str(m_o/m_s)+" +- "+str(100*math.sqrt((m_s_err/m_s)**2 + (m_o_err/m_o)**2))+" %",
##             fontsize=10)
##    
    plt.show()   

