"""
XAS fitting code

Main characteristics of Mn:
    L_I     769.1 #eV
    L_II    649.9 #eV
    L_III   638.7 #eV

    'MN1': [0, 2.114214, 2.281614 , 3.072244 , 3.372561], #6S5/2, 8P5/2, 6D9/2, 4D7/2, 4P5/2
    'MN2': [0, 1.1745, 1.7761810, 3.4154, 3.703336], #7S3, 5S2, 5D4, 5G6, 3P2
    'MN3': [0, 3.325802,  3.616334,  4.005595, 4.85701] #6S5/2, 4G11/2, 4P5/2, 4D7/2, 2I11/2

Main characteristics of Ti:
    Ti L I	560.9 #eV
    Ti L II	460.2 #eV
    Ti L III	453.8 #eV

Main characteristics of Co:
    Co L I	925.1 #eV
    Co L II	793.2 #eV
    Co L III	778.1 #eV



Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2014). NIST
Atomic Spectra Database (ver. 5.2), [Online]. Available:
http://physics.nist.gov/asd [2016, March 14]. National Institute of
Standards and Technology, Gaithersburg, MD.
"""
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import math
import copy
import random
import h5py
#FILE = "data/i09-75022.nxs"
#FILE = "data/i09-75371.nxs"
FILE = "data/i09-75372.nxs"
SHIFT = 0.1

def getProportion(spectra, intensity, energy):
    
    l = [0.1*i for i in range(10)]
    bestScore = float("inf")
    print 'checking proportion'
    for e in l:
        for r in l:
            for t in l:
                if e + r + t == 1:
                    thisSpect = []
                    for n, s in enumerate(spectra[0]):
                        thisSpect.append(s*e + spectra[1][n]*r + spectra[2][n]*t)
                    score = getScore(intensity, thisSpect)
                    if score < bestScore:
                        bestScore = score
                        proportion = (e, r, t)
                        spect = thisSpect

    print "Best proportion:",proportion
    
    return spect

def getData(filename):
    # data gathered
    data=open(filename).read()
    energy, intensity = [], []
    
    for line in data.split("\n"):
        this = line.split(";")
        if len(this) < 2: break
        e, i = this
        energy.append(float(e.replace(",",".")))
        intensity.append(float(i.replace(",",".")))

    return energy, intensity

def normalise(energy, i, E_L3):
    # Hand-waivy normalisation function 
    norm = []
##    for e in energy:
##        if e < 635:
##            a = -3.33*10**(-4)
##            b = 0.211
##            norm.append(.128 + a*e + b)
##        elif e < E_L3:
##            norm.append(.128)
##        else: norm.append(i[-1])

    A, B = E_L3, 0.8
    for e in energy:
        # fermi-dirac distribution-like normalisation function
        n = i[1] + (i[-1]-i[1])*1/(math.exp((A-e)/B)+1)
        norm.append(n-0.015)
        
    return norm, [intensity - norm[n] for n, intensity in enumerate(i)]

def normMN(energy, intensity, start, edge3, edge2):
    # Hand-waivy normalisation function
    E_L2 = 649.9 +SHIFT #eV
    E_L3 = 638.7 +SHIFT #eV
    
    norm = []
##    A, B = E_L3, 0.8
##    for e in energy:
##        # fermi-dirac distribution-like normalisation function
##        n = start + (edge2-start)*1/(math.exp((A-e)/B)+1)
##        norm.append(n)
    for e in energy:
        q = start
        if e > E_L3:
            q += edge3
        if e > E_L2:
            q += edge2
        norm.append(q)

    return norm, [i - norm[n] for n, i in enumerate(intensity)]

def scaleIntensity(intensities, energy1, energy2, mx):
    rtn = []
    for int_l in intensities:
        intensity = []
        ratio = mx/max(int_l)
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
                new = a*e1 + b
                if new < 0: new = 0
                intensity.append(new * ratio * 1.15)
            except: intensity.append(0)
        rtn.append(intensity)

    return rtn

def makeFit(energy, edges, ion, main_intensity):
    """
    An attempt to fit the absorbtion spectrum is made

    Level splitting for ground state (only relevant levels accounted):
        7S3: 0 eV
        5S2: 1.1745 eV
        5G6: 3.4154 eV
    """
    spect = [0.0 for e in energy]

    E = [638.7+SHIFT, 649.9+SHIFT] #eV, L3 -> L2
    IONS = {
        'MN1': [0, 2.114214, 2.281614 , 3.072244 , 3.372561], #6S5/2, 8P5/2, 6D9/2, 4D7/2, 4P5/2
        'MN2': [0, 1.1745, 1.7761810, 3.4154, 3.703336, 3.7058299, 3.784458, 3.907326, 4.065181, 4.109792 ,4.3099508,4.497480, 4.688076, 4.7568065 , 4.800700], #7S3, 5S2, 5D4, 5G6, 3P2
        'MN3': [0, 3.325802,  3.616334,  4.005595, 4.85701] #6S5/2, 4G11/2, 4P5/2, 4D7/2, 2I11/2
        }
    E_levels = IONS[ion]

    for n, e in enumerate(energy):
        for k, i in enumerate(edges):
            edge = E[k]
            for q, l in enumerate(E_levels):
                poep = math.exp(-(e-edge-l-SHIFT)**2/i[q][1])
                spect[n] += main_intensity*i[q][0]*poep
                
    return spect

def getScore(intensity, spect):
    s = 0
    for n, i in enumerate(intensity):
        diff = i - spect[n]
        s += diff**2
    return s

def hillClimb(energy, intensity, ion, iterations=100):
    """
    Hillclimbing algorithm.
    1 < intensity < 15, in steps of 1
    0.1 < lifetime < 2.5, in steps of 0.1
    """
    print "\nStarting hillblimber..."
    
    main_intensity = max(intensity)
    main = [[[0.21999999999999986, 0.26], [0.6800000000000002, 0.1], [0.3799999999999999, 1.3], [0.030000000000000027, 2.9000000000000012], [-0.01, 0.7999999999999999], [-0.01, 0.5], [-0.01, 0.4], [-0.01, 0.5], [-0.00999999999999999, 1.7000000000000004], [-0.00999999999999999, 2.3000000000000007], [0.09000000000000001, 3.0000000000000013], [0.0, 0.20000000000000004], [-0.01, 1.7000000000000004], [0.009999999999999997, 0.6], [0.04, 0.7999999999999999]], [[0.060000000000000005, 0.6], [-0.010000000000000004, 3.0000000000000013], [0.10999999999999999, 2.8000000000000007], [0.019999999999999997, 2.800000000000001], [-3.469446951953614e-18, 2.9000000000000012], [-0.010000000000000004, 2.3000000000000007], [-0.01, 0.7999999999999999], [-0.01, 1.0999999999999999], [-3.469446951953614e-18, 3.0000000000000013], [0.009999999999999997, 2.800000000000001], [-0.010000000000000004, 3.0000000000000013], [0.0, 1.2], [-0.01, 1.6000000000000003], [-0.01, 1.8000000000000005], [0.02, 1.7000000000000004]]]
    spect = makeFit(energy, main, ion, main_intensity)
    bestScore = getScore(intensity, spect)

    u = 0
    for i in range(iterations):
        newMain = copy.deepcopy(main)

        tup = random.choice(random.choice(newMain))
        
        if random.choice([True, False]):
            if random.choice([True, False]) and tup[0] < 1:
                tup[0] += 0.01
            elif random.choice([True, False]) and tup[0] > -0.005:
                tup[0] -= 0.01
        elif random.choice([True, False]):
            if random.choice([True, False]) and tup[1] < 3:
                tup[1] += 0.1
            elif random.choice([True, False]) and tup[1] > 0.25:
                tup[1] -= 0.1
        else: continue

        spect = makeFit(energy, newMain, ion, main_intensity)
        score = getScore(intensity, spect)
        try:
            evaluate = math.exp((bestScore-score)*50*i/3.0)
        except:
            if bestScore-score > 0: evaluate = 1
            else: evaluate = 0

        p = random.random()
        if evaluate > p:
            bestScore = score
            main = newMain
            if evaluate < 1: u += 1

        if i%50 == 0 and i != 0:
            print bestScore
            print "Annealing:",u
            u = 0
            
    print main    
    return spect

def do(energy, intensity):
    # Building fits with hillblimber
    spectra = []
    for i in ['MN1', 'MN2', 'MN3']:
        spect = hillClimb(energy, intensity_norm, i, 1500)
        spectra.append(spect)

    # Getting percentages
    bestSpect = getProportion(spectra, intensity_norm, energy)

    # Plotting
    plt.plot(energy, bestSpect)
    line, = plt.plot(energy, intensity_norm)
    plt.setp(line, linewidth=2)    
    line, = plt.plot(energy, [0.14+s for s in spectra[0]])
    plt.setp(line, linestyle='--')
    line, = plt.plot(energy, [0.07+s for s in spectra[1]])
    plt.setp(line, linewidth=2.5)
    line, = plt.plot(energy, [0.21+s for s in spectra[2]])
    plt.setp(line, linestyle='--')
    plt.plot(energy, [0 for e in energy])
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Energy (eV)')

    # Legend
    measurement = plt.Line2D([], [], color='blue',
                              markersize=10, label='Measurement data')
    MN1 = plt.Line2D([], [], color='green',
                              markersize=10, label='Mn simulation')
    MN2 = plt.Line2D([], [], color='red',
                              markersize=10, label='Mn 1+ simulation')
    MN3 = plt.Line2D([], [], color='cyan',
                              markersize=10, label='Mn 2+ simulation')
    plt.legend(handles=[measurement, MN1, MN2, MN3])
    plt.title("Manganese XAS spectrum, both measured and simulated.")

    # Show legend
    plt.show()

if __name__ == '__main__':
    print "Gathering measurements data"
    #Get data
    with h5py.File(FILE, 'r') as f:
        enrg = f["/entry1/sm5amp8/pgmenergy"]
        spect = f["/entry1/smpmamp39/smpmamp39"]
        i_0 = f["/entry1/sm5amp8/sm5amp8"]
        blacklist = []
        
        energy = [e for n, e in enumerate(enrg) if n not in blacklist]
        intensity = [s/i_0[n] for n, s in enumerate(spect) if n not in blacklist]

    # Normalisation determined
    norm, intensity_norm = normalise(energy, intensity, E_L3 = 638.7+SHIFT)

    # Manganese things
    f = open("Mn-valences/energies.csv").read()
    E_literature = [float(e) for e in f.split(",")]
    f = open("Mn-valences/Mn(II).csv").read()
    MNII = [float(i) for i in f.split(",")]
    f = open("Mn-valences/Mn(III).csv").read()
    MNIII = [float(i) for i in f.split(",")]
    f = open("Mn-valences/Mn(IV).csv").read()
    MNIV = [float(i) for i in f.split(",")]

    # Normalise
    MNII_norm, MNII = normMN(E_literature, MNII, 1.5, 3, 5)
    MNIII_norm, MNIII = normMN(E_literature, MNIII, 1.5, 14.5, 3)
    MNIV_norm, MNIV = normMN(E_literature, MNIV, 2, 15.5, 0)

    # Rescale
    MNII, MNIII, MNIV = scaleIntensity([MNII, MNIII, MNIV], energy, E_literature, max(intensity_norm))


    spect = hillClimb(energy, MNII, 'MN2', 10000)
    # Get optimal proportions
    #spect = getProportion([MNII, MNIII, MNIV], intensity_norm, energy)

    plt.plot(energy, [0 for e in energy])
    plt.plot(energy, MNII)
    plt.plot(energy, spect)
    plt.show()

