import h5py
import numpy as np
import matplotlib.pyplot as plt
import math


FILES = ["i09-75371.nxs"]
# script to get all files
##for i in range(110):
##    if i < 10:
##        FILES.append("i09-7500"+str(i)+".nxs")
##    elif i < 100:
##        FILES.append("i09-750"+str(i)+".nxs")
##    else:
##        FILES.append("i09-75"+str(i)+".nxs")
        
FOLDER = "/entry1"
MAXLEN = 50



def getData(f, folder, sets):
    print
    print "\t",folder
    
    theseSets = []
    for item in f[folder].items():
        if isinstance(item[1], h5py._hl.dataset.Dataset):
            theseSets.append(item[0])
    print len(theseSets), "datasets: ",
    for s in theseSets: print s,",",
    print
    
    for item in f[folder].items():
        if isinstance(item[1], h5py._hl.dataset.Dataset):
            if item[1].shape[0] == 1:
                if not isinstance(item[1].value[0], unicode):
                    print item[0], ": ",type(item[1].value[0])
                    
                else:
                    if len(item[1].value[0]) < MAXLEN:
                        print item[0], ": ",item[1].value[0]
                    else:
                        print item[0], ": ",item[1].value[0][:MAXLEN],"(...)"
            else:
                sets.append(item)
        else:
            sets = getData(f, item[1].name, sets)
    
    return sets

if __name__ == "__main__":
    for FILE in FILES:
##        try:
        with h5py.File(FILE, 'r') as f:
            print "=========================================="
            print "Data in",FILE
            print "=========================================="
            sets = []
            getData(f, "/entry1", sets)

            print
            for item in sets:
                print item[0], "in shape",item[1].shape

            energies = f["/entry1/default/pgmenergy"]
            spectrum = f["/entry1/default/sm5iamp8"]
            spectrum2 = f["/entry1/default/smpmiamp39"]
            print energies
##            plt.plot(energies, spectrum)
            plt.plot(energies, spectrum2)
            plt.show()
##        except: pass


##    print image

##    # Scale things 
##    bottom = f["/entry1/SXVB_ARPES/excitation_energy"][0]
##    top = bottom + 451 * f["/entry1/SXVB_ARPES/energy_step"][0]

##    # ARPES data (energy scale, angle scale, image)
##    angles = f["/entry1/SXVB_ARPES/angles"][0]
##    energies = f["/entry1/SXVB_ARPES/energies"][0]
##    image = f["/entry1/SXVB_ARPES/image_data"][0]
##    # No idear
##    external_io = f["/entry1/SXVB_ARPES/external_io_data"][0]
##    # Hard x-ray intensity spectrum, plot against energies-array
##    spectrum = f["/entry1/SXVB_ARPES/spectrum_data"][0]
    
def getkvalue(ang):
    m = 5.1*10**5 #eV/c^2
    hbar = 6.58*10**(-16) #eVs
    c = 3.0*10**8 #m/s
    a = 0.6 * 10**(-10) #m

    hv = 5.5*10**3 #eV
    phi = 4.6 #eV
    E_bind = 460.15 #eV

    E_kin = hv - phi - E_bind #eV
    
        
    theta = 2*math.pi*ang/360
    k = math.sqrt(2*m*E_kin)/(hbar*c)* math.sin(theta) * a / math.pi#1/m
    return k


##image = image.T
##plt.imshow(
##    image, origin="lower", aspect='auto', extent=(
##        getkvalue(angles.min()),getkvalue(angles.max()),bottom,top
##        )
##    )
##plt.show()

##plt.plot(energies, [math.log(s) for s in spectrum])
##plt.plot(energies, spectrum)
##
##plt.show()
