from scipy import misc
import matplotlib.pyplot as plt
import numpy as np
import cv2

MnII = cv2.imread('Mn(II)-fromMnO.bmp')
MnIII = cv2.imread('Mn(III)-fromLaMnO_3.bmp')
MnIV = cv2.imread('Mn(IV)-fromSrMnO_3.bmp')

def imgToData(img):
    nrPoints = len(img)
    data = img.transpose()[0]
    spectrum = [0 for q in data]
    for n, line in enumerate(data):
        l = []
        for i, point in enumerate(line):
            if point == 0:
                l.append(nrPoints-i)
        spectrum[n] = np.average(l)

    return spectrum

# Convert data
MnIISpect = imgToData(MnII)
MnIIISpect = imgToData(MnIII)
MnIVSpect = imgToData(MnIV)

def frange(start, stop, length):
    i = start
    l = []
    step = (stop-start)/float(length)
    for n in range(length):
        l.append(i)
        i += step
    return l

energy = frange(637, 662, 390)

print MnIISpect
print MnIIISpect
print MnIVSpect
print energy

# Plot data
plt.plot(energy, MnIISpect)
plt.plot(energy, MnIIISpect)
plt.plot(energy, MnIVSpect)
plt.show()
