from scipy import misc
import matplotlib.pyplot as plt
import numpy as np
import cv2

f = cv2.imread('124,0.TIF')
mov = cv2.VideoCapture('practice.mov')

def average(cap = mov, x = (0, 570), y = (0, 570)):
    i = 0
    intens = []
    while cap.grab():
        i += 1
        if i%20 == 0:
            print "Frame: ",i
        ret, frame = cap.read()
        this = clip(frame, 285, 315, 345, 365)
        intens.append(getIntensity(this))
        if i == 1:
            img = frame
        else:
            img = cv2.addWeighted(img,(i-1)/float(i),frame,1/float(i),0)

##    plt.imshow(frame)
##    plt.show()
##    plt.imshow(img)
##    plt.show()
##    print "Intensity: ",getIntensity(img)
##    new = clip(img, 285, 315, 345, 365)
##    plt.imshow(new)
##    plt.show()
##    print "Intensity: ",getIntensity(new)
##    plt.plot(intens)
##    plt.show()
    cap.release()
    cv2.destroyAllWindows()

    print "destroyed!"
    height, width = 200, 200
    video = cv2.VideoWriter('video4.avi',-1,15,(width,height))

    print "hallo"
    for i in range(150):
        print i+1
        this = clip(img,100+i, 300+i, 100, 300)
        video.write(this)

    print "hallo2"
    cv2.destroyAllWindows()
    video.release()
    
    return img

def reduceScale(img=f, reduction=6):

    #reduce colums
    main = []
    for i in img:
        sub = []
        l = []
        for n, e in enumerate(i): 
            if n%reduction == 0 and n != 0:
                sub.append(np.mean(l))
                l = []
            l.append(e)
        main.append(sub)
    main = np.array(main)

    #reduce rows
    new = []
    l = [0 for e in range(len(main[0]))]
    for n, i in enumerate(main):
        if n%reduction == 0 and n != 0:
            l = l/float(reduction)
            new.append(l)
            l = [0 for e in range(len(l))]
        l = np.add(l, i)
        
    return np.array(new)

def clip(img, xmin, xmax, ymin, ymax):
    l = []
    for n, i in enumerate(img):
        if n >= ymin and n < ymax:
            l.append(i[xmin:xmax])
    return np.array(l)

def getIntensity(img=f):
    l = []
    for i in img:
        l.append(np.mean(i))
    return np.mean(l)   

def plot(img):
    plt.imshow(img)
    plt.show()

if __name__ == '__main__':
    height, width = 200, 200
    intensity = 0
    video = cv2.VideoWriter('video2.avi',-1,4,(350,350))
    for n in range(751):
        for d in range(10):
            FILE = str(n)+','+str(d)+".TIF"
            try:
                f = cv2.imread(FILE)
                f = clip(f, 180, 530, 20, 370)
##                this = clip(f, 40, 50, 250, 260)
##                if intensity == 0:
##                    intensity = getIntensity(this)
##                ratio = getIntensity(this)/intensity
##                print "before",getIntensity(f)
##                new = []
##                for line in f:
##                    newl = []
##                    for l in line:
##                        newl.append(ratio*l)
##                    new.append(newl)
##                print "after",getIntensity(new)
                video.write(f)
                print "frame added"
            except: pass

    cv2.destroyAllWindows()
    video.release()
    print "Video ready"



