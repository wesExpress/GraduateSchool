from astropy.io import fits
from astropy import wcs
import numpy as np
import scipy as sp
import math as m
import sys as sys
import os
import os.path
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import Normalize

def Main():
    SetGlobals()
    FileRead()
    GetPhot()
    PhotConv()
    GetColor()

#############
# PLOTCOLOR #
#############

def PlotColor(color):
    global length
    global titleName 

    testing = False

    cdict = {'red':  ((0.0, 0.0, 0.0),   
                      (0.5, 0.5, 0.5),   
                      (1.0, 0.8, 0.8)),  

             'green': ((0.0, 0.0, 0.0),   
                       (0.5, 0.2, 0.2),   
                       (1.0, 0.0, 0.0)),  

             'blue':  ((0.0, 0.8, 0.8),   
                       (0.5, 0.5, 0.5),   
                       (1.0, 0.0, 0.0))  
       }

    plt.figure(1)
    if testing:
        norm = Normalize(vmin=0,vmax=1,clip=True)
        plt.imshow(used,cmap='Greys',norm=norm,interpolation='nearest')
        plt.xlim(0,length)
        plt.ylim(0,length)
    else:
        norm = Normalize(vmin=0.0,vmax=2.0,clip=True)
        plt.imshow(color,cmap='Greys',norm=norm,interpolation='nearest')
        plt.xlim(length/2 - length/4, length/2 + length/4)
        plt.ylim(length/2 - length/4, length/2 + length/4)

    cbar = plt.colorbar()
    cbar.ax.set_ylabel("B-I")

    plt.rcParams.update({'font.size':22})

    plt.xlabel("Pixels")
    plt.ylabel("Pixels")
    plt.title(titleName)
    plt.show(block=False)

############
# GETCOLOR #
############

def GetColor():
    global coordsB, coordsI
    global pa, inc
    global surfB, surfI
    global maxR, barLen
    global badPix
    global used
    global titleName

    xCen = length/2
    yCen = length/2
    
    intenValB = surfB[yCen,xCen]

    xCenB = xCen
    yCenB = yCen

    searchRange = int(length*0.05)

    difs = np.zeros((searchRange,searchRange))

    # find center
    for i in range(searchRange): 
        y = yCen -searchRange/2 + i
        for j in range(searchRange):
            x = xCen -searchRange/2 + j
            if surfB[y,x] < intenValB:
                intenValB = surfB[y,x]
                xCenB = x
                yCenB = y

    worldB = wcs.utils.pixel_to_skycoord(xCenB,yCenB,coordsB,origin=0)
    guessI = wcs.utils.pixel_to_skycoord(xCenB,yCenB,coordsI,origin=0)

    newSearchRange = searchRange
    oldCheck = worldB.separation(guessI)

    for y in range(int(newSearchRange)):
        yCoord = yCen - int(newSearchRange)/2 + y
        for x in range(int(newSearchRange)):
            xCoord = xCen - int(newSearchRange)/2 + x

            testI = wcs.utils.pixel_to_skycoord(xCoord,yCoord,coordsI,origin=0)
            newCheck = worldB.separation(testI)
            difs[y,x] = newCheck.hour
            if newCheck.hour < oldCheck.hour:
                print "new minimum separation: " +  str(newCheck)
                oldCheck = newCheck
                worldI = testI
                #print worldB
                #print worldI

    yCenI = worldI.to_pixel(coordsI,origin=0)[1]
    xCenI = worldI.to_pixel(coordsI,origin=0)[0]

    xDifB = int(round(xCenB - xCen))
    yDifB = int(round(yCenB - yCen))
    xDifI = int(round(xCenI - xCen))
    yDifI = int(round(yCenI - yCen))

    surfB = np.roll(surfB, -yDifB, axis=0)
    surfB = np.roll(surfB, -xDifB, axis=1)
    surfI = np.roll(surfI, -yDifI, axis=0)
    surfI = np.roll(surfI, -xDifI, axis=1)

    vertMove = False
    horizMove = False

    titleName = raw_input("Enter plot title: ")

    while True:
        print ' '
        print "Image centers: (x) " + str(xCen) + ", (y) " + str(yCen)
        print "B centers: (x) " + str(xCenB) + ", (y) " + str(yCenB)
        print "I centers: (y) " + str(int(xCenI)) + ", (y) " + str(int(yCenI))
        print ' '

        print "B offsets: " + "(x) " + str(xDifB) + ", (y) " + str(yDifB)
        print "I offsets: " + "(x) " + str(xDifI) + ", (y) " + str(yDifI)

        if vertMove:
            surfI = np.roll(surfI, yDifI, axis=0)
            vertMove = False
        elif horizMove:
            surfI = np.roll(surfI, xDifI, axis=1)
            horizMove = False

        pa = pa*m.pi/180.
        inc = inc*m.pi/180.

        barColor = 0
        barCounter = 0
    
        discColor = 0
        discCounter = 0

        areaColor = 0
        areaCounter = 0

        image = surfB - surfI
        used = np.zeros((length,length))
        for y in range(length):
            for x in range(length):
                if surfB[y,x] != badPix and surfI[y,x] != badPix:
                    newX = ((x-xCen)*np.cos(pa) + (y-yCen)*np.sin(pa))/np.cos(inc)
                    newY = (y-yCen)*np.cos(pa) - (x-xCen)*np.sin(pa)

                    cond = newX*newX + newY*newY
            
                    # bar color
                    if cond <= barLen*barLen:
                        barColor += image[y,x]
                        barCounter += 1
                    
                        #used[y,x] = 1
                    # area color
                    if cond <= maxR*maxR:
                        areaColor += image[y,x]
                        areaCounter += 1

                        #used[y,x] = 1
                        # disc color
                        if cond > barLen*barLen:
                            discColor += image[y,x]
                            discCounter += 1

                            #used[y,x] = 1
                        
        barColor = barColor/barCounter
        discColor = discColor/discCounter
        areaColor = areaColor/areaCounter
            
        print ' '
        print "Bar color: " + str(round(barColor,4))
        print "Disc color: " + str(round(discColor,4))
        print "Area color: " + str(round(areaColor,4))
        print ' '
        
        PlotColor(image)

        while True:
            cont = raw_input("Good? (y) or (n): ")
            if cont == 'y':
                sys.exit()
            elif cont == 'n':
                while True:
                    print "Up (w) down (s) left (a) right (d): "
                    key = raw_input()
                    if key == 'w':
                        yDifI = 1
                        vertMove = True
                        break
                    elif key == 's':
                        yDifI = -1
                        vertMove = True
                        break
                    elif key == 'a':
                        xDifI = 1
                        horizMove = True
                        break
                    elif key == 'd':
                        xDifI = -1
                        horizMove = True
                        break
                    else:
                        print "Unrecognize input. Try again."
                        print ""
                    break
            else:
                print "Unrecognized unput. Try again."
                print ""
            plt.close(1)
            break

    imageName = "color.fits"
    if os.path.isfile(imageName):
        os.remove(imageName)
    fits.writeto(imageName,image)

############
# PHOTCONV #
############

def PhotConv():
    global imageDataB, imageDataI
    global surfB, surfI
    global length
    global zeroB, zeroI, exptB, exptI, pix, extB, extI
    global badPix

    surfB = np.zeros((length,length))
    surfI = np.zeros((length,length))

    for y in range(length):
        for x in range(length):
            meanB = 0
            meanI = 0

            meanRange = 5
            if y > meanRange and x > meanRange and y < length - meanRange and x < length - meanRange:
                for i in range(-2,2):
                    for j in range(-2,2):
                        meanB += imageDataB[y + i, x + j]
                        meanI += imageDataI[y + i, x + j]

            meanB = meanB/(meanRange*meanRange)
            meanI = meanI/(meanRange*meanRange)

            if imageDataB[y,x] > 0 and imageDataI[y,x] > 0 and imageDataB[y,x] < 3*meanB and imageDataI[y,x] < 3*meanI:
                surfB[y,x] = -2.5*np.log10(imageDataB[y,x]/(exptB*pix**2.)) + zeroB - extB
                surfI[y,x] = -2.5*np.log10(imageDataI[y,x]/(exptI*pix**2.)) + zeroI - extI
            elif meanB == 0:
                surfB[y,x] = badPix
            elif meanI == 0:
                surfI[y,x] = badPix
            else:
                surfB[y,x] = badPix
                surfI[y,x] = badPix

###########
# GETPHOT #
###########

def GetPhot():
    global zeroB, zeroI, exptB, exptI, pix, extB, extI, pa, inc, maxR, barLen

    params = np.loadtxt("calibs.txt")
    
    zeroB = params[0]
    zeroI = params[1]
    exptB = params[2]
    exptI = params[3]
    pix = params[4]
    extB = params[5]
    extI = params[6]
    maxR = params[7]
    pa = params[8]
    inc = params[9]
    barLen = params[10]

    pix = pix*2

    maxR = maxR/pix
    barLen = barLen/pix

    print ' '
    print "Bar length: " + str(round(barLen,4)) + " pixels."
    print "25.5 mag/arcsec^-2 radius: " + str(round(maxR,4)) + " pixels."
    print ' '

############
# FILEREAD #
############

def FileRead():
    global imageDataB, imageDataI, coordsB, coordsI
    global length

    while True:
        base = raw_input("Enter base name (q to quit): ")

        if base == 'q':
            sys.exit("Quitting.")

        try:
            dumB = fits.open(base + "_b.fits")
            dumI = fits.open(base + "_i.fits")
            break

        except IOError:
            print "No files matching this base. Try again."

    imageB = base + "_b_new.fits"
    imageI = base + "_i_new.fits"

    hduListB = fits.open(imageB)
    imageDataB = hduListB[0].data
    hduListI = fits.open(imageI)
    imageDataI = hduListI[0].data

    length = len(imageDataB)

    coordsB = wcs.WCS(hduListB[0].header)
    coordsI = wcs.WCS(hduListI[0].header)

    hduListB.close()
    hduListI.close()

##############
# SETGLOBALS #
##############

def SetGlobals():
    global length
    global imageDataB, imageDataI, coordsB, coordsI
    global surfB, surfI
    global zeroB, zeroI, exptB, exptI, pix, extB, extI, maxRad, inc, pa, barLen
    global badPix

    badPix = -5

#####################################################
Main()
