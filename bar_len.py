import numpy as np
import matplotlib.pyplot as plt

def main():
    global cont
    setParams()
    readData()
    if cont:
        barLen()

##############
# SET PARAMS #
##############

def setParams():
    global pix, n, disk_start, r, pa, pa_err, ecc, ecc_err, pa_avg, pa_avg_err, ecc_avg, ecc_avg_err, inc

    pix = 0.228

#############
# READ DATA #
#############

def readData():
    global n, r, ecc, ecc_err, pa, pa_err, cont

    while True:
        try:
            in_file = raw_input('Enter file name (q to quit): ')
            print ''
            if in_file == 'q':
                cont = False
                break
            f = open(in_file, 'r')
            cont = True
            break
        except IOError:
            print 'Error: no such file in this directory. Try again.'
            print ''
        except TypeError:
            print 'Error: no such file in this directory. Try again.'
            print ''

    if cont:
        lines = f.readlines()
        n = len(lines)
        f.close()

        r_col = 0
        ecc_col = 5
        ecc_err_col = 6
        pa_col = 7
        pa_err_col = 8

        r = np.empty(shape=0)
        ecc = np.empty(shape=0)
        ecc_err = np.empty(shape=0)
        pa = np.empty(shape=0)
        pa_err = np.empty(shape=0)

        for line in lines:
                p = line.split()
                r = np.append(r, float(p[r_col]))
                ecc = np.append(ecc, float(p[ecc_col]))
                ecc_err = np.append(ecc_err, float(p[ecc_err_col]))
                pa = np.append(pa, float(p[pa_col]))
                pa_err = np.append(pa_err, float(p[pa_err_col]))

###########
# BAR LEN #
###########

def barLen():
    global ecc, ecc_err, r, pix, n

    r = r*pix

    fig = plt.figure(1)
    ax1 = fig.add_subplot(211)

    ecc_max = 0
    pa_max_ecc = 0
    bar_pa = 0

    ax1.plot(r,ecc)
    fig.show()
    fig.canvas.draw()

    ecc_ind = 0
    pa_ind = 0

    
    while True:
        r_start = raw_input("Enter init fitting radius: ")
        try:
            r_start = float(r_start)
        except ValueError:
                print 'Make sure input is a number.'
                print ''

        r_end = raw_input("Enter max fitting radius: ")
        try:
                r_end = float(r_end)
        except ValueError:
                print 'Make sure input is a number.'
                print ''

        if isinstance(r_start,float) or isinstance(r_start,int):
                if isinstance(r_end, float) or isinstance(r_end,int):
                        break
        else:
                print 'Make sure input is a number.'
                print ''

    for i in range(n):
        if r[i] >= r_start and r[i] <= r_end:
            if ecc[i] > ecc_max:
                ecc_max = ecc[i]

    for i in range(n-1):
        if r[i] >= r_start and r[i] <= r_end:
            if ecc[i] == ecc_max:
                bar_maxe = r[i]
                pa_max_ecc = pa[i]
                ecc_ind = i

    for i in range(n):
        if r[i] <= r_end and r[i] >= bar_maxe:
            if abs(pa[i] - pa_max_ecc) > 5:
                bar_pa = r[i]
                pa_ind = i
                break

    ecc_err = (abs(r[ecc_ind] - r[ecc_ind-1]) + abs(r[ecc_ind]-r[ecc_ind+1]))/2.
    pa_err = (abs(r[pa_ind] - r[pa_ind-1]) + abs(r[pa_ind]-r[pa_ind+1]))/2.

    x1 = [bar_maxe, bar_maxe]
    y = [0, 1]
    x2 = [bar_pa, bar_pa]

    while True:
        ax1.plot(x1,y,'r')
        ax1.plot(x2,y,'g')
        fig.show()
        fig.canvas.draw()

        print "The max ecc bar length is " + str(round(bar_maxe,2)) + " +/- " + str(round(ecc_err,2)) + " arcsec."
        print "The PA discontin. bar length is " + str(round(bar_pa,2)) + " +/-" + str(round(pa_err,2)) + " arcsec."
        print ''

        bar_maxe = bar_maxe/pix
        bar_pa = bar_pa/pix

        print "The max ecc bar length is " + str(round(bar_maxe,2)) + " pixels."
        print "The PA discontin. bar length is " + str(round(bar_pa,2)) + " pixels."

        end = raw_input('To end, enter e: ')
        if not isinstance(end,basestring):
            print 'Enter a string.'
        elif end == 'e':
            break

################################################################
main()
