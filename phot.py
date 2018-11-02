import numpy as np
import matplotlib.pyplot as plt

def main():
    ReadData()
    #PlotIntenProf()
    PlotSurfBright()

def ReadData():
    global n, r, inten, inten_err, surf, cen_surf, scale_len, ax1, fig

    while True:
        try:
            in_file = raw_input('Enter file name: ')
            print ''
            f = open(in_file, 'r')
            break
        except IOError:
            print 'Error: no such file in this directory. Try again.'
            print ''
        except TypeError:
            print 'Error: no such file in this directory. Try again.'
            print ''

    lines = f.readlines()
    n = len(lines)
    f.close()

    r_col = 0
    inten_col = 1
    inten_err_col = 2

    r = np.empty(shape=0)
    inten = np.empty(shape=0)
    inten_err = np.empty(shape=0)

    for line in lines:
        p = line.split()
        r = np.append(r, float(p[r_col]))
        inten = np.append(inten, float(p[inten_col]))
        inten_err = np.append(inten_err, float(p[inten_err_col]))

def PlotIntenProf():
    global n, r, inten, inten_err, surf, cen_surf, scale_len, ax1, fig

    fig = plt.figure(1)
    ax1 = fig.add_subplot(211)

    while True:
        ax1.scatter(r,inten)
        fig.show()
        fig.canvas.draw()

        end = raw_input('To quit enter q: ')
        if not isinstance(end,basestring):
            print 'Enter a letter.'
        elif end == 'q':
            break

def FitPhot():
    global n, r, inten, inten_err, surf, cen_surf, scale_len, ax1, fig

    while True:
        fit_start = GetValue("fit starting value")
        fit_end = GetValue("fit ending value")

        fit_r_d = []
        fit_surf_d = []

        for i in range(n):
            if r[i] >= fit_start and r[i] <= fit_end:
                fit_r_d.append(r[i])
                fit_surf_d.append(surf[i])
        fit_r = np.array(fit_r_d)
        fit_surf = np.array(fit_surf_d)

        A = np.vstack([fit_r, np.ones(len(fit_r))]).T
        try:
            b, m = np.linalg.lstsq(A,fit_surf)[0]
        except ValueError:
            return
        except TypeError:
            return

        cen_surf = m
        scale_len = 1.086/b

        fit_y = cen_surf + 1.086*fit_r/scale_len
        #print fit_r
        #print fit_y
        #print fit_surf

        ax1.plot(fit_r,fit_y,'r')
        #plt.gca().invert_yaxis()
        fig.show()
        fig.canvas.draw()

        while True:
            end = raw_input("Good fit? ")
            if end == 'y':
                inc = GetValue("inclination")
                if inc == 'q':
                    return
                inc = inc*np.pi/180.0
                cen_surf += 2.5*np.log10(1/np.cos(inc))
                mag = cen_surf - 2.5*np.log10(2*np.pi*scale_len*scale_len) - 2.5*np.log10(np.cos(inc))
                print ''
                print "The central surface brightness is " + str(round(cen_surf,2)) + " mag arcsec^-2"
                print "The disk scale length is " + str(round(scale_len,2)) + " arcsec."
                print "The total magnitude is " + str(round(mag,2)) + " mag"
                print ''
                return
            else:
                break


def GetValue(val_in):
    while True:
        val = raw_input("Enter " + val_in + " (q to quit): ")
        if val == 'q':
            return 'q'
            break
        try:
            val = float(val)
            return val
        except ValueError:
            print 'Make sure input is a number.'
            print ''

def PlotSurfBright():
    global n, r, inten, inten_err, surf, cen_surf, scale_len, ax1, fig

    numArgs = 4
    phot = np.empty(numArgs)
    words = []
    words.append("zero point")
    words.append("exposure time")
    words.append("plate scale")
    words.append("extinction")

    for i in range(numArgs):
        try:
            phot[i] = float(GetValue(words[i]))
            cont = True
        except ValueError:
            cont = False
            break

    if cont:
        print ''
        r = r*phot[2]

        surf = -2.5*np.log10(inten/(phot[1]*phot[2]*phot[2])) + phot[0] - phot[3]
        #print surf

        fig = plt.figure(1)
        ax1 = fig.add_subplot(211)
        
        ax1.scatter(r,surf)
        plt.gca().invert_yaxis()
        fig.show()
        fig.canvas.draw()

        while True:
            end = raw_input('To quit enter q, to fit phot enter f: ')
            if not isinstance(end,basestring):
                print 'Enter a letter.'
            elif end == 'q':
                break
            elif end == 'f':
                FitPhot()

###################################################

main()
