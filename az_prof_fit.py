import numpy as np
import math as m
import sys as sys
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

# prompts user for input                       #
# checks that the file exists in the directory #
#                                              #
# input data file must be in format:           #
# -rows:    radii                              #
# -columns: azimuthal divisions                #

while True:
    try:
        print ''
        prof_file = raw_input('Enter file name: ')
        print ''
        f = open(prof_file, 'r')
        break
    except IOError:
        print 'Error: no such file in this directory.'
        print ''
    except TypeError:
        print 'Error: no such file in this directory.'
        print ''

lines = f.readlines()
n = len(lines)
f.close()
raddiv = n

# reads in the data into a [raddiv,angdiv] array #

for line in lines:
    columns = line.split()
    numbers = [float(n) for n in columns]
    angdiv = len(numbers)

prof = np.empty((raddiv,angdiv))

dum = 0
for line in lines:
    columns = line.split()
    numbers = [float(n) for n in columns]
    for j in range(angdiv):
        prof[dum,j] = numbers[j]
    dum = dum+1

# creates the phi array #

phi = np.linspace(0, 2*m.pi, num=angdiv)
phi = phi*(180./m.pi)

#####################################################################################
# plots the data and waits for user input                                           #
#                                                                                   #
# button press inputs:                                                              #
# -'g' => begins prof_fit routine                                                   #
# -'n' => goes to next profile                                                      #
# -'r' => goes to specified radius                                                  #
# -'x' => exits plotting loop                                                       #
#####################################################################################

def fit():
    global phi, prof, rad

    prof_fit(phi, prof[rad,:], rad)

def goto():
    global rad

    while True:
        try:
            rad = raw_input('Go to which radius? ')
            print ''
            rad = int(rad)
            break
        except ValueError:
            print 'Error: must be an integer.'
            print ''

def on_key_event(event):
    global rad, phi, prof, plot_num

    key = event.key
    sys.stdout.flush()
            
    if key in 'xngr':
        pass
    else:
        print 'Error: not a recognized input (xngr).'
        print ''

    if key == 'x':
        plt.close('all')

    if key == 'n':
        rad = rad + 1
        
        if rad < raddiv:
            print rad
            print ''

            prof_plot()
        else:
            plt.close('all')

    if key == 'g':
        fit()

    if key == 'r':
        goto()

        if rad < raddiv:
            rad = rad - 1
        else:
            print 'Error: radius must be less than ' + str(raddiv)
            print ''
            goto()
        
        prof_plot()
        
def prof_plot():
    global phi, rad, prof, plot_num

    plot_num = plot_num + 1
    fig = plt.figure(plot_num)
    plt.close(plot_num - 1)
    fig.canvas.mpl_connect('key_press_event', on_key_event)
    plt.plot(phi, prof[rad,:], 'black')
    plt.show(plot_num)

def rad_in():
    global rad

    while True:
        try:
            rad = raw_input('Begin at which radius? ')
            rad = int(rad)
            break
        except ValueError:
            print 'Error: must be an integer.'
            print ''

##############################################
# fits gaussians to user specified range     #
# plots the fit over the profile             #
#                                            #
# At start:                                  #
# -'y' => continue routine                   #
# -'n' => leave routine, return to prof_plot #
#                                            #
# After fit:                                 #
# -'y' => redo fit                           #
# -'n' => keep fit, return to prof_plot      #
#                                            #
##############################################

def gauss_fit():
    popt, pcov = curve_fit(gauss, x_fit, y_fit, p0)
    return popt, pcov

def prof_fit(x,y,n):
    print ''
    print '!!!!!!!!!!!!!!!'
    print '! Profile Fit !'
    print '!!!!!!!!!!!!!!!'
    print ''

    # checks that the inputs are numbers       #
    # checks that the first number is smaller  #

    while True:
        try:
            x1_c = raw_input('Enter starting angle: ')
            x2_c = raw_input('Enter ending angle: ')
            print ''

            x1 = float(x1_c)
            x2 = float(x2_c)
        except ValueError:
            print 'Error: one or both inputs are not numbers.'
            print ''
        except RuntimeError:
            print 'Something went wrong. Exiting.'
            print ''
            return

        try:
            test1 = x1
            test2 = x2
        except UnboundLocalError:
            print 'Inputs can only be numbers. Exiting.'
            print ''
            return

        if x1 > x2:
            print 'Error: first angle must be the smaller angle.'
        else:
            break

    try:
        cen_guess = (x1 + x2)/2.
    except UnboundLocalError:
        print 'Make sure inputs are ONLY numbers.'
        print ''
        return

    # allocates the x and y arrays for fitting #

    dum = len(x)
    x_fit = []
    y_fit = []

    for i in range(dum):
        if x[i] >= x1 and x[i] <= x2:
            x_fit.append(x[i])
            y_fit.append(y[i])

    # the function being fit #

    def gauss(x_in, cen, sig, a, c):
        return a*(np.exp(-(x_in - cen)**2./(2.*sig**2.))) + c

    # fit the range with the specified function #

    p0 = [cen_guess,1.,20,1]
    try:
        popt, pcov = curve_fit(gauss, x_fit, y_fit, p0)
    except RuntimeError:
        print 'Fit failed. Try a different fitting range.'
        print ''
        return
    except ValueError:
        print 'Fit failed. Try a different fitting range.'
        print ''
        return
    except TypeError:
        print 'Fit failed. Try a larger fitting range.'
        print ''
        return

    print ''

    try:
        centroid = round(popt[0], 2)
        cen_err = round(np.sqrt(pcov[0,0]), 2)
        sigma = round(popt[1], 2)
        sig_err = round(np.sqrt(pcov[1,1]), 2)
        normal = round(popt[2], 2)
        norm_err = round(np.sqrt(pcov[2,2]), 2)
        offset = round(popt[3], 2)
        off_err = round(np.sqrt(pcov[3,3]), 2)
        max_val = round(max(y_fit),2)
    except TypeError:
        print 'Something went wrong. Try fitting with different range.'
        print ''
        return

    print "The centroid is " + str(centroid) + " +/- " + str(cen_err) + " degrees."
    print "The sigma is " + str(sigma) + " +/- " + str(sig_err) + "."
    print "The normalization is " + str(normal) + "."
    print "The offset is " + str(offset) + " +/- " + str(off_err) + " counts."
    print "The peak is " + str(max_val) + " counts."
    print ''

    # plots the fit in a new figure                  #
    # waits for input to redo or return to prof_plot #
    
    fit_plot = plot_num + 1
    fig2 = plt.figure(fit_plot)
    fig2.canvas.mpl_connect('key_press_event', on_key_event)
    plt.plot(x,y, 'black')
    plt.plot(x,gauss(x,*popt), 'r')
    plt.show(fit_plot)

    while True:
        redo = raw_input('Redo? (yn): ')
        if redo in 'yn':
            print ''
            break
        else:
            print 'Error: not a recognized input (yn).'
            print ''

    if redo == 'y':
        plt.close(fit_plot)
        fit()
    if redo == 'n':
        plt.close(fit_plot)
        return

#######################
# begins the plotting #
#######################

# takes in the radius to start plotting                 
# checks that the input is an integer                   
# checks that the radius is smaller than the max radius 

rad_in()

while True:
    if rad < raddiv:
        rad = rad - 1
        break
    else:
        print 'Error: radius must be less than ' + str(raddiv)
        print ''
        rad_in()

print ''
plot_num = 1
fig = plt.figure(plot_num)
prof_plot()
