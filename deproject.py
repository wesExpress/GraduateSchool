import numpy as np
import matplotlib.pyplot as plt

def main():
    global cont
    setParams()
    readData()
    if cont:
        deproject()

##############
# SET PARAMS #
##############

def setParams():
    global pix, n, disk_start, r, pa, pa_err, ecc, ecc_err, x_0, y_0, pa_avg, pa_avg_err, ecc_avg, ecc_avg_err, inc

    pix = 0.228

#############
# READ DATA #
#############

def readData():
    global n, r, ecc, ecc_err, pa, pa_err, x_0, y_0, cont

    while True:
        try:
            in_file = raw_input('Enter file name ( q to quit): ')
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
        x_0_col = 9
        y_0_col = 10

        r = np.empty(shape=0)
        ecc = np.empty(shape=0)
        ecc_err = np.empty(shape=0)
        pa = np.empty(shape=0)
        pa_err = np.empty(shape=0)
        x_0 = np.empty(shape=0)
        y_0 = np.empty(shape=0)

        for line in lines:
                p = line.split()
                r = np.append(r, float(p[r_col]))
                ecc = np.append(ecc, float(p[ecc_col]))
                ecc_err = np.append(ecc_err, float(p[ecc_err_col]))
                pa = np.append(pa, float(p[pa_col]))
                pa_err = np.append(pa_err, float(p[pa_err_col]))
                x_0 = np.append(x_0, float(p[x_0_col]))
                y_0 = np.append(y_0, float(p[y_0_col]))  

        for i in range(n):
                if pa[i] < -90:
                        pa[i] = pa[i] + 180.
                elif pa[i] > 90:
                        pa[i] = pa[i] - 180.

#############
# DEPROJECT #
#############

def deproject():
    global pa_avg, pa_avg_err, ecc_avg, ecc_avg_err, inc, disk_start, r, pa, pa_err, ecc, ecc_err, pix
    
    r = r*pix

    fig = plt.figure(1)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    while True:
        ax1.scatter(r,pa)
        ax2.scatter(r,ecc)
        fig.show()
        fig.canvas.draw()

        while True:
            disk_start = raw_input("Starting disk range? ")
            try:
                disk_start = float(disk_start)
            except ValueError:
                print 'Make sure input is a number.'
                print ''

            disk_end = raw_input("Ending disk range? ")
            try:
                disk_end = float(disk_end)
            except ValueError:
                print 'Make sure unput is a number.'
                print ''

            print ''
            if isinstance(disk_start,float) or isinstance(disk_start,int):
                if isinstance(disk_end, float) or isinstance(disk_end,int):
                        break
            else:
                print 'Make sure input is a number.'
                print ''

        pa_avg = 0
        pa_avg_err = 0
        ecc_avg = 0
        ecc_avg_err = 0
    
        cnt = 0

        ind = 0
        for i in range(n):
            if r[i] >= disk_start:
                ind = i
                break

        for i in range(n):
            if r[i] >= disk_start and r[i] <= disk_end:
                if abs(pa[i] - pa[ind]) < 20:
                    pa_avg += pa[i]
                    pa_avg_err += pa_err[i]
                    ecc_avg += ecc[i]
                    ecc_avg_err += ecc_err[i]
                    cnt += 1

        bad = -10000.
        if cnt != 0: 
            pa_avg = pa_avg/cnt
            pa_avg_err = pa_avg_err/cnt
            ecc_avg = ecc_avg/cnt
            ecc_avg_err = ecc_avg_err/cnt
            inc_err = np.arccos(ecc_avg_err)
        else:
            pa_avg = bad
            pa_avg_err = bad
            ecc_avg = bad
            ecc_avg_err = bad
            inc_err = bad

        cond = pa_avg != bad and pa_avg_err != bad and ecc_avg != bad and ecc_avg_err != bad and inc_err != bad
        if cond:
            x = [disk_start, disk_end]
            y_pa = [pa_avg,pa_avg]
            y_ecc = [ecc_avg,ecc_avg]

            ax1.plot(x,y_pa,'r')
            ax2.plot(x,y_ecc,'r')
            fig.show()
            fig.canvas.draw()
    
            end = raw_input("Good range? ")
            if end == 'y':
                inc = np.arccos(1 - ecc_avg)
                inc = inc*180./np.pi

                print "The position angle is " + str(round(pa_avg,2)) + " +/- " + str(round(pa_avg_err,2)) + " degrees."
                print "The inclination is " + str(round(inc,2)) + " +/- " + str(round(inc_err,2)) + " degrees." 
                
                plt.close()

                break
            else:
                ax1.cla()
                ax2.cla()
        else:
            print "Bad fit"
            print ""

###################################################################
main()
