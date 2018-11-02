import numpy as np
import math as m
import sys as sys
import matplotlib.pyplot as plt

#########################################
# change these based on individual runs #
#########################################

raddiv = int(sys.argv[2])
angdiv = int(sys.argv[1])
print(raddiv, angdiv)

####################
# R and Phi arrays #
####################

rstart = int(sys.argv[3])
rmax = int(sys.argv[4])
print(rstart, rmax)
rad = np.linspace(rstart, rmax+1, num=raddiv)
#print(rad)
rad_trim = np.delete(rad,-1)

phi = np.linspace(0, 2*m.pi, num=angdiv)
phi = phi*(180./m.pi)

####################
# read in the data #
####################

inname1 = sys.argv[5]
inname2 = sys.argv[6]
print(inname1, inname2)

f1 = open(inname1, 'r')
f2 = open(inname2, 'r')
lines1 = f1.readlines()
lines2 = f2.readlines()
n1 = len(lines1)
n2 = len(lines2)
f1.close()
f2.close()

ang_prof = []
rad_prof = []

for line in lines1:
    columns1 = line.split()
    numbers1 = [float(n) for n in columns1]
    ang_prof.append(numbers1)
for line in lines2:
    columns2 = line.split()
    numbers2 = [float(n) for n in columns2]
    rad_prof.append(numbers2)

#################
# plot the data #
#################

switch = int(sys.argv[7])
print(switch)
plt.figure()

if switch == 1:
    for i in range(raddiv):
        plt.plot(phi,ang_prof[i], 'k')

elif switch == 2:
    for i in range(angdiv):
        if i & 2 == 0:
            plt.plot(rad_trim,rad_prof[i])

#plt.legend()
plt.ylabel('Intensity (ADU)')
plt.xlabel('Azimuthal Angle (degrees)')
axes = plt.gca()
axes.set_xlim([-10, 370])
plt.show()
