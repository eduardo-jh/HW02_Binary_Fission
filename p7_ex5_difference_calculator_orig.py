#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW2 - Problem 7. Bacteria growth, difference calculator
https://mathinsight.org/bacteria_growth_initial_model_exercises Exercise 5

Created on Thu Jan 21 12:11:46 2021
@author: eduardo
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# Bacteria density data
# B = np.array([0.028, 0.047, 0.082, 0.141, 0.240, 0.381])  # Exercise 5
B = np.array([0.022, 0.036, 0.060, 0.101, 0.169, 0.266])  # From problem 2

steps = len(B)  # Adjust the length of vectors to the number of steps
dB = np.zeros(steps)
Bgraph = np.zeros(steps)
dt = 16  # time interval

t = np.linspace(0, (steps-1)*dt, steps)  # actual time vector

for i in range(1, steps):
    dB[i] = B[i] - B[i-1]  # compute the increment between time steps
    Bgraph[i] = B[i-1]
print(B, Bgraph)

# Perform a linear regression with the data
slope, intercept, r_value, p_value, std_err = stats.linregress(Bgraph, dB)

# Figure 1, plotting dB vs B
plt.figure(1)
plt.plot(Bgraph, dB, 'bx', Bgraph, slope*Bgraph, 'r-')  # plot and linear eq.
plt.legend(['data', 'linear regression $R^2$=%.3f' % r_value**2], loc='best')
plt.xlabel('B')
plt.ylabel('dB')
plt.savefig('p7_bacteria_linear.png', dpi=300, bbox_inches='tight')
print('slope =', slope, 'intercept =', intercept)

# Generate an exponential equation ('exact solution')
tdouble = np.log(2)/np.log(1+slope)*dt
K = np.log(2)/tdouble
Bexp = B[0] * np.exp(K*t)
print('tdouble =', tdouble, 'K =', K)

# Make 'predictions' using the analytical solution to the linear dynamical system,
# (also an exponential equation) in the form B(t) = B[0]*R^t with R>1
# we don't need to know the previous value, each calculation is only dependant of the time 't'
Bmodel = B[0]*pow(slope+1, t/dt)
print("The population after %d steps is: %.3f" % (steps, Bmodel[-1]))

# Figure 2, plotting B vs t
plt.figure(2)
plt.plot(t, B, 'bx', t, Bmodel, 'r-', t, Bexp, 'k+')  # Plot data vs exponential growth eq.
plt.legend(['data',
            'numerical B=%g$\cdot$(1+%.4f)$^t$' % (B[0], slope),
            'exact B=%g$\cdot$exp(%.4f$\cdot$t)' % (B[0], K)],
           loc='best')
plt.xlabel('Time (minutes)')
plt.ylabel('Bacteria population')
plt.savefig('p7_bacteria_%dsteps.png' % steps, dpi=300, bbox_inches='tight')
