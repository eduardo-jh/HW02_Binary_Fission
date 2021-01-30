#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW2 - Problem 7. Bacteria growth, difference calculator
https://mathinsight.org/bacteria_growth_initial_model_exercises Exercise 5

Created on Thu Jan 21 12:11:46 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm  # allows linear regression without intercept

# Bacteria density data
# B = [0.028, 0.047, 0.082, 0.141, 0.240, 0.381]  # Exercise 5
B = np.array([0.022, 0.036, 0.060, 0.101, 0.169, 0.266])  # From problem 2

steps = len(B)  # Adjust the length of the vectors with the number of steps
dB = np.zeros(steps)
dt = 16  # time interval

t = np.linspace(0, (steps-1)*dt, steps)  # actual time vector

for i in range(1, steps):
    dB[i] = B[i] - B[i-1]  # compute the increment between time steps

# Perform a linear regression with B and dB, then plot (dB vs B)
model = sm.OLS(dB, B)  # No intercept by default, force through the origin
results = model.fit()
slope = results.params[0]  # grow rate of population
print("slope=", slope)

# Figure 1, plotting dB vs B
plt.figure(1)
plt.plot(B, dB, 'kx', B, slope*B, 'b-')
plt.legend(['data', 'linear regression $R^2$=%.3f' % results.rsquared], loc='best')
plt.xlabel('B')
plt.ylabel('dB')
plt.savefig('p7_bacteria_linear2.png', dpi=300, bbox_inches='tight')

# Generate an exponential equation ('exact solution')
tdouble = np.log(2)/np.log(1+slope)*dt
print('tdouble =', tdouble)
K = np.log(2)/tdouble
Bexp = B[0] * np.exp(K*t)

# Make 'predictions' using the analytical solution to the linear dynamical system,
# (also an exponential equation) in the form B(t) = B[0]*R^t with R>1
# we don't need to know the previous value, each calculation is only dependant of the time 't'
Bmodel = B[0]*pow(slope+1, t/dt)
print("The population after %d steps is: %.3f" % (steps, Bmodel[-1]))

# Figure 2, plotting P (from data and model) vs t
plt.figure(2)
plt.plot(t, B, 'bx', t, Bmodel, 'r-', t, Bexp, 'k+')
plt.legend(['data',
            'numerical B=%g$\cdot$(1+%.4f)$^t$' % (B[0], slope),
            'exact B=%g$\cdot$exp(%.4f$\cdot$t)' % (B[0], K)],
           loc='best')
plt.xlabel('Time (minutes)')
plt.ylabel('Bacteria population')
plt.savefig('p7_bacteria_%dsteps2.png' % steps, dpi=300, bbox_inches='tight')
