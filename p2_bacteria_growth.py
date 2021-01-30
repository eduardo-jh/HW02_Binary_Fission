#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW2 - Problem 2. Bacteria growth
https://mathinsight.org/bacteria_growth_initial_model

Created on Wed Jan 20 11:17:17 2021
@author: eduardo
"""
import matplotlib.pyplot as plt
import numpy as np

B_0 = 0.022  # initial bacteria population
dt = 16  # delta
N_t = 6  # Number of time steps
C = 2/3  # Relative growth rate, approx. slope of linear regression
R = 5/3  # The value of R in the exact solution
B = np.zeros(N_t)  # initialize the bacteria population
B[0] = B_0

t = np.linspace(0, (N_t-1)*dt, N_t)  # timestep array @ dt=16 minutes

# Use growth equation in 'function iteration form' to update the bacteria
# population, in this case we need to know the population in previus step
for n in range(len(t)-1):
    B[n+1] = B[n] + C*B[n]

# Make 'predictions' using the exponential equation, this is the analytical 
# solution to the linear dynamical system, in the form B(t) = B[0]*R^t with R>1
# we don't need to know the previous value, each calculation is only dependant of the time 't'
Bexp = B[0] * pow(R, t/dt)

# Plot the difference form solution (numerical) vs exponential (exact)
plt.plot(t, B, 'bo', t, Bexp, 'r-')
plt.legend(['numerical (fuc. iter. form)', 'exact (exponential sol.)'], loc='best')
plt.xlabel('Time (minutes)')
plt.ylabel('Bacteria population')
plt.savefig('p2_bacteria_%dsteps.png' % N_t, dpi=300, bbox_inches='tight')

print(t, B)
