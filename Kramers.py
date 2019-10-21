#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 11:55:30 2019

@author: giorgio
"""

''' Kramers Project Software and Computing '''




from Kram_Functions import Boltz, V, Calc_min_V, Prob_fuga, Prob_fin
import random
import numpy as np
import logging
import re
import sys
from scipy.misc import derivative
#import matplotlib.pyplot as plt
#import seaborn as sns






#print of presentation 
print('\n\n')   
print('''          ******Welcome to Kramers Simulation******\n
          Thankyou for having filled params.txt with the Values required for the Simulation!\n
          At the end of the integration please run Graphic.py to plot the Numerical Results obtained
          \n          *****************************************''')



'''Define parameters of the Simulation'''


N       = 800    #Numero passi di integrazione
dt      = 0.1    #time Step 


f = open('params.txt') # Open file on read mode
lines = f.read().split("\n") # Create a list containing all lines
f.close()


EPS     = re.search(r"[+-]?[0-9]+\.[0-9]+", lines[11])
GAMMA   = re.search(r"[+-]?[0-9]+\.[0-9]+", lines[12])
NUM_SIM = re.search(r"\d+", lines[13])
A       = re.search(r"[+-]?[0-9]+\.[0-9]+", lines[14])
B       = re.search(r"[+-]?[0-9]+\.[0-9]+", lines[15])
S0      = re.search(r"[+-]?[0-9]+\.[0-9]+", lines[16])






while True:
    try:
        eps = float(EPS.group())
    
    except ValueError:
        logging.error("Sorry, epsilon should be a float digit.")
        sys.exit()
    
    if eps < 0.3 or eps>1:
        logging.warning("epsilon should be in range (0.3, 1) to see some effects")
        sys.exit()
    
    else:
        break





while True:
    try:
        gamma = float(GAMMA.group())
    
    except ValueError:
        logging.error("gamma should be a float digit.")
        sys.exit()
    
    if gamma < 0 or gamma>1:
        logging.warning("gamma should be in range (0, 1)")
        sys.exit()
    
    else:
        break





while True:
    try:
        num_sim = int(NUM_SIM.group(0))
    
    except ValueError:
        logging.error(" N should be an int digit.")
        sys.exit()
    
    if num_sim < 1000:
        logging.warning(" The error goes as 1/sqrt(N) insert N higher then 1000")
        sys.exit()
    
    else:
        break





while True:
    try:
        a = float(A.group())
    
    except ValueError:
        logging.error("Sorry, a should be a float digit.")
        sys.exit()
    
    if a < 0:
        logging.warning("a should be greater then 0")
        sys.exit()
    
    else:
        break





while True:
    try:
        b = float(B.group())
    
    except ValueError:
        logging.error("Sorry, b should be a float digit.")
        sys.exit()
    
    if a < 0:
        logging.warning("b should be greater then 0")
        sys.exit()
    
    else:
        break




while True:
    try:
        s0 = float(S0.group())
    
    except ValueError:
        logging.error("Sorry, seed should be a float digit.")
        sys.exit()
    
    if s0 < 0:
        logging.warning("seed should be greater then 0")
        sys.exit()
    
    else:
        break



print('\n          **************')
print('          Parameters of Simulation:\n')
print('\n          --> eps:  ', eps)
print('          --> gamma:', gamma)
print('          --> N:    ', num_sim)
print('\n          --> a:  ', a)
print('          --> b:  ', b)
print('\n          --> seed: ', s0)
print('          **************\n\n')



KT = Boltz(eps, gamma)

y_min1, y_min2 = Calc_min_V(a,b)

r1, r2 = Prob_fuga(y_min1, y_min2, KT, a, b)

p1_teor, p2_teor = Prob_fin(N, r1, r2)






'''Build white noise with gaussian distribution of mean 0 and variance 1
    with the seed from params.txt'''


random.seed(s0)
mu = 0
sigma = 1

csi = np.zeros((num_sim, N))

for t in range (0, num_sim):
    for k in range(0, N):
        csi[t,k] = random.gauss(mu, sigma)
    



        





''' Integer the motion for different simulations '''

x0 = y_min1          #Initial position in left minimum
p0 = 0               #Initial impulse equal to 0

    
x = np.zeros((num_sim, N))   #Matrices containing different particles in the rows
p = np.zeros((num_sim, N))   #and their time evolutions in colomns
U = np.zeros((num_sim, N))

frac_sx = np.zeros((num_sim, N))  #fraction of particles in the left side
frac_dx = np.zeros((num_sim, N))  #fraction of particles in the right side

t = np.zeros(N)      #time


for i in range(0, num_sim):
    
    x[i,0] = x0
    
    p[i,0] = p0

    U[i,0] = V(x[i,0], a, b)
    

    for j in range(1, N):
    
        V_diff = derivative(V, x[i,j-1], dx=1e-6, args=(a,b))
        
            
        p[i,j] = p[i,j-1] - gamma*p[i,j-1]*dt - V_diff*dt + eps*csi[i,j]*pow(dt, 0.5)
        x[i,j] = x[i,j-1] + p[i,j]*dt
            
        U[i,j] = V(x[i,j], a, b)
        
        t[j] = j*dt

        if x[i,j]<0:
            frac_sx[i,j] = 1
            frac_dx[i,j] = 0
            
        elif x[i,j]>0:
            frac_sx[i,j] = 0
            frac_dx[i,j] = 1
            
    
        

    


np.savetxt('Posizioni.txt', x[:,:], header='Positions of single Particles')
np.savetxt('Impulsi.txt',   p[:,:], header='Impulses of single Particles')
np.savetxt('LeftFraction.txt',  frac_sx[:,:], header='Particles Fraction on the left side')
np.savetxt('RightFraction.txt', frac_dx[:,:], header='Particles Fraction on the right side')
np.savetxt('Time.txt', t[:], header='Time steps')




print('''\n\n          ******End of Simulation!******\n\n          If you want to plot Numerical results please run: Graphic.py\n
          Else visit the following txt files for other manipulations:\n
          --> Posizioni.txt
          --> Impulsi.txt
          --> LeftFraction.txt
          --> RightFraction.txt
          --> Time.txt\n''')






