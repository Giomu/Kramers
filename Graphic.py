#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:59:52 2019

@author: giorgio
"""

''' Graphical part of Kramers.py '''



from Kramers import a, b, eps, gamma, N, num_sim, r1, r2, y_min1, y_min2, p1_teor, p2_teor
import numpy as np
from Kram_Functions import V, Boltz
import matplotlib.pyplot as plt





x = np.loadtxt('Posizioni.txt')
p = np.loadtxt('Impulsi.txt')
t = np.loadtxt('Time.txt')
frac_sx = np.loadtxt('LeftFraction.txt')
frac_dx = np.loadtxt('RightFraction.txt')




KT = Boltz(eps, gamma)
bi = 40   #number of bins





''' Build histograms of positions '''

 


#theoretical rho of positions

sub_int=10000

M = np.linspace(-3, 3, sub_int+1)

rho_teor = np.exp(-V(M, a, b)/KT) #theoretical distribution not normalized

amp_int = abs(M[1] - M[0])        #interval lenght

trap = 0.5*amp_int*(rho_teor[0] + rho_teor[-1])   #Trapezi method of normalization
trap += sum(rho_teor[1:-1]) * amp_int             #normalization Trapezi

plt.plot(M, rho_teor[:]/trap, color='black', label='Theoretical')



#sperimental rho of positions

for s in range(100, N, 150):
    
    
    num_x = [0]*(bi+1)
    
    xmax=np.max(x[:,s])
    xmin=np.min(x[:,s])
    
    deltax= (xmax-xmin)/bi
    xvar = np.linspace(xmin, xmax, (bi+1))
    
    
    for z in range(0, num_sim):
        
        k = (x[z,s]-xmin)/deltax
        
        num_x[int(k)] += 1 
        
    
    distr = num_x[:]/(num_sim*deltax)
    
    plt.step(xvar[:], distr[:], label='T=%i' %s)



plt.title('Distribution of Positions')
plt.xlabel('Positions')
plt.ylabel('rho')
plt.grid(True)
plt.legend(loc='best')
plt.show()


#*********************##******************************







''' Build histogram of impulses '''



#theoretical rho of impulses


rho_p_teor = np.zeros(N)

f = np.linspace(-5, 5, N)

for h in range(1, N):
        
    rho_p_teor[h] = np.exp(-(f[h])*(f[h])/(2*KT))/np.sqrt(2*np.pi*KT)
    

    
#Numerical rho of impulses

for s in range(10, 150, 20):
    
    
    num_p = [0]*(bi+1)
    
    pmax=np.max(p[:,s])
    pmin=np.min(p[:,s])
    
    deltap= (pmax-pmin)/bi
    pvar = np.linspace(pmin, pmax, (bi+1))
   
    for z in range(0, num_sim):
        
        a = (p[z,s]-pmin)/deltap
        
        num_p[int(a)] += 1 
        

    plt.step(pvar[:], num_p[:]/(num_sim*deltap), label='T=%i' %s)   



plt.plot(f, rho_p_teor[:], color='black', label='Theoretical')
plt.title('Impulse Distribution')
plt.ylabel('rho')
plt.legend(loc='best')
plt.grid(True)
plt.show()



#****************************##***************************************





''' Changements of probabilities in Time '''



somma_var_sx = np.zeros(N)    #Inizializations
somma_var_dx = np.zeros(N)

p1_teor_var = np.zeros(N)
p2_teor_var = np.zeros(N)




for o in range(0, N):
    
    #Sperimental probabilities in time
    somma_var_sx[o] = np.sum(frac_sx[:,o])/num_sim
    somma_var_dx[o] = np.sum(frac_dx[:,o])/num_sim
    
    #PTheoretical probabilities in time
    p1_teor_var[o]  = np.exp(-t[o]*(r1+r2)) + r2*(1-np.exp(-t[o]*(r1+r2)))/(r1+r2)
    p2_teor_var[o]  = r1*(1-np.exp(-t[o]*(r1+r2)))/(r1+r2)
    






plt.plot(t[1:], somma_var_sx[1:], label='Left Numerical')
plt.plot(t[:], p1_teor_var[:],  label='Theoretical in time')
plt.title('Probability that is left in time')
plt.xlabel('Time')
plt.grid(True)
plt.legend(loc='best')
plt.show()



plt.plot(t[1:], somma_var_dx[1:], label='Right Numerical')
plt.plot(t[:], p2_teor_var[:],  label='Theoretical in time')
plt.title('Probability that is right in time ')
plt.xlabel('Time')
plt.grid(True)
plt.legend(loc='best')
plt.show()








#*************************#********************


''' Probabilities at End of Simulation '''



sinistra = np.zeros(num_sim)   #inizializations
destra   = np.zeros(num_sim)


##Popolations at final time    
for h in range(0 , num_sim):
    
    sinistra[h] = frac_sx[h,N-1]
    destra[h]   = frac_dx[h,N-1]



distribuzioni = [np.sum(sinistra), np.sum(destra)] #height of histogram

p1_sper = np.sum(sinistra)/num_sim  #sperimental prob sx at end simulation
p2_sper = np.sum(destra)/num_sim    #sperimental prob dx at end simulation

centri = [y_min1, y_min2]                          #histogram bases



plt.bar(centri, height = distribuzioni, width=1 ,color='b')
plt.xlabel('Positions')
plt.ylabel('Number of Particles')
plt.title('Particles ripartition at Simulations end')
plt.grid(True)
plt.text(2,170, 'p1_theor=%1.3f' %p1_teor, color = 'r')
plt.text(2,130, 'p2_theor=%1.3f' %p2_teor, color = 'r')
plt.text(2,90, 'p1_sper=%1.3f' %p1_sper, color = 'b')
plt.text(2,50, 'p2_sper=%1.3f' %p2_sper, color = 'b')
plt.show()






















