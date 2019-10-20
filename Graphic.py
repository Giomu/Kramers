#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:59:52 2019

@author: giorgio
"""

''' Parte grafica di Kramers.py '''



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
bi = 40   #numero di bins





''' Creiamo gli Istogrammi delle posizioni '''

 


#rho teorica delle posizioni

sub_int=10000

M = np.linspace(-3, 3, sub_int+1)

rho_teor = np.exp(-V(M, a, b)/KT) #distribuzione teoricaa non normalizzata

amp_int = abs(M[1] - M[0])        #ampiezza dell'intervallino

trap = 0.5*amp_int*(rho_teor[0] + rho_teor[-1])   #normalizzazione col metodo dei Trapezi
trap += sum(rho_teor[1:-1]) * amp_int             #normalizzazione Trapezi

plt.plot(M, rho_teor[:]/trap, color='black', label='Teorica')



#rho sperimentale delle posizioni

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



plt.title('Distribuzione delle particelle')
plt.xlabel('Posizione')
plt.ylabel('rho')
plt.grid(True)
plt.legend(loc='best')
plt.show()


#*********************##******************************







''' Creo gli istogrammi degli impulsi '''



#rho teorica degli impulsi


rho_p_teor = np.zeros(N)

f = np.linspace(-5, 5, N)

for h in range(1, N):
        
    rho_p_teor[h] = np.exp(-(f[h])*(f[h])/(2*KT))/np.sqrt(2*np.pi*KT)
    

    
#rho numerica degli impulsi

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



plt.plot(f, rho_p_teor[:], color='black', label='Teorica')
plt.title('Distribuzione degli Impulsi')
plt.ylabel('rho')
plt.legend(loc='best')
plt.grid(True)
plt.show()



#****************************##***************************************





''' Andamenti delle probabilità al variare del tempo '''



somma_var_sx = np.zeros(N)    #inizializzazioni
somma_var_dx = np.zeros(N)

p1_teor_var = np.zeros(N)
p2_teor_var = np.zeros(N)




for o in range(0, N):
    
    #Probabilità sperimentali nel tempo
    somma_var_sx[o] = np.sum(frac_sx[:,o])/num_sim
    somma_var_dx[o] = np.sum(frac_dx[:,o])/num_sim
    
    #Probabilità teoriche nel tempo
    p1_teor_var[o]  = np.exp(-t[o]*(r1+r2)) + r2*(1-np.exp(-t[o]*(r1+r2)))/(r1+r2)
    p2_teor_var[o]  = r1*(1-np.exp(-t[o]*(r1+r2)))/(r1+r2)
    






plt.plot(t[1:], somma_var_sx[1:], label='Numerica sinistra')
plt.plot(t[:], p1_teor_var[:],  label='Teorica nel tempo')
plt.title('Probabilità che sia a sinistra nel tempo')
plt.xlabel('Tempo')
plt.grid(True)
plt.legend(loc='best')
plt.show()



plt.plot(t[1:], somma_var_dx[1:], label='Numerica destra')
plt.plot(t[:], p2_teor_var[:],  label='Teorica nel tempo')
plt.title('Probabilità che sia a destra nel tempo ')
plt.xlabel('Tempo')
plt.grid(True)
plt.legend(loc='best')
plt.show()








#*************************#********************


''' Andamenti delle probabilità a fine simulazione '''



sinistra = np.zeros(num_sim)   #inizializzazioni
destra   = np.zeros(num_sim)


##popolazioni ai tempi finali    
for h in range(0 , num_sim):
    
    sinistra[h] = frac_sx[h,N-1]
    destra[h]   = frac_dx[h,N-1]



distribuzioni = [np.sum(sinistra), np.sum(destra)] #altezza dell'istogramma

p1_sper = np.sum(sinistra)/num_sim  #prob sperimentale a fine simulazione sx
p2_sper = np.sum(destra)/num_sim    #prob sperimentale a fine simulazione dx

centri = [y_min1, y_min2]                          #base dell'istogramma



plt.bar(centri, height = distribuzioni, width=1 ,color='b')
plt.xlabel('Positions')
plt.ylabel('Number of Particles')
plt.title('Particles ripartition at end Simulation')
plt.grid(True)
plt.text(2,170, 'p1_teor=%1.3f' %p1_teor, color = 'r')
plt.text(2,130, 'p2_teor=%1.3f' %p2_teor, color = 'r')
plt.text(2,90, 'p1_sper=%1.3f' %p1_sper, color = 'b')
plt.text(2,50, 'p2_sper=%1.3f' %p2_sper, color = 'b')
plt.show()






















