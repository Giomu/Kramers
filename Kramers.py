#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 11:55:30 2019

@author: giorgio
"""

''' Kramers Project Software and Computing '''





import random
import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as plt
#import seaborn as sns









'''Definisco i parametri della simulazione'''

N       = 800                 #Numero di passi dell'integrazione
dt      = 0.1                 #Incrementino temporale

eps     = 0.4                 #modula l'effetto stocastico
gamma   = 0.1                 #Attrito
KT      = eps*eps/(2*gamma)   #kT (epsilon=1)

num_sim = 10000                #numero delle simulazioni da svolgere
                               #ovvero il numero di particelle



'''Calcolo il Rumore bianco con distribuzione gaussiana a 
media nulla e varianza 1'''

mu = 0
sigma = 1

csi = np.zeros((num_sim, N))

for t in range (0,num_sim):
    for k in range(0, N):
        csi[t,k] = random.gauss(mu, sigma)
    





'''Definisco il Potenziale'''

a = np.sqrt(2)    #Parametri del potenziale asimmetrico per individuare
b = np.sqrt(2)    #i punti dei due minimi delle due buche


#Definiamo l'espressione del potenziale
def V(y):
    return ((y*y-a*y)*(y*y+b*y))


#Definiamo la derivata seconda del potenziale
def V_II(y):
    return (12*y*y + 6*y*(b-a) - 2*a*b)



y_max = 0                                                 #max di V(y)

y_min1 = (1/8)*(3*(a-b) - np.sqrt(9*(a-b)*(a-b)+32*a*b))  #min di sx
y_min2 = (1/8)*(3*(a-b) + np.sqrt(9*(a-b)*(a-b)+32*a*b))  #min di dx
 
print('\nminimo di sinistra:',y_min1,'\nminimo di destra:',y_min2)





''' Calcoliamo le probabilità di fuga dalle due buche r1 e r2 '''

r1 = (np.sqrt(V_II(y_min1))/(2*np.pi))*np.exp((-V(0)+V(y_min1))/KT)

r2 = (np.sqrt(V_II(y_min2))/(2*np.pi))*np.exp((-V(0)+V(y_min2))/KT)





''' Le probabilità che al tempo finale la pallina si trovi a sx o dx '''

p1_teor = np.exp(-N*(r1+r2)) + r2*(1-np.exp(-N*(r1+r2)))/(r1+r2)

p2_teor = r1*(1-np.exp(-N*(r1+r2)))/(r1+r2)






''' Integro il moto per le diverse simulazioni '''

x0 = y_min1          #posizione iniziale nel minimo di sinistra
p0 = 0               #impulso iniziale fermo

    
x = np.zeros((num_sim, N))   #matrici contenenti le particelle nelle righe
p = np.zeros((num_sim, N))   #e gli istanti temporali nelle colonne
U = np.zeros((num_sim, N))

frac_sx = np.zeros((num_sim, N))   
frac_dx = np.zeros((num_sim, N))

t = np.zeros(N)


for i in range(0, num_sim):
    
    x[i,0] = x0
    
    p[i,0] = p0

    U[i,0] = V(x[i,0])
    

    for j in range(1, N):
    
        V_diff = derivative(V, x[i,j-1], dx=1e-6)
        
            
        p[i,j] = p[i,j-1] - gamma*p[i,j-1]*dt - V_diff*dt + eps*csi[i,j]*pow(dt, 0.5)
        x[i,j] = x[i,j-1] + p[i,j]*dt
            
        U[i,j] = V(x[i,j])
        
        t[j] = j*dt

        if x[i,j]<0:
            frac_sx[i,j] = 1
            frac_dx[i,j] = 0
            
        elif x[i,j]>0:
            frac_sx[i,j] = 0
            frac_dx[i,j] = 1
            
    
        

'''plt.plot(x[i,:], U[i,:])
plt.title('Potenziale')
plt.show() '''      


#np.savetxt('posizioni_asi.txt', x[:,:],newline='''\n\n\n\n\n\nPARTICELLA %j
#           \n\n\n\n\n''',header="""Posizioni delle singole particelle""")
#np.savetxt('impulsi_asi.txt', p[:,:], newline='''\n\n\n\n\n\nPARTICELLA %j
#           \n\n\n\n\n''',header='Impulsi delle singole particelle')
#np.savetxt('potenziali_asi.txt', U[:,:], newline='''\n\n\n\n\n\nPARTICELLA %j
#           \n\n\n\n\n''',header='Potenziali delle singole particelle')









''' Creiamo gli Istogrammi delle posizioni '''

bi = 40    


#rho teorica delle posizioni

sub_int=10000

M = np.linspace(-3, 3, sub_int+1)

rho_teor = np.exp(-V(M)/KT) #distribuzione teoricaa non normalizzata

amp_int = abs(M[1] - M[0]) #ampiezza dell'intervallino

trap = 0.5*amp_int*(rho_teor[0] + rho_teor[-1])   #normalizzazione
trap += sum(rho_teor[1:-1]) * amp_int             #normalizzazione

plt.plot(M, rho_teor[:]/trap, color='black', label='Teorica')



#rho sperimentale delle posizioni

for s in range(100, N, 150):
    
    
    num_x = [0]*(bi+1)
    
    xmax=np.max(x[:,s])
    xmin=np.min(x[:,s])
    
    deltax= (xmax-xmin)/bi
    mino = np.linspace(xmin, xmax, (bi+1))
    
    
    for z in range(0, num_sim):
        
        k = (x[z,s]-xmin)/deltax
        
        num_x[int(k)] += 1 
        

    #plt.plot(mino[:], num_x[:]/(num_sim*deltax), label='T=%i' %s)   
     
    ##da qui fino ai plot modifico le x e le y in modo che gli elementi degli
    ##array si duplichino e la y venga sfalzata in modo da costruire
    ##manualmente le successioni di puunti per gli istogrammi.
    ##la linea sopra serve nel caso non mi servisse l'isto ma solo la curva
    ##in tal caso cancella da qui fino ai plot e riabilita la riga sopra
    
    distr = num_x[:]/(num_sim*deltax)
    
    x_histo = [0]*2*len(mino)
    y_histo = [0]*2*len(distr)
    
    for c in range(1, (2*len(mino)), 2):
        
        mino = np.insert(mino, c, np.zeros(1))
        distr = np.insert(distr, c, np.zeros(1))
    
    #print('\n\nMino piena di zeri', mino[:], '\n\ndistr piena di zeri', distr[:])
    


    for h in range(0, len(x_histo), 2):
        
        x_histo[h] = mino[h]
        x_histo[h+1] = mino[h]
        
        y_histo[h] = distr[h]
        y_histo[h+1] = distr[h]
        
    
    
    y_histo.pop((len(y_histo)-1))
    y_histo.pop((len(y_histo)-1))
    y_histo.insert(0,0)
    y_histo.append(0)
    
    
    plt.plot(x_histo[:], y_histo[:], label='T=%i' %s)
    




plt.title('Distribuzione delle particelle')
plt.xlabel('Posizione')
plt.ylabel('rho')
plt.grid(True)
plt.legend(loc='best')
plt.show()










''' Creo gli istogrammi degli impulsi a mano '''



#Creiamo la rho teorica degli impulsi

avg_p    = np.zeros(N)
rho_p_teor = np.zeros(N)



f = np.linspace(-5, 5, N)

for h in range(1, N):
    
    avg_p[h] = 0
        
    rho_p_teor[h] = np.exp(-(f[h] - avg_p[h])*(f[h] - 
              avg_p[h])/(2*KT))/np.sqrt(2*np.pi*KT)
    
    
#creiamo l'istogramma numerico

for s in range(10, 150, 20):
    
    
    num_p = [0]*(bi+1)
    
    pmax=np.max(p[:,s])
    pmin=np.min(p[:,s])
    
    deltap= (pmax-pmin)/bi
    gino = np.linspace(pmin, pmax, (bi+1))
   
    for z in range(0, num_sim):
        
        a = (p[z,s]-pmin)/deltap
        
        num_p[int(a)] += 1 
        

    #plt.plot(gino[:], num_p[:]/(num_sim*deltap), label='T=%i' %s)   
    
    ##vale la stessa identica cosa scritta per le posizioni al paragrafo sopra
    
    distr_imp = num_p[:]/(num_sim*deltap)
    
    x_p_histo = [0]*2*len(gino)
    y_p_histo = [0]*2*len(distr_imp)



    for l in range(1, (2*len(gino)), 2):
        
        gino = np.insert(gino, l, np.zeros(1))
        distr_imp = np.insert(distr_imp, l, np.zeros(1))
        
        
        
    for o in range(0, len(x_p_histo), 2):
        
        x_p_histo[o] = gino[o]
        x_p_histo[o+1] = gino[o]
        
        y_p_histo[o] = distr_imp[o]
        y_p_histo[o+1] = distr_imp[o]
        
        
        
    y_p_histo.pop((len(y_p_histo)-1))
    y_p_histo.pop((len(y_p_histo)-1))
    y_p_histo.insert(0,0)
    y_p_histo.append(0)

    plt.plot(x_p_histo[:], y_p_histo[:], label='T=%i' %s)





plt.plot(f, rho_p_teor[:], color='black', label='Teorica')
plt.title('Distribuzione degli Impulsi')
plt.ylabel('rho')
plt.legend(loc='best')
plt.grid(True)
plt.show()








''' Andamenti delle popolazioni nelle buche di sinistra e destra '''

sinistra = np.zeros(num_sim)
destra   = np.zeros(num_sim)


somma_var_sx = np.zeros(N)
somma_var_dx = np.zeros(N)

p1_teor_var = np.zeros(N)
p2_teor_var = np.zeros(N)



##popolazioni ai tempi finali    
for h in range(0 , num_sim):
    
    sinistra[h] = frac_sx[h,N-1]
    destra[h]   = frac_dx[h,N-1]

Sx = np.sum(sinistra)
Dx = np.sum(destra)

print('\n\na sinistra ci sono:',Sx,'particelle\na destra ci sono:',
      Dx,'particelle\n\n')




##popolazioni al variare del tempo
for o in range(0, N):
    
    somma_var_sx[o] = np.sum(frac_sx[:,o])/num_sim
    somma_var_dx[o] = np.sum(frac_dx[:,o])/num_sim
    p1_teor_var[o]  = np.exp(-t[o]*(r1+r2)) + r2*(1-np.exp(-t[o]*(r1+r2)))/(r1+r2)
    p2_teor_var[o]  = r1*(1-np.exp(-t[o]*(r1+r2)))/(r1+r2)
    
    



distribuzioni = [Sx, Dx]

Prob_Sper_Sx = Sx/(num_sim)
Prob_Sper_Dx = Dx/(num_sim)



p1_sper = Sx/num_sim
p2_sper = Dx/num_sim


centri = [y_min1, y_min2]





##Andamenti delle probabilità al variare del tempo histogramma'''

for S in range(10, N, 150):
    
    plt.bar(centri, height=[np.sum(frac_sx[:,S-1]),np.sum(frac_dx[:,S-1])], 
                            width=0.85, alpha=0.1, 
                            label='%i T' %S )

plt.xlabel('Posizione')
plt.ylabel('Numero di particelle')
plt.title('Particelle a sinistra e destra al variare del tempo')
plt.legend(loc='best')    
plt.grid(True)    
plt.show()

    

##Anda,enti delle probabilità a fine simulazione '''

plt.bar(centri, height = distribuzioni, width=0.35 ,color='#0504aa', alpha=0.7)
plt.xlabel('Posizione')
plt.ylabel('Numero di particelle')
plt.title('Particelle a sinistra e destra a fine simulazione')
plt.grid(True)
plt.text(-0.25,160, 'p1_teor=%1.3f' %p1_teor, color = 'r')
plt.text(-0.25,130, 'p2_teor=%1.3f' %p2_teor, color = 'r')
plt.text(-0.25,90, 'p1_sper=%1.3f' %p1_sper, color = 'b')
plt.text(-0.25,50, 'p2_sper=%1.3f' %p2_sper, color = 'b')
plt.show()






##Andamenti delle probabilità al variare del tempo '''

'''for z in range(0, N, 100):
    
    plt.bar(t[z], height=somma_var_sx[z], width=4, color='b', 
            alpha=0.35)
    plt.bar(t[z], height=p1_teor_var[z], width=4, color='r', 
            alpha=0.35)

plt.xlabel('Tempo')
plt.title('Probabilità nella buca di Sinistra')
plt.grid(True)
#plt.legend(loc='best')    
plt.show()




for z in range(0, N, 100):
    
    plt.bar(t[z], height=somma_var_dx[z], width=4, color='b', 
            alpha=0.35)
    plt.bar(t[z], height=p2_teor_var[z], width=4, color='r', 
            alpha=0.35)

plt.xlabel('Tempo')
plt.title('Probabilità nella buca di Destra')
plt.grid(True)
#plt.legend(loc='best')    
plt.show()'''


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






























