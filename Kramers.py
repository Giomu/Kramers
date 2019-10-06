#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 11:55:30 2019

@author: giorgio
"""

''' Kramers Project Software and Computing '''





import random
import numpy as np
import logging
from scipy.misc import derivative
import matplotlib.pyplot as plt
#import seaborn as sns





''' Metti un bel messaggio di presentazione della simulazione spiegando
    davvero tutto quello che fa il programma '''



'''Definisco i parametri della simulazione'''
'''
N       = 800                 #Numero di passi dell'integrazione
dt      = 0.1                 #Incrementino temporale

eps     = 0.4                 #modula l'effetto stocastico
gamma   = 0.1                 #Attrito
KT      = eps*eps/(2*gamma)   #kT (epsilon=1)

num_sim = 1000                #numero delle simulazioni da svolgere
'''                               #ovvero il numero di particelle


N       = 800    #Numero passi di integrazione
dt      = 0.1    


while True:
    try:
        eps = float(input("PLEASE ENTER EPSILON: "))
    
    except ValueError:
        logging.error("Sorry, number should be a float digit.")
        continue
    
    if eps < 0.3 or eps>1:
        logging.warning('''Sorry but epsilon should be in range (0.3, 1) to see 
                            some effects.. Try again!''')
        continue
    
    else:
        break

print('eps: ', eps)



while True:
    try:
        gamma = float(input("PLEASE ENTER GAMMA: "))
    
    except ValueError:
        logging.error("Sorry, number should be a float digit.")
        continue
    
    if gamma < 0 or gamma>1:
        logging.warning('''Sorry but gamma should be in range (0, 1).. Try again!''')
        continue
    
    else:
        break

print('gamma: ', gamma)



while True:
    try:
        num_sim = int(input("PLEASE ENTER NUMBER OF PARTICLES N: "))
    
    except ValueError:
        logging.error("Sorry, number should be an int digit.")
        continue
    
    if num_sim < 1000:
        logging.warning('''Sorry but the error goes as 1/sqrt(N) 
                            give me a number bigger then 1000''')
        continue
    
    else:
        break

print('N: ', num_sim)






def Boltz(eps, gamma):
    return eps*eps/(2*gamma)

KT = Boltz(eps, gamma)
#'''





'''Calcolo il Rumore bianco con distribuzione gaussiana a 
media nulla e varianza 1'''

mu = 0
sigma = 1

csi = np.zeros((num_sim, N))

for t in range (0, num_sim):
    for k in range(0, N):
        csi[t,k] = random.gauss(mu, sigma)
    
#verifica quella cosa del vecchio progr




'''Definisco il Potenziale'''

#richiedili da imput e si potrebbe verificare il tipo cioè che effettivamente siano numeri e non altro
a = 2#np.sqrt(2)    #Parametri del potenziale asimmetrico per individuare
b = 1#np.sqrt(2)    #i punti dei due minimi delle due buche



def V(y):
    
    ''' Questa funzione crea il Potenziale V(x) a partire dai 
        valori di a e b forniti da input'''
    
    return ((y*y-a*y)*(y*y+b*y))



def V_II(y):
    
    ''' Questa funzione calcola la derivata seconda del Potenziale V(x) '''
    
    return (12*y*y + 6*y*(b-a) - 2*a*b)


#********QUI HA SENSO CREARE UNA CLASS KRAM CHE CALCOLI y_min1,2 r1,2 pi_teor????????????????*********


##mi conviene mettere tutte queste funzioni in un altro file e chiamare solo le funzioni???


#y_max = 0                                                 #max di V(y)

def Calc_min_V(a, b):   
    
    ''' Questa funzione calcola i minimi a sinistra e destra
        del Potenziale V(x) printandone la posizione sull'asse x'''
    
    y_minsx = (1/8)*(3*(a-b) - np.sqrt(9*(a-b)*(a-b)+32*a*b))  #min di sx
    y_mindx = (1/8)*(3*(a-b) + np.sqrt(9*(a-b)*(a-b)+32*a*b))  #min di dx
    print('\nminimo di sinistra:',y_minsx,'\nminimo di destra:',y_mindx)
    
    return y_minsx, y_mindx


y_min1, y_min2 = Calc_min_V(a,b)
    


#potremmo printare il grafico del potenziale




def Prob_fuga(y_min1, y_min2):
    
    ''' Questa funzione calcola le probabilità di fuga
        dalle due buche di sinistra e di destra P_f1 e P_f2'''
    
    P_f1 = (np.sqrt(V_II(y_min1))/(2*np.pi))*np.exp((-V(0)+V(y_min1))/KT)
    P_f2 = (np.sqrt(V_II(y_min2))/(2*np.pi))*np.exp((-V(0)+V(y_min2))/KT)
    
    return P_f1, P_f2

r1, r2 = Prob_fuga(y_min1, y_min2)


#Si puo mettere in una funzione?
if a == b:
    assert r1 == r2

if a > b:
    assert r1 > r2



def Prob_fin(N, r1, r2):
    
    ''' Questa funzione calcola la probabilità che a fine simulazione,
        ovvero al tempo finale, la particella si trovi nella buca di
        sinistra o di destra rispettivamente'''
    
    p1_fin = np.exp(-N*(r1+r2)) + r2*(1-np.exp(-N*(r1+r2)))/(r1+r2)
    p2_fin = r1*(1-np.exp(-N*(r1+r2)))/(r1+r2)
    
    return p1_fin, p2_fin

p1_teor, p2_teor = Prob_fin(N, r1, r2)



#non so bene come farlo girare in una funzione
'''Verify that the sum of the 2 probabilities are equal to one'''
assert p1_teor, p2_teor >= 0
assert p1_teor + p2_teor == 1



#************     fino a qui    **************



''' Integro il moto per le diverse simulazioni '''

x0 = y_min1          #posizione iniziale nel minimo di sinistra
p0 = 0               #impulso iniziale fermo

    
x = np.zeros((num_sim, N))   #matrici contenenti le particelle nelle righe
p = np.zeros((num_sim, N))   #e gli istanti temporali nelle colonne
U = np.zeros((num_sim, N))

frac_sx = np.zeros((num_sim, N))  #fraction of particles in the left side
frac_dx = np.zeros((num_sim, N))  #fraction of particles in the right side

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
        
    
    distr = num_x[:]/(num_sim*deltax)
    
    plt.step(mino[:], distr[:], label='T=%i' %s)


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
        

    plt.step(gino[:], num_p[:]/(num_sim*deltap), label='T=%i' %s)   



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






























