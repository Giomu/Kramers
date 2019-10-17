#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 18:41:58 2019

@author: giorgio
"""

''' Functions of Kramers.py '''


import numpy as np








def Boltz(eps, gamma):
    
    ''' Questa funzione calcola la costante di Boltzmann
        a partire dai valori di epsilon e gamma forniti da input '''
        
    return eps*eps/(2*gamma)





def V(y, a, b):
    
    ''' Questa funzione crea il Potenziale V(x) a partire dai 
        valori di a e b forniti da input'''
    
    return ((y*y-a*y)*(y*y+b*y))





def V_II(y, a, b):
    
    ''' Questa funzione calcola la derivata seconda del Potenziale V(x) '''
    
    return (12*y*y + 6*y*(b-a) - 2*a*b)





def Calc_min_V(a, b):   
    
    ''' Questa funzione calcola i minimi a sinistra e destra
        del Potenziale V(x) printandone la posizione sull'asse x'''
    
    y_minsx = (1/8)*(3*(a-b) - np.sqrt(9*(a-b)*(a-b)+32*a*b))  #min di sx
    y_mindx = (1/8)*(3*(a-b) + np.sqrt(9*(a-b)*(a-b)+32*a*b))  #min di dx
    print('\nminimo di sinistra:',y_minsx,'\nminimo di destra:',y_mindx)
    
    return y_minsx, y_mindx





def Prob_fuga(y_min1, y_min2, KT, a, b):
    
    ''' Questa funzione calcola le probabilità di fuga
        dalle due buche di sinistra e di destra P_f1 e P_f2'''
    
    P_f1 = (np.sqrt(V_II(y_min1, a, b))/(2*np.pi))*np.exp((-V(0, a, b)+V(y_min1, a, b))/KT)
    P_f2 = (np.sqrt(V_II(y_min2, a, b))/(2*np.pi))*np.exp((-V(0, a, b)+V(y_min2, a, b))/KT)
    
    return P_f1, P_f2





def Prob_fin(N, r1, r2):
    
    ''' Questa funzione calcola la probabilità che a fine simulazione,
        ovvero al tempo finale, la particella si trovi nella buca di
        sinistra o di destra rispettivamente'''
    
    p1_fin = np.exp(-N*(r1+r2)) + r2*(1-np.exp(-N*(r1+r2)))/(r1+r2)
    p2_fin = r1*(1-np.exp(-N*(r1+r2)))/(r1+r2)
    
    return p1_fin, p2_fin








































