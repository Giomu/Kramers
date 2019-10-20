#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 18:41:58 2019

@author: giorgio
"""

''' Functions of Kramers.py '''


import numpy as np








def Boltz(eps, gamma):
    
    ''' \n\nThis function computes the Boltzmann Constant given epsilon and gamma
    as parameters.\n\n-->epsilon: modulation of stocasthic noise given from params.txt
-->  gamma: friction of the system given from params.txt '''
        
    return eps*eps/(2*gamma)





def V(y, a, b):
    
    ''' \n\nThis function returns the Potential V(y) = (y**2 - ay)(y**2 + by)
    given a and b from imput.\n
    --> a: Alfa given from params.txt
    --> b: Beta given from params.txt'''
    
    return ((y*y-a*y)*(y*y+b*y))





def V_II(y, a, b):
    
    ''' This function returns the second derivative of the Potential 
    given a and b from imput. \n
    --> a: Alfa given from params.txt
    --> b: Beta given from params.txt'''
    
    return (12*y*y + 6*y*(b-a) - 2*a*b)





def Calc_min_V(a, b):   
    
    ''' This function returns the left and right minimums of the potential V(x)
    printing their positions on the x-axes. \n
    --> a: Alfa given from params.txt
    --> b: Beta given from params.txt'''
    
    y_minsx = (1/8)*(3*(a-b) - np.sqrt(9*(a-b)*(a-b)+32*a*b))  #min di sx
    y_mindx = (1/8)*(3*(a-b) + np.sqrt(9*(a-b)*(a-b)+32*a*b))  #min di dx
    print('\nLeft minimum:',y_minsx,'\nRight minimum:',y_mindx)
    
    return y_minsx, y_mindx





def Prob_fuga(y_min1, y_min2, KT, a, b):
    
    ''' this function returns the left and right probability of escape r1 and r2 
    from the two holes of the potential V(x). \n
    --> y_min1: left minimum of the potential V(x) returned from Calc_min_V(a,b)
    --> y_min2: right minimum of the potential V(x) returned from Calc_min_V(a,b)
    --> KT    : Boltzmann Constant returned from Boltz(eps, gamma)
    --> a: Alfa given from params.txt
    --> b: Beta given from params.txt'''
    
    P_f1 = (np.sqrt(V_II(y_min1, a, b))/(2*np.pi))*np.exp((-V(0, a, b)+V(y_min1, a, b))/KT)
    P_f2 = (np.sqrt(V_II(y_min2, a, b))/(2*np.pi))*np.exp((-V(0, a, b)+V(y_min2, a, b))/KT)
    
    return P_f1, P_f2





def Prob_fin(N, r1, r2):
    
    ''' This function returns the probability that at the end of the simulation
    the particle is in the left or right hole of the potential V(x). \n
    --> N: is the number of time step it is fixed for all the simulation
    --> r1: Probability of left escape returned from Prob_fuga(y_min1, y_min2, KT, a, b)
    --> r2: Probability of right escape returned from Prob_fuga(y_min1, y_min2, KT, a, b)'''
    
    p1_fin = np.exp(-N*(r1+r2)) + r2*(1-np.exp(-N*(r1+r2)))/(r1+r2)
    p2_fin = r1*(1-np.exp(-N*(r1+r2)))/(r1+r2)
    
    return p1_fin, p2_fin








































