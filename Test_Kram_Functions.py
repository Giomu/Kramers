#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 11:56:57 2019

@author: giorgio
"""

''' Tests of Kramers.py '''



from hypothesis import given
import hypothesis.strategies as st
from Kram_Functions import Boltz, V, V_II, Calc_min_V, Prob_fuga, Prob_fin
import numpy as np








def test_Boltz_1():
    
    ''' Tests that eps=1, gamma=1 returns Boltz(eps, gamma) == 0.5 '''
    
    assert Boltz(1, 1) == 0.5
        


def test_Boltz_2():
    
    ''' Tests that eps=0, gamma=0.1 returns Boltz(eps, gamma) == 0 '''
    
    assert Boltz(0, 0.1) == 0



@given(gamma=st.floats(min_value=0, allow_infinity=False, exclude_min=True))
def test_Boltz_3(gamma):
    
    ''' Tests that eps=0, gamma contained in (0, +inf) returns Boltz(eps, gamma) == 0 '''
    
    assert Boltz(0, gamma) == 0







def test_V_1():
    
    ''' Tests that y=1, a=b=0 returns V(y, a, b) == 1 '''
    
    assert V(1, 0, 0) == 1
    


def test_V_2(): 
    
    ''' Tests that y=0, a=b=1 returns V(y, a, b) == 0 '''
    
    assert V(0, 1, 1) == 0
    


@given(a=st.floats(allow_nan=False, allow_infinity=False), b=st.floats(allow_nan=False, allow_infinity=False))
def test_V_3(a, b):
    
    ''' Tests that for every a and b positive floats, y=0, returns V(y, a, b) == 0 '''
    
    assert V(0, a, b) == 0







@given(a1=st.floats(-50, 50), b1=st.floats(-50, 50))
def test_V_II_1(a1, b1):
    
    ''' Tests that for every a and b contained in [-50, 50], 
        V_II(y, a, b) returns -2*a*b as expected'''
    
    print('a1: ', a1, '    b1: ', b1)
    assert V_II(0, a1, b1) == -2*a1*b1
        


@given(y=st.floats(-5000, 5000))
def test_V_II_2(y):
    
    ''' Tests that for y contained in [-5000, 5000] and a=b=0,
        V_II(y, a, b) returns 12*y^2 as expected'''
    
    print('y: ', y)
    
    assert V_II(y, 0, 0) == 12*y*y




    


def test_Calc_min_V_1():
    
    ''' Tests that for a=b=sqrt(2), Calc_min_V(a, b) returns (-1, 1) as expected '''
    
    assert Calc_min_V(np.sqrt(2), np.sqrt(2)) == (-1.0, 1.0)




@given(a=st.floats(1, 150), b=st.floats(1, 150))
def test_Calc_min_V_2(a, b):
    
    ''' Tests that for a and b contained in [1, 150], Calc_min_V(a, b) returns
        negative left min and a positive right min as expected '''
    
    print(a, b)
    x = Calc_min_V(a, b)
    print(x)
    
    assert x[0] < 0 and x[1] > 0
    



@given(a=st.floats(1, 5000))
def test_Calc_min_V_3(a):
    
    ''' Test that given a and b if they are the same then Calc_min_V(a, b)
        returns two symmetrical min of the potential '''
    
    b = a
    print('a: ',a,'  b: ', b)
    x = Calc_min_V(a, b)
    if a == b:
        assert abs(x[0]) == abs(x[1])




@given(a2=st.floats(1, 500))
def test_Calc_min_V_4(a2):
    
    ''' Test that if a=b then Calc_min_V(a, b) returns two min in
        +\-sqrt(32*a*b)/8 as expected'''
    
    b2 = a2
    print(a2, b2)
    x = Calc_min_V(a2, b2)
    assert x[1] == np.sqrt(32*a2*b2)/8 and x[0] == -np.sqrt(32*a2*b2)/8








def test_Prob_fuga_1():
    
    ''' Tests that given y_min1=-1, y_min2=1, KT=0.8, a=b=sqrt(2),
        Prob_fuga(y_min1, y_min2, KT, a, b) returns the two probabilities 
        r1 and r2 as 0.12897247163525316 as expected'''
    
    assert Prob_fuga(-1.0, 1.0, 0.8000000000000002, np.sqrt(2), np.sqrt(2)) == (0.12897247163525316, 0.12897247163525316)    
     
 
    
def test_Prob_fuga_2():
    
    ''' Tests that for y_min1=y_min2=a=b=0 and KT=1,
        Prob_fuga(y_min1, y_min2, KT, a, b) returns the two probabilities 
        r1 and r2 as 0.0 as expected'''
    
    assert Prob_fuga(0, 0, 1, 0, 0) == (0, 0)



def test_Prob_fuga_3():
    
    '''Tests that for y_min1=y_min2=a=b=KT=1,
        Prob_fuga(y_min1, y_min2, KT, a, b) returns the two probabilities 
        r1=r2 as sqrt(10)/(2*pi) as expected'''
    
    assert Prob_fuga(1, 1, 1, 1, 1) == (np.sqrt(10)/(2*np.pi), np.sqrt(10)/(2*np.pi))    







@given(value = st.integers(200, 2500))
def test_Prob_fin_1(value):
    
    ''' Tests that for N contained in [200, 2500] and
        r1=r2=0.12897247163525316, Prob_fin(N, r1, r2) returns as
        p_fin1=p_fin2=0.5 as expected'''
    
    print(value)
    assert Prob_fin(value, 0.12897247163525316, 0.12897247163525316 ) == (0.5, 0.5)



@given(N=st.integers(200, 20000))
def test_Prob_fin_2(N):
    
    ''' Tests that for N contained in [200, 20000], r1=1, r2=0,
    Prob_fin(N, r1, r2) returns pfin_1=exp(-N) and pfin_2=(1-exp(-N)) as expected'''
    
    print(N)
    assert Prob_fin(N, 1, 0) == (np.exp(-N), (1-np.exp(-N)))


    
@given(N=st.integers(200, 20000))
def test_Prob_fin_3(N):
    
    ''' Tests that for N contained in [200, 20000], r1=0 and r2=1,
    Prob_fin(N, r1, r2) returns p_fin1=1 and p_fin2=0 as expected'''
    
    print(N)
    assert Prob_fin(N, 0, 1) == (1, 0)
















