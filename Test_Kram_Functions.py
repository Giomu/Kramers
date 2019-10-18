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
        
    assert Boltz(1, 1) == 0.5
        

def test_Boltz_2():
    
    assert Boltz(0, 0.1) == 0


@given(gamma=st.floats(min_value=0, allow_infinity=False, exclude_min=True))
def test_Boltz_3(gamma):
    
    assert Boltz(0, gamma) == 0






def test_V_1():
    
    assert V(1, 0, 0) == 1
    

def test_V_2():    
    
    assert V(0, 1, 1) == 0
    

@given(a=st.floats(allow_nan=False, allow_infinity=False), b=st.floats(allow_nan=False, allow_infinity=False))
def test_V_3(a, b):
    
    assert V(0, a, b) == 0






@given(a1=st.floats(-50, 50), b1=st.floats(-50, 50))
def test_V_II_1(a1, b1):
    print('a1: ', a1, '    b1: ', b1)
    assert V_II(0, a1, b1) == -2*a1*b1
        

@given(y=st.floats(-5000, 5000))
def test_V_II_2(y):
    print('y: ', y)
    
    assert V_II(y, 0, 0) == 12*y*y




    

def test_Calc_min_V_1():
    
    assert Calc_min_V(np.sqrt(2), np.sqrt(2)) == (-1.0, 1.0)


#test left min is negative, while right min is positive
@given(a=st.floats(1, 150), b=st.floats(1, 150))
def test_Calc_min_V_2(a, b):
    
    print(a, b)
    x = Calc_min_V(a, b)
    print(x)
    
    assert x[0] < 0 and x[1] > 0
    

#Test that given a and b if they are the same then the 2 min of the potential are simmetrical
@given(a=st.floats(1, 5000))
def test_Calc_min_V_3(a):
    b = a
    print('a: ',a,'  b: ', b)
    x = Calc_min_V(a, b)
    if a == b:
        assert abs(x[0]) == abs(x[1])


#test that if a=b then all works fine
@given(a2=st.floats(1, 500))
def test_Calc_min_V_4(a2):
    b2 = a2
    print(a2, b2)
    x = Calc_min_V(a2, b2)
    assert x[1] == np.sqrt(32*a2*b2)/8 and x[0] == -np.sqrt(32*a2*b2)/8







def test_Prob_fuga_1():
    
    assert Prob_fuga(-1.0, 1.0, 0.8000000000000002, np.sqrt(2), np.sqrt(2)) == (0.12897247163525316, 0.12897247163525316)    
     
    
def test_Prob_fuga_2():
    
    assert Prob_fuga(0, 0, 1, 0, 0) == (0, 0)


def test_Prob_fuga_3():
    
    assert Prob_fuga(1, 1, 1, 1, 1) == (np.sqrt(10)/(2*np.pi), np.sqrt(10)/(2*np.pi))    






@given(value = st.integers(200, 2500))
def test_Prob_fin_1(value):
    
    print(value)
    assert Prob_fin(value, 0.12897247163525316, 0.12897247163525316 ) == (0.5, 0.5)


@given(N=st.integers(200, 20000))
def test_Prob_fin_2(N):
    
    print(N)
    assert Prob_fin(N, 1, 0) == (np.exp(-N), (1-np.exp(-N)))

    
@given(N=st.integers(200, 20000))
def test_Prob_fin_3(N):
    
    print(N)
    assert Prob_fin(N, 0, 1) == (1, 0)
















