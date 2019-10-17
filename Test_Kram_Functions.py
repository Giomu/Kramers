#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 11:56:57 2019

@author: giorgio
"""

''' Tests of Kramers.py '''


from hypothesis import given
import hypothesis.strategies as st
from Kram_Functions import Boltz, V, Calc_min_V, Prob_fuga, Prob_fin
import numpy as np


def test_Boltz_1():
        
    assert Boltz(1, 1) == 0.5
        




def test_V_1():
    
    assert V(1, 0, 0) == 1
    

def test_V_2():    
    
    assert V(0, 1, 1) == 0
    



#verifica qualche proprietà della derivata seconda
#ad esempio potresti verificare che integrandola si ottenga la derivata prima di V
#oppure che integrandola due volte ottieni proprio V
#magari anzichè farlo simpolicamente potresti farlo proprio numerico
#per qualche punto o magari fallo ocn hypothesis, non saprei.
    

def test_Calc_min_V_1():
    
    assert Calc_min_V(np.sqrt(2), np.sqrt(2)) == (-1.0, 1.0)





def test_Prob_fuga_1():
    
    assert Prob_fuga(-1.0, 1.0, 0.8000000000000002, np.sqrt(2), np.sqrt(2)) == (0.12897247163525316, 0.12897247163525316)    
    
    


@given(value = st.integers(200, 2500))
def test_Prob_fin_1(value):
    
    print(value)
    assert Prob_fin(value, 0.12897247163525316, 0.12897247163525316 ) == (0.5, 0.5)




















