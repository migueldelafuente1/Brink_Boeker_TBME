#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 17:57:34 2018

  * Generator of matrix elements in Antoine format

@author: miguel
"""
from __future__ import division, print_function

import Input_Output as io
import BM_brackets as bm
import Matrix_Elements as me

if __name__ == "__main__":
    
    bm.fact_build(200)
#    
#    n1 = 0
#    l1 = 0
#    j1 = 1
#    n2 = 0
#    l2 = 0
#    j2 = 1
#    J=1
#    T=0
#    
#    bra = me.QN_2body_jj_Coupling(n1,l1,j1, n2,l2,j2, J,T)
#    
#    ket = bra
#    parameters = [[0. for i in range(7)],[0. for i in range(7)]]
#    
#    parameters[0][0], parameters[1][0] = 1.3204, 1.3204
#    parameters[0][1], parameters[1][1] = 0.7,	1.4
#    parameters[0][2], parameters[1][2] = 714.085,	-72.612
#    parameters[0][3], parameters[1][3] = -443.085, -44.79
##    
##    bra.T = 1
##    ket.T = 1
##    ket.J = 0
##    ket.J = 0
#    print('\t BBME =',me.Brink_Boeker_ME(bra, ket, parameters))
#    
#    bra.T = 0
#    ket.T = 0
#    ket.J = 1
#    ket.J = 1
#    print('\t BBME =',me.Brink_Boeker_ME(bra, ket, parameters))
#    
#    bra.T = 1
#    ket.T = 1
#    ket.J = 1
#    ket.J = 1
#    print('\t BBME =',me.Brink_Boeker_ME(bra, ket, parameters))

    
    
    io.run_antoine_output('INPUT_P.txt')
    print(" The program has ended without any incidences.")