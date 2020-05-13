#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 17:56:08 2018

@author: miguel
"""
from __future__ import division, print_function
import numpy as np
import BM_brackets as bm
#from BM_brackets import fact
from sympy.physics.wigner import wigner_9j
from sympy import evalf

class QN_2body_jj_Coupling:
    
    ## Quantum Numbers for a 2 body Wave Function in jj coupling 
    ## with Couplign to total angular momentum J and total isospin T 
    
    # j1 and j2 are fractions and they are defined as 2*j
    def __init__(self, n1,l1,j1, n2,l2,j2, J,T):
         self.n1 = n1
         self.l1 = l1
         self.j1 = j1   
         self.n2 = n2
         self.l2 = l2
         self.j2 = j2
         self.J  = J
         self.T  = T
         
    def norm(self):
        # Norm of a 2 body antisimetriced wave function
        delta = 0
        
        if((self.n1 == self.n2) and (self.l1 == self.l2) 
            and(self.j1 == self.j2)):
            delta = 1
        
        return np.sqrt(1-delta*((-1)**(self.T+self.J))) / (1+delta)


def B_coeff(n,l, n_q,l_q, p, b_param):
    
    ## SHO normalization coefficients for WF
    
    # parity condition
    if(((l+l_q)%2)!=0):
        return 0
    
    Const = 0.5*((bm.fact[n]+bm.fact[n_q] + 
                   bm.fact[2*(n+l+1)]+bm.fact[2*(n_q+l_q+1)]) -
                    (bm.fact[n+l+1] + bm.fact[n_q+l_q+1]))
    Const += (bm.fact[2*(p+1)] - bm.fact[p+1])
    
    Const = ( ((-1)**(p-(l+l_q)/2)) * np.exp(Const) /
             ((b_param**3) * (2**(n+n_q))) )
    
    aux_sum = 0.
    major = min(n, (p-int((l+l_q)/2)) )
    minor = max(0, (p-int((l+l_q)/2)-n_q) )
    
    for k in range(minor,major+1):
        aux_sum += np.exp( (bm.fact[l+k+1] + bm.fact[p-k+1-int((l-l_q)/2)]) -
                    (bm.fact[k]+bm.fact[n-k]+bm.fact[2*(l+k+1)]+
                     bm.fact[p-int((l+l_q)/2)-k]+
                     bm.fact[n_q-p+k+int((l+l_q)/2)]+bm.fact[2*(p-k+1)+l_q-l]))
                     
                    
    
    return Const * aux_sum


def Talmi_integral(p, b_param, mu_param):
    ## Talmi integral for Gaussian potential shape 
    return (b_param**3) / ((1 + ((b_param/mu_param)**2) )**(p + 1.5))


def reduced_ME(n,l,n_q,l_q, b_param,mu_param):
    
    ## evaluate <nl|| V(r)||n_q,l_q> in the relative coordinates
    
    if(((l+l_q)%2)!=0):
        return 0
    
    aux_sum = 0.
    for p in range(int((l+l_q)/2), int((l+l_q)/2)+n+n_q +1):
        aux_sum += (B_coeff(n,l,n_q,l_q, p, b_param) * 
                    Talmi_integral(p, b_param, mu_param))
    
    return aux_sum


## Conditions of angular momentum and energy for the Matrix Elements
def angular_condition(l,L,lam):
    if abs(l-L)<=lam:
        if lam<=l+L:
            return True
    else:
        return False
def condicion(n,l,N,L, comparar):
    if(2*n + l + 2*N + L == comparar):
        return True
    else:
        return False
    
    
def Moshinsky_Transformation_ME(n_a, l_a, n_b, l_b, 
                                n_c, l_c, n_d, l_d, lamda, b_param, mu_param):
    
    ## This function applies the transformation of the lambda coupled 
    ## 2 body wave functions from the coordinate system to the relative + 
    ## mass center coordinates using the  Moshinsky transformation
    
    # Are l_a,l_b,l_c and l_d fullfiling the angular momentum conservation ?
    if ((not angular_condition(l_a,l_b,lamda)) or 
        (not angular_condition(l_c,l_d,lamda))):
        return 0
    # BM Bracket parity restriction
    if((l_c+l_d-l_a-l_b)%2 != 0):
        return 0
    
    Sum = 0.
        
    rho = 2*(n_a + n_b) + l_a + l_b
    # N,L and n are restricted to non-negative and Energy-resticted values,
    # althought, these variables can minimice the iterations following this order:
    for N in range(np.floor_divide(rho,2)  +1):
        for n in range(np.floor_divide(rho,2) - N  +1):            
            for L in range(0,rho - 2*(n + N)  +1):
                
                # l is N,L,n dependent
                l = rho - L - 2*(N+n)
                
                if (l >= 0):
                    
                    if(angular_condition(l,L,lamda)):
                        
                        if condicion(n,l,N,L, rho):
                            
                            # Perform the BM Brackets and the reduced Matrix 
                            # Element n' is energetically dependent of rho 
                            # and n l',L' and N' must be l,L & N
                            
                            BMB_bra = bm.BM_Bracket(n,l,N,L, 
                                                    n_a,l_a,n_b,l_b, lamda)
                            if(abs(BMB_bra) > 1.e-15 ):
                                
                                n_q = ( (n_c+n_d) - (n_a+n_b) + 
                                       int((l_c+l_d-l_a-l_d)/2) + n)
                                BMB_ket = bm.BM_Bracket(n_q,l,N,L, 
                                                        n_c,l_c,n_d,l_d, lamda)
                                if(abs(BMB_ket) > 1.e-15 ):
                                    ## BM Brackets are non-cero, perform <nl|| V ||n'l'>
                                    m_element = reduced_ME(n,l,n_q,l,
                                                           b_param,mu_param)
                                    
#                                    print("BMB1=",BMB_bra," BMB2=",BMB_ket,
#                                          "  <V(r)>=",m_element)
                                    Sum += m_element * BMB_bra * BMB_ket

                            
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
    
    return Sum
    
    
def non_antisimetrized_BB_ME(bra, ket, parameters):
    
    # Passing of parameters
    n_a, n_b = bra.n1, bra.n2
    n_c, n_d = ket.n1, ket.n2
    l_a, l_b = bra.l1, bra.l2
    l_c, l_d = ket.l1, ket.l2
    j_a, j_b = bra.j1, bra.j2
    j_c, j_d = ket.j1, ket.j2
    
    J, T = bra.J, bra.T
    
    # Adjust the oscillator lenght to the reducced mass of the proton
    b_param = parameters[0][0] * np.sqrt(2)
    
    Sum = 0.
    L_max = min((l_a+l_b),(l_c+l_d))
    L_min = max(abs(l_a-l_b),abs(l_c-l_d))
    
    for L in range(L_min, L_max+1):
        for S in range(0,2): 
            
            try: 
                # j atribute are defined as 2*j
               w9j_bra = wigner_9j(l_a,1/2,j_a/2, l_b,1/2,j_b/2,
                                   L,S,J, prec=None)#.evalf()
            except ValueError:
                w9j_bra = 0
            
            if(abs(w9j_bra) > 1.e-15 ):
                try:
                    w9j_ket = wigner_9j(l_c,1/2,j_c/2, l_d,1/2,j_d/2, 
                                        L,S,J, prec=None)#.evalf()
                except ValueError:
                    w9j_ket = 0
                
                if(abs(w9j_ket) > 1.e-15 ):
                        # j atribute are defined as 2*j
                    recoupling = ((2*S+1) * (2*L+1) *
                                  np.sqrt((j_a+1)*(j_b+1)*(j_c+1)*(j_d+1)) *
                                   w9j_bra * w9j_ket)
                    
                    # Sum of gaussians and projection operators
                    for i in range(2):
                        
                        mu_param = parameters[i][1]
                        
                        # Radial Part for Gauss Integral (L == lambda)
                        Radial = Moshinsky_Transformation_ME(n_a, l_a, n_b, l_b, 
                                n_c, l_c, n_d, l_d, L, b_param, mu_param)
                        
#                        print('Radial['+ str(i)+ ']: '+str(Radial))
                        # Exchange Part
                        Exchange_Energy = ( parameters[i][2] +                  # Wigner
                                           (parameters[i][3]*((-1)**(L))) +     # Majorana
                                           (parameters[i][4]*((-1)**(S))) +     # Barlett
                                           (parameters[i][5]*((-1)**(T))))      # Heisenberg
                        
#                        print('Exchange['+ str(i)+ ']: ',Exchange_Energy, '  couping=',recoupling)
                        # Add up
                        Sum += recoupling * Radial * Exchange_Energy
                
    
    
    return Sum


def Brink_Boeker_ME(bra, ket, parameters):
    
    # Passing of parameters
    n_a, n_b = bra.n1, bra.n2
    n_c, n_d = ket.n1, ket.n2
    l_a, l_b = bra.l1, bra.l2
    l_c, l_d = ket.l1, ket.l2
    j_a, j_b = bra.j1, bra.j2
    j_c, j_d = ket.j1, ket.j2
    
    J,J_q, T,T_q = bra.J, ket.J, bra.T, ket.T 
         
    if(J != J_q):
        return 0
    if(T != T_q):
        return 0

    #  Isospin and total coupled angular momentum 
    #  must be odd it the states are in the same orbit:
    if(((n_a == n_b) and (l_a == l_b) and (j_a == j_b)) or
       ((n_c == n_d) and (l_c == l_d) and (j_c == j_d))):
        if((J + T)%2 != 1):
            return 0

    # construct the exchange ket
    exchange_ket = QN_2body_jj_Coupling
    exchange_ket.n2 = ket.n1
    exchange_ket.n1 = ket.n2
    exchange_ket.l2 = ket.l1
    exchange_ket.l1 = ket.l2
    exchange_ket.j2 = ket.j1
    exchange_ket.j1 = ket.j2
    exchange_ket.J = ket.J
    exchange_ket.T = ket.T
    
    fase = (-1)**(((j_c+j_d)/2) - (J+T) + 1)   
    
    return ( bra.norm() * ket.norm() * (2*J + 1) *
            ( non_antisimetrized_BB_ME(bra, ket, parameters) - 
             (fase * non_antisimetrized_BB_ME(bra, exchange_ket, parameters))))
