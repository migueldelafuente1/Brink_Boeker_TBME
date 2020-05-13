#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 17:55:30 2018

@author: miguel
"""
from __future__ import division, print_function
import numpy as np

from sympy.physics.wigner import racah
from sympy import evalf

def fact_build(max_order_fact):
    # construct a global factorial base
    global fact
    
    fact = [np.log(1)]
    for i in range(1,max_order_fact):
        fact.append(fact[i-1] + np.log(i))
    
    # !! logarithm values for the factorials !!



def matrix_r2(n,l,N,L, n_q,l_q,N_q,L_q, lamda):
    if((n_q<0) or (l_q<0) or (N_q<0) or(L_q<0)):
        return 0
    
    if(n_q == (n-1)):
        # relations 1,3,4
        if(l_q == l):
            if(N_q == N):
                if(L_q == L):
                    # rel 1
                    return 0.5*np.sqrt(n*(n+l+0.5))
                else:
                    return 0
            else:
                return 0
        elif(l_q == (l+1)):
            if(N_q == (N-1)):
                if(L_q == (L+1)):
                    # rel 3
                    try:
                        racah_aux = (racah(l,(l+1),L,(L+1),1,lamda)).evalf()
                    except AttributeError:
                        racah_aux = 0
                        
                    return (((-1)**(lamda+L+l))*
                            np.sqrt(n*N*(l+1)*(L+1)) * racah_aux)
                            
                else:
                    return 0
            elif(N_q == N):
                if(L_q == (L-1)):
                    # rel 4
                    try:
                        racah_aux = (racah(l,(l+1),L,(L-1),1,lamda)).evalf()
                    except AttributeError:
                        racah_aux = 0
                        
                    return (((-1)**(lamda+L+l))*
                            np.sqrt(n*L*(l+1)*(N+L+0.5)) * racah_aux)
                            
                else:
                    return 0
            else: 
                return 0
            
        else:
            return 0
            
        
    elif(n_q == n):
        # relations 2,5,6
        if(l_q == (l-1)):
            if(N_q == (N-1)):
                if(L_q == (L+1)):
                    # rel 5
                    try:
                        racah_aux = (racah(l,(l-1),L,(L+1),1,lamda)).evalf()
                    except AttributeError:
                        racah_aux = 0
                        
                    return (((-1)**(lamda+L+l))*
                            np.sqrt(N*l*(n+l+0.5)*(L+1))*racah_aux)
                            
                else:
                    return 0
            elif(N_q == N):
                if(L_q == (L-1)):
                    # rel 6
                    try:
                        racah_aux = (racah(l,(l-1),L,(L-1),1,lamda)).evalf()
                    except AttributeError:
                        racah_aux = 0
                        
                    return (((-1)**(lamda+L+l))*
                            np.sqrt(L*l*(n+l+0.5)*(N+L+0.5))*racah_aux)
                            
                else:
                    return 0
            else:
                return 0
                
        elif(l_q == l):
            if(N_q == (N-1)):
                if(L_q == L):
                    # rel 2
                    return 0.5*np.sqrt(N*(N+L+0.5))
                else:
                    return 0
            else:
                return 0
        else:
            return 0
    
    else:
        return 0
    
    return ValueError


def A_coeff(l1,l,l2,L, x):
    # Coefficient for the BMB_00
    const = (0.5 * (fact[l1+l+x+1]+fact[l1+l-x]+fact[l1+x-l]-
                    fact[l+x-l1]))
    
    const += (0.5 * (fact[l2+L+x+1]+fact[l2+L-x]+fact[l2+x-L]-
                    fact[L+x-l2]))
    
    const = np.exp(const)
    aux_sum = 0.
    
    # limits for non negative factorials
    # q is non negative
    c1 = l1 - l
    c2 = l2 - L
    c3 = -x - 1

    c4 = l + l1
    c5 = L + l2

    major = min(c4,c5);
    minor = max(max(max(max(c1,c2),c3),x),0)
    
    for q in range(minor, major+1):
        
        if( ((l+q-l1)%2) == 0 ):
            numerator = (fact[l+q-l1] + fact[L+q-l2])
            
            denominator = ((fact[int((l+q-l1)/2)] 
                + fact[int((l+l1-q)/2)] + fact[q-x] + fact[q+x+1] +
                fact[int((L+q-l2)/2)] + fact[int((L+l2-q)/2)]))
            
            aux_sum += (((-1)**((l+q-l1)/2)) * 
                np.exp(numerator - denominator))
    
    return const * aux_sum


def BM_Bracket00(n,l,N,L, l1,l2, lamda):
    # Limit of the recurrence relation   

    const = ((fact[l1] + fact[l2] + fact[n+l] + fact[N+L]) -
            (fact[2*l1] + fact[2*l2] + fact[n] + fact[N] +
             fact[2*(n+l)+1] + fact[2*(N+L)+1]))
    const += ((np.log(2*l+1)+np.log(2*L+1)) - ((l+L)*np.log(2)))
    
    aux_sum = 0.
    
    major = min((l+l1),(L+l2))
    minor = max(abs(l-l1),abs(L-l2))
    for x in range(minor, major+1):
        try:
            racah_aux = racah(l,L,l1,l2, lamda,x).evalf()
        except AttributeError:
            racah_aux = 0
        
        aux_sum += ((2*x+1)*A_coeff(l1,l,l2,L, x)*racah_aux)
            
    return np.exp(0.5*const) * aux_sum * ((-1)**(n+l+L-lamda))


def BM_Bracket(n,l,N,L, n1,l1,n2,l2, lamda):
    
    # Non-negative conditions over constants
    if((n<0) or (l<0) or (N<0) or(L<0)):
        return 0
    if((n1<0) or (l1<0) or (n2<0) or (l1<0) or (lamda<0)):
        return 0
    
    # Energy condition
    if((2*(n1+n2)+l1+l2) != (2*(n+N)+l+L)):
        return 0
    # Angular momentum conservation
    if((abs(l1-l2) > lamda) or ((l1+l2) < lamda)):
        return 0
    if((abs(l-L) > lamda) or ((l+L) < lamda)):
        return 0    
    
    # RECURRENCE RELATIONS
    # there is only 6 non-zero combinations of n'l'N'L' 
    if(n1 == 0):
        if(n2 == 0):
            # BMB00
            return BM_Bracket00(n,l,N,L, l1,l2, lamda)
        else:
            # permutate the n1 l1 with n2 l2
            fase = ((-1)**(L-lamda))
            aux_sum = 0.
            aux_sum += fase * (matrix_r2(n,l,N,L ,n-1,l,N,L, lamda)*
                    BM_Bracket(n-1,l,N,L, n2-1,l2,n1,l1, lamda))
            aux_sum += fase * (matrix_r2(n,l,N,L ,n,l,N-1,L, lamda)*
                    BM_Bracket(n,l,N-1,L, n2-1,l2,n1,l1, lamda))
            aux_sum += fase * (matrix_r2(n,l,N,L, n-1,l+1,N-1,L+1, lamda)*
                    BM_Bracket(n-1,l+1,N-1,L+1,  n2-1,l2,n1,l1, lamda))
            aux_sum += fase * (matrix_r2(n,l,N,L, n-1,l+1,N,L-1, lamda)*
                    BM_Bracket(n-1,l+1,N,L-1, n2-1,l2,n1,l1, lamda))
            aux_sum += fase * (matrix_r2(n,l,N,L, n,l-1,N-1,L+1, lamda)*
                    BM_Bracket(n,l-1,N-1,L+1, n2-1,l2,n1,l1, lamda))
            aux_sum += fase * (matrix_r2(n,l,N,L, n,l-1,N,L-1, lamda)*
                    BM_Bracket(n,l-1,N,L-1, n2-1,l2,n1,l1, lamda))

            
            return  np.sqrt(1./(n2*(n2+l2+0.5))) * aux_sum
            
    else:
        # normal case
        aux_sum = 0.
        aux_sum += (matrix_r2(n,l,N,L, n-1,l,N,L, lamda)*
                    BM_Bracket(n-1,l,N,L, n1-1,l1,n2,l2, lamda))
        aux_sum += (matrix_r2(n,l,N,L, n,l,N-1,L, lamda)*
                    BM_Bracket(n,l,N-1,L, n1-1,l1,n2,l2, lamda))
        aux_sum += (matrix_r2(n,l,N,L, n-1,l+1,N-1,L+1, lamda)*
                    BM_Bracket(n-1,l+1,N-1,L+1, n1-1,l1,n2,l2, lamda))
        aux_sum += (matrix_r2(n,l,N,L, n-1,l+1,N,L-1, lamda)*
                    BM_Bracket(n-1,l+1,N,L-1, n1-1,l1,n2,l2, lamda))
        aux_sum += (matrix_r2(n,l,N,L, n,l-1,N-1,L+1, lamda)*
                    BM_Bracket(n,l-1,N-1,L+1, n1-1,l1,n2,l2, lamda))
        aux_sum += (matrix_r2(n,l,N,L, n,l-1,N,L-1, lamda)*
                    BM_Bracket(n,l-1,N,L-1, n1-1,l1,n2,l2, lamda))
        
        return np.sqrt(1./(n1*(n1+l1+0.5))) * aux_sum
        
    return ValueError
