#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 17:57:34 2018

@author: miguel
"""
from __future__ import division, print_function
import numpy as np
import BM_brackets as bm
import Matrix_Elements as me

def read_Antoine(index):
    
    ## returns the Quantum numbers from string Antoine's format
    aux = int(index)
    if(aux == 1):
        return[0,0,1]
    else:
        n = int((aux)/1000)
        l = int((aux - (n*1000))/100)
        j = int(aux - (n*1000) - (l*100))# is 2j 
        return [n,l,j]
    
    return [0,0,0]


def read_input_file(filename):
    ## Pasing of arguments to perform the program
    ##
    ## Form of the File requiered:
    
    #============== ...
    #INPUT FILE (use 1 tabulation to separate name from arguments)
    # do not change ...
    #============== ...
    #Interaction	Brink_Boeker
    #File_Name	bb_SP
    #A_Mass	2
    #H_Omega	23.8
    #b_length	1.3206
    #Q_Numbers	001	101	103
    #QN_Energies
    #Core		0	0
    #CENTRAL_mu	1.0	0.0
    #CENTRAL_Wign	-1.0	0.0
    #CENTRAL_Majo	0.0	0.0
    #CENTRAL_Barl	0.0	0.0
    #CENTRAL_Heis	0.0	0.0
    
    ## 
    ## if there are properties unnused, just remove all spaces in  the line
    
    INPUT = [np.array(map(str,line.split())) for line in open(filename) ]
    
    for i in range(4):
        INPUT.remove(INPUT[0])
    
    Q_Numbers = []
    Q_Numbers_str = []
    Q_Energies = []
    parameters = [[0. for i in range(7)],[0. for i in range(7)]]
    
    # Read Lines and asign/convert the parameters
    for line in INPUT:
        if line[0] == 'Interaction':
            Interaction = line[1]
            
        elif line[0] == 'File_Name':
            File_Name = line[1]
            
        elif line[0] == 'A_Mass':
            if len(line) == 1:
                A = 1
            else:
                A = float(line[1])
                
        elif line[0] == 'H_Omega':
            if len(line) == 1:
                HBAR_OMEGA = 41. / (A**(1./3))
            else:
                HBAR_OMEGA = float(line[1])
                
        elif line[0] == 'b_length':
            if len(line) == 1:
                b_param = 1.010 / (A**(1./6))
            else:
                b_param = float(line[1])
            parameters[0][0] = b_param
            parameters[1][0] = b_param
        
        elif line[0] == 'Q_Numbers':
            for  i in range(1,len(line)):
                if(line[i] != ''):
                    Q_Numbers.append(read_Antoine(line[i]))
                    Q_Numbers_str.append(line[i])
                else:
                    pass
                
        elif line[0] == 'QN_Energies':
            if len(line)>1:
                for  i in range(1,len(line)):
                    Q_Numbers.append(float(line[i]))
            else:
                for QN in Q_Numbers:
                    Q_Energies.append(HBAR_OMEGA*(2*QN[0]+QN[1]+1.5))
        
        elif line[0] == 'Core':
            if len(line)>1:
                Core = [ int(line[1]), int(line[2])]
            else:
                Core = [ 0 , 0 ]
                
        elif line[0] == 'BB_mu':
            parameters[0][1] = float(line[1])
            parameters[1][1] = float(line[2])
            
        elif line[0] == 'BB_Wign':
            parameters[0][2] = float(line[1])
            parameters[1][2] = float(line[2])
        elif line[0] == 'BB_Majo':
            parameters[0][3] = float(line[1])
            parameters[1][3] = float(line[2])
        elif line[0] == 'BB_Barl':
            parameters[0][4] = float(line[1])
            parameters[1][4] = float(line[2])
        elif line[0] == 'BB_Heis':
            parameters[0][5] = float(line[1])
            parameters[1][5] = float(line[2])
        
    # Selecting argumets by their applications
    Parameters = [Interaction,A,HBAR_OMEGA, parameters]
    Parameters_Output = [File_Name,Q_Numbers_str,Q_Energies, Core]
    Wave_Functions = [Q_Numbers]
        
    return [Wave_Functions, Parameters, Parameters_Output]




def run_antoine_output(filename_input):
    
    ## Main run procedure to produce a Antoine Format Matrix element file
    
    ## Import the data
    data = read_input_file(filename_input)
    
    Wave_Functions = data[0]    
    Parameters = data[1]
    Parameters_Output = data[2]

    Interaction = Parameters[0]
    parameters  = Parameters[3]
    
    File_Name = Parameters_Output[0]
    Aux_Name = File_Name.split('_')
    Shell = Aux_Name[1]
    Q_Numbers_str = Parameters_Output[1]
    Q_Energies = Parameters_Output[2]
    Core = Parameters_Output[3]
    
    Q_Numbers = Wave_Functions[0]
    
    number = len(Q_Numbers)
    
    ## Write Header
    f = open(File_Name, 'w')
    f.write('  '+Interaction+' INT. ('+Shell+' SHELL)\n')
    
    f.write(' 1 %d' %  number)
    for QN in Q_Numbers_str:
        f.write(' '+QN)
    f.write('\n   ')
    for QE in Q_Energies:
        f.write(' '+str(QE))
    f.write('\n')
    f.write('    0    '+str(Core[0])+'    '+str(Core[1])+
            '    0.333333    0.000000\n')
    
    ## Perform all the levels
    for ia in range(number):
        a=[Q_Numbers[ia][0],Q_Numbers[ia][1],
           Q_Numbers[ia][2],Q_Numbers_str[ia]]
        for ib in range(ia, number):
            b=[Q_Numbers[ib][0],Q_Numbers[ib][1],
               Q_Numbers[ib][2],Q_Numbers_str[ib]]
            
            Bra = me.QN_2body_jj_Coupling(a[0],a[1],a[2],b[0],b[1],b[2],
                                          0,0)
            j_bra_min=int(abs(a[2]-b[2])/2)
            j_bra_max=int((a[2]+b[2])/2)
            
            for ic in range(ia, number):
                c=[Q_Numbers[ic][0],Q_Numbers[ic][1],
                   Q_Numbers[ic][2],Q_Numbers_str[ic]]
                
                if(ia == ic):
                    id_min = max(ic,ib)
                else:
                    id_min = ic
                for i_d in range(id_min,number):
                    d=[Q_Numbers[i_d][0],Q_Numbers[i_d][1],
                       Q_Numbers[i_d][2],Q_Numbers_str[i_d]]
                    
                    Ket = me.QN_2body_jj_Coupling(c[0],c[1],c[2],d[0],d[1],d[2],
                                                  0,0)
                    
                    j_ket_min=int(abs(c[2]-d[2])/2)
                    j_ket_max=int((c[2]+d[2])/2)
                    
                    if((j_ket_min <= j_bra_max) and 
                       (j_bra_min <= j_ket_max)):
                        f.write(' 0 1 '+a[3]+" "+b[3]+" "+c[3]+" "+d[3]+" "+
                            str(max(j_bra_min,j_ket_min))+" "+
                            str(min(j_bra_max,j_ket_max))+"\n")
                        
                        ## T-J loops
                        for T in range(2):
                            for J in range(max(j_bra_min,j_ket_min),
                                           min(j_bra_max,j_ket_max)+1):
                                
                                Bra.J = J
                                Bra.T = T
                                Ket.J = J
                                Ket.T = T
                                
#                                print(a[3]+" "+b[3]+" "+c[3]+" "+d[3]+" JT: "+str(J)+" "+str(T))
                                
                                numeric=me.Brink_Boeker_ME(Bra, Ket, parameters)
                                
                                if(numeric < 0):
                                    if(abs(numeric)<1e-17):
                                        f.write('   0.00000')
                                    else:
                                        f.write('  '+str(round(numeric,5)))
                                else:
                                    if(abs(numeric)<1e-17):
                                        f.write('   0.00000')
                                    else:
                                        f.write('   '+str(round(numeric,5)))
                                    
                            f.write('\n')
                        
    
    f.close()    
    
