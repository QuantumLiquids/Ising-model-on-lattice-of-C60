#!/usr/bin/env python
# -*- encoding:utf-8 -*- 


# Calculate the density of states using Wang-Landau MC method
# Reference: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.2050
# 10.1103/PhysRevE.64.056101

# The sites are labeled from 1 to 60

# Acknowledgement: Wen-Jie Jiang help for the grammar

from math import exp
from math import sqrt
from math import log
import random
import numpy as np
import time



beta = 1 # reciprocal of tempreture
L=60        # The number of total sites
J = 1        # 1 for AF, and -1 for FM

# Para List
flat_para = 0.8
sample_time = 2000000
min_factor = 1+pow(10,-8)

def getEnergy(S):
    Energy = 0
    
    for i in range(5):
        Energy = Energy + J*S[i]*S[(i+1)%5] # the first layer
        Energy = Energy + J*S[i]*S[3*i+6]# bond between the 1st and 2nd layer
        Energy = Energy + J*S[3*i+5]*S[3*i+6] + J*S[3*i+6]*S[3*i+7]+J*S[3*i+7]*S[(3*i+8-5)%15+5]# the 2nd layer
        Energy = Energy + J*S[3*i+5]*S[4*(i+5)]+ J*S[3*i+7]*S[4*i+21]# bond between the 2nd and 3rd layer

        # The other half part can get from reflection
        Energy = Energy + J*S[59-i]*S[59-((i+1)%5)] # the last layer, with PBC
        Energy = Energy + J*S[59-i]*S[59-(3*i+6)]# bond between the 5st and 4nd layer
        Energy = Energy + J*S[59-(3*i+5)]*S[59-(3*i+6)] + J*S[59-(3*i+6)]*S[59-(3*i+7)]+J*S[59-(3*i+7)]*S[59-((3*i+8-5)%15+5)]# the 4nd layer
        Energy = Energy + J*S[59-(3*i+5)]*S[59-(4*(i+5))]+ J*S[59-(3*i+7)]*S[59-(4*i+21)]# bond between the 4nd and 3rd layer

    # the middle layer
    for i in range(20,39):
        Energy = Energy + J*S[i]*S[i+1]

    Energy = Energy + J*S[39]*S[20]        
        
    return Energy


if __name__ == "__main__":
    print("==============================================")
    print("Para List:")
    print("flat factor:\t ", flat_para)
    print("sample time uppper bound: \t", sample_time)
    print("modification factor lower bound:\t",min_factor)
    print("===============================================")
    
    En= [ i for i in range(-90,91,2)] # possible values for energy
    g = np.ones(len(En))              # initial density of state
    f = 2.7182818284590455348848081484903 # e

    time_start = time.time()
    for j in range(30):
        H = np.zeros(len(En))
        print('step',j,'f=',  f, end='')
        S = [random.randint(0,1)*2-1 for i in range(L)] # initial a random state
        S_flip=[0 for l in range(L)]
        for l in range(L):
               S_flip[l]=S[l]
        for i in range(sample_time):
            tag=0
            E1=getEnergy(S)
            
            k = random.randint(0,L-1)
            S_flip[k]=-S_flip[k]

            E2=getEnergy(S_flip)

            index1=int((E1+90)/2)
            index2=int((E2+90)/2)
            if(g[index1] >= g[index2] or (g[index1] < g[index2] and random.random()<=g[index1]/g[index2])):
               S[k]=-S[k]
               E=E2
            else:
               E=E1
               S_flip[k]=-S_flip[k]

            indexE=int((E+90)/2)
            if g[indexE]>pow(10,20):
                g=g/sum(g)*pow(10,18)
                
            g[indexE]=g[indexE]*f
            H[indexE]=H[indexE]+1

            if i>1000 and np.min(H[np.nonzero(H)])>=flat_para*np.mean(H[np.nonzero(H)]) :
                tag=1
                print('\t sweep_time=', i+1,end='')
                break
        if(tag==0):
            print('\t sweep_time=', i+1,end='')
            
            
        g=g/sum(g)*1152921504606846976
        print('\t GS degeneracy:', round(g[12]),end='')

        Z=0 # partition function
        for i in range(len(En)):
            Z=Z+g[i]*exp(-beta*En[i])

        print('\t ln(Z)/N=', log(Z)/L)
        if(f<=min_factor):
            break
        f=sqrt(f)


    g=g/sum(g)*1152921504606846976; # 1152921504606846976.0=2^60

    Z=0 # partition function
    for i in range(len(En)):
        Z=Z+g[i]*exp(-beta*En[i])

    print("===============================================")
    print('ln(Z)/N = \t', log(Z)/L)
    
    print('The possible energy level List: ')
    print(En)
    print('Corresponding energy density : ')
    print(np.rint(g))
    time_end = time.time()
    print("time: %f s" %(time_end-time_start))


    
