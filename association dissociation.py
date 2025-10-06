#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 14:10:51 2021

@author: ravikishan
"""

import numpy as np
import matplotlib.pyplot as plt
import math
n=1000000
rb0=1
rd0=2
nb=100
t=0
dna=1000
d1=dna
protein=1
time=[]
Nb=[]

for i in range (n):
    rb=rb0*(dna-(protein*nb))
    #print(dna)
    #print("Np",nb)
    rd=rd0*(protein*nb)
    pb=rb/(rb+rd)
    pd=rd/(rb+rd)
    #print("pb",pb)
   # print("pd",pd)
    r1=np.random.uniform(0,1)
    #print()
    #print()
    if r1<pb:
        nb=nb+1
        
        #print("DNA",dna)
    else:
        nb=nb-1
        
        #print("DNA2",dna)
    r2=np.random.uniform(0,1)
    t=t-(math.log(r2)/(rb+rd))
    Nb.append(nb)
    time.append(t)

plt.plot(time,Nb)
plt.ylabel("No. of bound protein")
plt.xlabel("Time")
plt.title("rb0=1 & rd0=2")