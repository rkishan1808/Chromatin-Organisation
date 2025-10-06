#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 17:26:09 2021

@author: ravikishan
"""

import numpy as np
import matplotlib.pyplot as plt
import math
r=1
m=1000
n=1000
nb=700
rb0=1
rd0=1
dna=100
protein=1
t=0
time=np.zeros((m,n))
Nb= np.zeros((m,n))
time1=np.zeros((m,n))
Nb1= np.zeros((m,n))
def kmcavg(a):
    for i in range (m):
        dna=100
        t=0
        nb=a
        for j in range(n):
            rb=rb0*((dna-(protein*nb)))
            rd=rd0*(protein*nb)
            pb=rb/(rb+rd)
            pd=rd/(rb+rd)
            r1=np.random.uniform(0,1)
            if r1<pb:
                nb=nb+1
            else:
                nb=nb-1
                if dna<100:
                    dna=dna+protein
            r2=np.random.uniform(0,1)
            t=t-(math.log10(r2)/(rb+rd))
            time[i][j]=t
            Nb[i][j]=nb
    time1=np.mean(time,axis=0)
#print(Nb)
    Nb1=np.mean(Nb,axis=0)
    return time1,Nb1


a=10
t1,n1=kmcavg(a)
plt.plot(t1,n1)
#plt.ylim(0,500)
plt.ylabel("Avg no. of cohesin bound")
plt.xlabel("Time")
plt.title("rb=rd=r r=1")


a=30
t1,n1=kmcavg(a)
plt.plot(t1,n1)
#plt.ylim(0,500)
plt.ylabel("Avg no. of cohesin bound")
plt.xlabel("Time")
plt.title("rb=rd=r r=1")


a=50
t1,n1=kmcavg(a)
plt.plot(t1,n1)
#plt.ylim(0,500)
plt.ylabel("Avg no. of cohesin bound")
plt.xlabel("Time")
plt.title("rb=rd=r r=1")


a=70
t1,n1=kmcavg(a)
plt.plot(t1,n1)
#plt.ylim(0,500)
plt.ylabel("Avg no. of cohesin bound")
plt.xlabel("Time")
plt.title("rb=rd=r r=1")

a=90
t1,n1=kmcavg(a)
plt.plot(t1,n1)
#plt.ylim(0,500)
plt.ylabel("Avg no. of cohesin bound")
plt.xlabel("Time")
plt.title("Avg")

a=0
t1,n1=kmcavg(a)
plt.plot(t1,n1)
#plt.ylim(0,500)
plt.ylabel("Avg no. of cohesin bound")
plt.xlabel("Time")
plt.title("rb=rd=r r=1")

a=100
t1,n1=kmcavg(a)
plt.plot(t1,n1)
#plt.ylim(0,500)
plt.ylabel("Avg no. of cohesin bound")
plt.xlabel("Time")
plt.title("rb0=1, rd0=1")