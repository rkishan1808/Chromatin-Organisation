#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 17:24:18 2021

@author: ravikishan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 02:16:47 2021

@author: ravikishan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 16:24:14 2021

@author: ravikishan
"""

#sliding
import numpy as np
import matplotlib.pyplot as plt
import math


rb0=1
rd0=1
rs0=1
n_bound=0
nb_free=0
nf_free=0

t=0
length=1000
dna=np.zeros(length)
time=[]
Nb=[]
l_loop=[]
sb=[]
sf=[]
available_sites=0
#dna[0]=2
#dna[999]=2
x=0
y=0
z=0
count1=0
L1=0
NL=0
n_loop=[]
loop=[]
sub_loop=[]
list_Avg_loop=[]
Avg_loop=0
Ls_avg=0
list_Ls_avg=[]

total_list_Ls_avg=[]
total_list_Avg_loop=[]
total_sub_loop=[]
total_l_loop=[]
total_n_loop=[]
total_time=[]
total_n_bound=[]

for j2 in range(100):
    dna=np.zeros(length)
    n_loop=[]
    loop=[]
    sub_loop=[]
    list_Avg_loop=[]
    time=[]
    Nb=[]
    l_loop=[]
    sb=[]
    sf=[]
    list_Ls_avg=[]
    n_bound=0
    t=0
    L1=0
    for i in range(10000):
        nb_free=0
        nf_free=0
        count1=0
        if n_bound>0:
            #print("sb",sb)
            #print("sf",sf)
            #print(n_bound)
            for j in range(n_bound):
                if sb[j]-1>=0 and dna[sb[j]-1]==0:
                    nb_free+=1
                if sf[j]+1<length and dna[sf[j]+1]==0:
                    nf_free+=1
         
        '''Available site counting'''
        available_sites=0
        for j2 in range(length-2):
            if dna[j2]==0 and dna[j2+1]==0:
                available_sites += 1
        
        nfree=nb_free + nf_free
        #print("n_free",nfree)
        rb=rb0*(available_sites)
        rd=rd0*n_bound
        rs=rs0*nfree
        pb=rb/(rb+rd+rs)
        pd=rd/(rb+rd+rs)
        ps=rs/(rb+rd+rs)
        r1=np.random.uniform(0,1)
        '''Binding'''
        if r1<=pb: 
           
            r2=np.random.randint(0,(length-2))
        
            #print("r2",r2)
            if dna[r2]==0 and dna[r2+1]==0:
               
                dna[r2]=1
                #print("x1",x1)
               
                dna[r2+1]=1
               
                sb.append(r2)
                #print("sb",sb)
                sf.append(r2+1)
                #print("sf",sf)
                #print("dna",dna)
                #print("sb",len(sb))
                #print("sf",len(sf))
                n_bound=n_bound+1
                #print("bound",n_bound)
                #print("nf_free",nf_free)
                
                
                
                         #print("nf_free",nf_free)
                    
            
            '''Dissociation'''
            
        elif pb<r1<=(pb+pd): 
            #len3=len(sf)
            r3=np.random.randint(0,n_bound)
            
            dna[sf[r3]]=0
            dna[sb[r3]]=0
            sf.remove(sf[r3])
            sb.remove(sb[r3])
            
            n_bound=n_bound-1
            #print("r3",r3)
            
            '''sliding'''
        else:
            
            r4=np.random.randint(0,n_bound)
            '''Forward'''
            if sf[r4]<length-1 and dna[sf[r4]+1]==0:
                
                dna[sf[r4]+1]=1
                dna[sf[r4]]=0
                sf[r4]+=1
               
               
                
            if sb[r4]>0 and dna[sb[r4]-1]==0:
                dna[sb[r4]-1]=1
                dna[sb[r4]]=0
                
                sb[r4]-=1
                    
                 
                   
                    #print("n",new)
            #loop=(r2+2)-(r2-1)
        #print("sb",sb)
        #print("sf",sf)
        '''Looping'''
        L1=0
        NL=0
        sub_NL=0
        loop=[]
        Ls=0
        for i3 in range (0,n_bound):
            if sf[i3]-sb[i3]>1:
                Ls+=sf[i3]-sb[i3]-1
            for i4 in range(sb[i3]+1,sf[i3]):
                if i4 in (loop):
                    if i4==(sb[i3]+1):
                        sub_NL+=1
                if i4 not in (loop):
                    loop.append(i4)
                    if i4==(sb[i3]+1):
                        NL+=1
               
                  
        L1=len(loop)
            
        if NL==0:
            Avg_loop=0
        else:
            Avg_loop=L1/NL
            #print("NL",NL)
            #print("loop",loop)
            #print("L1",L1)
        if n_bound==0:
            Ls_avg=0
        else:
            Ls_avg=(Ls/n_bound)
        
        list_Ls_avg.append(Ls_avg)
            
        list_Avg_loop.append(Avg_loop)
        sub_loop.append(sub_NL)
        l_loop.append(L1)
        n_loop.append(NL)
        #print("DNA",dna)
        r7=np.random.uniform(0,1)
        t=t-(math.log(r7)/(rb+rd+rs))
        #l_loop.append(loop)
        time.append(t)
        Nb.append(n_bound)
        
        '''if i%1==0:
            file=open("6position.xyz","a")
            file.write("{}\n{}\n".format(length,i))
            for i2 in range(length):
                file.write("{} {} {} {}\n".format(int(dna[i2]),x+i2,y,z))'''
        for i3 in range(length):
            if dna[i3]==1:
                count1+=1
                #print(i3)
    
        if n_bound!=count1/2:
            print("Error", i)
            print("Nbound",n_bound)
            print("Count",count1)
            print("sb",sb)
            print("sf",sf)
            
    total_list_Ls_avg.append(list_Ls_avg)
    total_list_Avg_loop.append(list_Avg_loop)
    total_sub_loop.append(sub_loop)
    total_l_loop.append(l_loop)
    total_n_loop.append(n_loop)
    total_time.append(time)
    total_n_bound.append(Nb)
        
mean_avg_ls=np.mean(total_list_Ls_avg,axis=0)
mean_avg_loop_length=np.mean(total_list_Avg_loop,axis=0)
mean_subloops=np.mean(total_sub_loop,axis=0)
mean_total_loop_length=np.mean(total_l_loop,axis=0)
mean_no_of_loops=np.mean(total_n_loop,axis=0)
mean_total_time=np.mean(total_time,axis=0)
mean_n_bound=np.mean(total_n_bound,axis=0)


plt.plot(mean_total_time,mean_n_bound)
plt.ylabel("Avg no. of bound cohesin")
plt.xlabel("Time")
plt.title("rb0=1,rd0=1,rs0=1 & length=1000")
plt.show()
  
plt.plot(mean_total_time,mean_avg_ls)
plt.xlabel("Time")
plt.ylabel("Avg length of each loops when subloops included ")
plt.title("rb0=1,rd0=1,rs0=1 & length=1000")
plt.show()

plt.plot(mean_total_time,mean_avg_loop_length)
plt.ylabel("avg loop length only longer loop")
plt.xlabel("Time")
plt.title("rb0=1,rd0=1,rs0=1 & length=1000")
plt.show()

plt.plot(mean_total_time,mean_total_loop_length)
plt.xlabel("Time")
plt.ylabel("Avg loop lenth extruded")
plt.title("rb0=1,rd0=1,rs0=1 & length=1000")
plt.show()

plt.plot(mean_total_time,mean_no_of_loops)
plt.xlabel("Time")
plt.ylabel("Avg no. of loops")
plt.title("rb0=1,rd0=1,rs0=1 & length=1000")
plt.show()
