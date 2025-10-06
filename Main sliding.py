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


rb0=3
rd0=5
rs0=100
n_bound=0
nb_free=0
nf_free=0

t=0
length=100
dna=np.zeros(length)
p=np.zeros((length,length))
count=0
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
dna[30]=2
dna[70]=2
dna[55]=2
for i in range(100000):
    
    count1=0
    
     
    '''Available site counting'''
    available_sites=0
    for j2 in range(length-2):
        if dna[j2]==0 and dna[j2+1]==0:
            available_sites += 1
    
    
    #print("n_free",nfree)
    rb=rb0*(available_sites)
    rd=rd0*n_bound
    rs=rs0*n_bound
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
        #if NL>0:
    if NL>0:
        
        Avg_loop=L1/NL
        #print("NL",NL)
        #print("loop",loop)
        #print("L1",L1)
    
    #Ls_avg=Ls/n_bound
    
    #list_Ls_avg.append(Ls_avg)
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
    
    if i%1==0:
        for j in range(n_bound):
            p[sb[j]][sf[j]]+=1
            count+=n_bound
            p[sf[j]][sb[j]]= p[sb[j]][sf[j]]
    
    '''if i%10==0:
        file=open("CTCF2_position.xyz","a")
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
        
#print("Sb",sb)   
#print("sf",sf)
#print("L1",L1)
f=open("NOr_With_CTCF.txt","w")   
for pi in range(length):
    for pj in range(length):
        p[pi][pj]=((p[pi][pj])/count)
        p[sf[j]][sb[j]]= p[sb[j]][sf[j]]
        f.write("{} {} {}\n".format(pi,pj,p[pi][pj]))
f.close()
plt.imshow(p, cmap='BuPu', interpolation='nearest')
plt.xlabel("site j")
plt.ylabel("site i")
plt.title("With CTCF")
plt.colorbar()
plt.show()
plt.savefig("xyz.eps")
    
'''plt.plot(time,Nb)
plt.xlabel("Time")
plt.ylabel("No. of cohesin bound")
plt.title("rb0=1, rd0=1, rs0=1 & length=1000")
plt.show()

plt.plot(time,l_loop)
plt.xlabel("Time")
plt.ylabel("Length of loop")
plt.title("rb0=1, rd0=1, rs0=1 & length=1000")
plt.show()

plt.plot(time,n_loop)
plt.xlabel("Time")
plt.ylabel("No. of loops")
plt.title("rb0=1, rd0=1, rs0=1 & length=1000")
plt.show()

plt.plot(time,sub_loop)
plt.xlabel("Time")
plt.ylabel("No. of sub loops")
plt.title("rb0=1, rd0=1, rs0=1 & length=1000")
plt.show()

plt.plot(time,list_Avg_loop)
plt.xlabel("Time")
plt.ylabel("Avg length of each loops")
plt.show()

plt.plot(time,list_Ls_avg)
plt.xlabel("Time")
plt.ylabel("Avg loop length with sub_loops")
#plt.plot(time,l_loop)
#plt.imshow(dna, cmap='hot', interpolation='nearest')
#plt.xlabel("sites [ic+n]")
#plt.ylabel("sites [jc-n]")
#plt.title("Probability of Loop extrusion b/w two convergent CTCF") 
#plt.colorbar()
#plt.show()
'''
            