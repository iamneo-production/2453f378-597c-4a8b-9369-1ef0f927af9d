import numpy as np
import rasterio
import pandas as pd
import random2
import rasterio
from osgeo import gdal, osr
import rasterio.plot
import pyproj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# GA Package 
class GA_pack: 
    # First we create a constructor for this class 
    # and add members to it, here models 

    # A normal print function 
    ###############################################################
    ###############################################################
    ############ Initial Encoding of Chromosome ###################
    ################Variable Length Encoding#######################
    ###############################################################  
    def submatrix( matrix, row, col):
        i=j=0
        opm = [[0 for x in range(row)] for y in range(col)]
        while i< row:
            while j< col:
                opm[[i],[j]]=matrix[[i],[j]]
                j +=1
            i +=1
        return (opm)
   
    def encode(chrom_num, u_lim, l_lim, row, col):
        mpool=[]
        i=0
        while i < chrom_num: 
            c=random2.randrange(l_lim,u_lim,1)
            j=0
            a =[]
            while j < (c*2): 
                if j%2 > 0:
                    a.append(int(random2.randrange(0,col,1)))
                else:
                    a.append(int(random2.randrange(0,row,1)))
                j += 1 
            mpool.append(a) 
            i +=1    
        print(mpool) 
        return (mpool)
    
    ###############################################################
    ###############################################################
    ############ Fitness Evaluation ###############################
    ###############################################################
    ###############################################################
    def fitcal(mpool, chrom_num, band, row, col, img):
        i=0
        red=rasterio.open(img).read(1)
        green=rasterio.open(img).read(2)
        blue=rasterio.open(img).read(3)
        nir=rasterio.open(img).read(4)
        intradis=[]
        interdis=[]
        fitpool=[]
        ###############################################################
        ###############################################################
        ############ Partition Matrix Generation ######################
        ###############################################################
        ###############################################################        
        while i<chrom_num :
            partition=np.zeros((row,col))
            chrom=[]
            chrom=mpool[i]
            j=0
            while j< row:
               k=0
               while k< col:
                   l=0
                   eudis=[]
                   while l < len(mpool[i]):
                      #print(red[[j],[k]])
                      temp=np.sqrt(((red[[j],[k]]-red[[chrom[l]],[chrom[l+1]]])**2)+((green[[j],[k]]-green[[chrom[l]],[chrom[l+1]]])**2)+((blue[[j],[k]]-blue[[chrom[l]],[chrom[l+1]]])**2)+((nir[[j],[k]]-nir[[chrom[l]],[chrom[l+1]]])**2)) 
                      eudis.append(temp)
                      #print(temp)
                      l = l+2
                      #print(l)
                      if l>len(mpool[i]):
                          break
                   
                   partition[[j],[k]]=eudis.index(min(eudis))
                   #print(eudis)
                   eudis=0
                   k +=1
                   #print(k)
                   if k>col:
                       break
               j +=1
               #print(j)
               if j > row:
                   break
            print(partition)
            ###############################################################
            ###############################################################
            ############ Intra Cluster Distence ###########################
            ###############################################################
            ############################################################### 
            count=0
            peudis=[]
            while count < len(chrom)/2:
                cluster=np.where(partition == count) 
                #if len(cluster[1]>1):
                ncl=len(cluster[1])
                nc=0
                peudis1=[]
                while nc < ncl:
                       nc1=0
                       while nc1 < ncl:
                           if nc1 != nc:
                               temp=np.sqrt(((red[[cluster[0][nc]],[cluster[1][nc]]]-red[[cluster[0][nc1]],[cluster[1][nc1]]])**2)+((green[[cluster[0][nc]],[cluster[1][nc]]]-green[[cluster[0][nc1]],[cluster[1][nc1]]])**2)+((blue[[cluster[0][nc]],[cluster[1][nc]]]-blue[[cluster[0][nc1]],[cluster[1][nc1]]])**2)+((nir[[cluster[0][nc]],[cluster[1][nc]]]-nir[[cluster[0][nc1]],[cluster[1][nc1]]])**2))
                               peudis1.append(temp)
                            #print(temp)
                           nc1 +=1
                       nc +=1
                print(peudis1)
                #try:
                edist=sum(peudis1)/len(peudis1)
               # except ZeroDivisionError: 
                peudis.append(edist)
                #print (peudis)
                count +=1
    
            tintra=max(peudis)
            #print(tintra)
            intradis.append(tintra) 
            ###############################################################
            ###############################################################
            ############ Inter Cluster Distence ###########################
            ###############################################################
            ############################################################### 
            count=0
            peudis=[]
            while count < len(chrom):
                count1=0
                while count1 < len(chrom):
                    if count1 != count:
                        temp=np.sqrt(((red[[chrom[count1]],[chrom[count1+1]]]-red[[chrom[count]],[chrom[count+1]]])**2)+((green[[chrom[count1]],[chrom[count1+1]]]-green[[chrom[count]],[chrom[count+1]]])**2)+((blue[[chrom[count1]],[chrom[count1+1]]]-blue[[chrom[count]],[chrom[count+1]]])**2)+((nir[[chrom[count1]],[chrom[count1+1]]]-nir[[chrom[count]],[chrom[count+1]]])**2))
                        peudis.append(temp)
                    count1 +=2
                count +=2
            tinter=min(peudis)
            interdis.append(tinter)
            

            i +=1
            #print(intradis)
            if i> chrom_num:
                break
        fitpool.append(intradis)
        fitpool.append(interdis)
        return (fitpool)
                   
                 
    ###############################################################
    ###############################################################
    ############ Non Dominated Sorting ############################
    ###############################################################
    ###############################################################                     
    
    def nondomsort(chrom_num, fit):
       i=0
       ###############################################################
       ############ Find Domination count(Np) & Sp ###################
       ###############################################################
       Sp=[]
       print(chrom_num)
       while i<chrom_num:
           Np=0
           j=0
           spt=[]
           
           while j<chrom_num:
               if i !=j:
                   if (fit[0][i] < fit[0][j]) and (fit[1][i]>fit[1][j]):
                       #print(Np)
                       Np+=1
                       spt.append(j)
               j+=1
           Sp.append(Np)
           Sp.append(spt)
           i+=1
       return(Sp)
                       
                
            
        