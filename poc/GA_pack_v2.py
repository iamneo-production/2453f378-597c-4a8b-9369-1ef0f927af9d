# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 21:34:05 2019

@author: Ramen Pal
"""

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
import statistics 
import matplotlib.pyplot as plt
from numpy import inf
from GA_help import GA_help as ghelp

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
        #print(mpool) 
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
        print('I am in fitcal')
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
            #print(partition)
            ###############################################################
            ###############################################################
            ############ Intra Cluster Distence ###########################
            ###############################################################
            ############################################################### 
            print('Partition matrix is done')
            count=0
            peudis=[]
            print('I am in Intra Cluster')
            while count < len(chrom)/2:
                #print(count)
                cluster=np.where(partition == count) 
                print(cluster)
                #if len(cluster[1]>1):
                ncl=len(cluster[1])
                print(ncl)
                #print(ncl)
                nc=0
                r=[]
                g=[]
                b=[]
                n=[]
                sd=[]
                while nc < ncl:
                    m1=cluster[0][nc]
                    m2=cluster[1][nc]
                    r.append(red[m1][m2])
                    g.append(green[m1][m2])
                    b.append(blue[m1][m2])
                    n.append(nir[m1][m2])
                    nc +=1
                         
                print(r)
                sd.append(np.std(r))
                sd.append(np.std(g))
                sd.append(np.std(b))
                sd.append(np.std(n))
                edist=sum(sd)/4
                peudis.append(edist)
                count +=1
    
            tintra=sum(peudis)/len(peudis)
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
       #print(chrom_num)
       front=[]
       while i<chrom_num:
           Np=0
           j=0
           spt=[]
           while j<chrom_num:
               if i !=j:
                   if (fit[i][0] < fit[j][0]) and (fit[i][1] > fit[j][1]):
                       Np+=1
                   if (fit[i][0] > fit[j][0]) and (fit[i][1] < fit[j][1]):
                       spt.append(j)
               j+=1
           Sp.append(Np)
           Sp.append(spt)
           i+=1  
       hfront=[]
       #print(Sp)
       count=0
       hfront, Sp=ghelp.front(Sp, chrom_num)
       count=len(hfront)
       #print(hfront)
       #print(hfront)
       front.append(hfront)
       while count<chrom_num:
           hfront1, Sp =ghelp.front_help(Sp, hfront)
           #print(hfront1)      
           front.append(hfront1)
           count +=len(hfront1)
           #print(count)
           if count>=chrom_num:
               break
       front2=[]
       front2=[z for z in front if z != []]
       print('I am in nondom sort')
       print(front)
       print('front2')
       print(front2)
       return( front2,Sp )
       
    def plot_op(fitpool, Sp):
        x = np.asarray(fitpool[0])
        y = np.asarray(fitpool[1])
        #plt.plot(p_front[0], p_front[1])
        j=0
        intra=[]
        inter=[]
        while j<len(Sp[0]):
            intra.append(fitpool[0][0])
            inter.append(fitpool[1][0])
            j+=1 
        #x_pareto = np.asarray(intra)
        #print(x_pareto)
        #y_pareto = np.asarray(inter)
        #plt.plot(x_pareto, y_pareto, color='r')
        #plt.show()
        #plt.hold(True)
        i=0 
        #plt.scatter(x, y)
        plt.xlabel('Intra-cluster')
        plt.ylabel('Inter-cluster')
        j=0
        intra=[]
        inter=[]
        while j<len(Sp[i]):
            intra.append(fitpool[0][j])
            inter.append(fitpool[1][j])
            j+=1
        print(intra)
        print(inter)
        x_pareto = np.asarray(intra)
        y_pareto = np.asarray(inter)
        #while i<len(Sp):
            #df = pd.DataFrame(Sp[i])
            #df.sort_values(0, inplace=True)
            #pareto_front = df.values
            #j=0
            #intra=[]
            #inter=[]
            #while j<len(Sp[i]):
                #intra.append(fitpool[0][j])
                #inter.append(fitpool[1][j])
                #j+=1
            #print(intra)
            #print(inter)
            #x_pareto = np.asarray(intra)
            #y_pareto = np.asarray(inter)
            #ptf=np.column_stack((x_pareto, y_pareto)).T
            #pareto_front_df = pd.DataFrame(ptf)
            #pareto_front_df.sort_values(0, inplace=True)
            #ptf = pareto_front_df.values
            #x_pareto = ptf[:, 0]
            #y_pareto = ptf[:, 1]
        plt.scatter(x_pareto, y_pareto, color='r')
            #plt.hold(True)
            #i+=1
        plt.show()
    
    def crowd_dist(fitpool, front, nobj, fit_type):
        
        
        #felem=[]
        obj=0
        front_hold=front
        numfront=len(front)
        crowd_dist=[0.000]*len(fitpool[0])
        INF=float(inf)
        while obj<nobj:
            temp=[]
            temp=fitpool[obj]
            i=0
            temp2=[]
            while i<len(front):
                temp2.append(temp[front[i]])
                i+=1
            if fit_type[obj]==0:
                front=ghelp.insertion_sort_rev(temp2,front_hold)
            else:
                front=ghelp.insertion_sort(temp2,front_hold)
            crowd_dist[front[0]]=INF
            crowd_dist[front[numfront-1]]=INF
            j=0
            INF=float(inf)
            while j < numfront:
                if crowd_dist[front[j]]!=INF and max(temp2)!=min(temp2):
                    crowd_dist[front[j]]+= abs(((temp[front[j-1]]-temp[front[j+1]])/(max(temp2)-min(temp2))))
                if crowd_dist[front[j]]!=INF and max(temp2)==min(temp2):
                    crowd_dist[front[j]]+=temp[front[j-1]]/max(temp2)
                j +=1
            obj+=1
        return crowd_dist, front
    
    def part_matrix(mpool, chrom_num, band, row, col, img):
        red=rasterio.open(img).read(1)
        green=rasterio.open(img).read(2)
        blue=rasterio.open(img).read(3)
        nir=rasterio.open(img).read(4)
        ###############################################################
        ###############################################################
        ############ Partition Matrix Generation ######################
        ###############################################################
        ###############################################################        
        partition=np.zeros((row,col))
        chrom=[]
        chrom=mpool
        j=0
        while j< row:
            k=0
            while k< col:
                l=0
                eudis=[]
                while l < len(mpool):
                    #print(red[[j],[k]])
                    temp=np.sqrt(((red[[j],[k]]-red[[chrom[l]],[chrom[l+1]]])**2)+((green[[j],[k]]-green[[chrom[l]],[chrom[l+1]]])**2)+((blue[[j],[k]]-blue[[chrom[l]],[chrom[l+1]]])**2)+((nir[[j],[k]]-nir[[chrom[l]],[chrom[l+1]]])**2)) 
                    eudis.append(temp)
                    #print(temp)
                    l = l+2
                    #print(l)
                    if l>len(mpool):
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
        return(partition)