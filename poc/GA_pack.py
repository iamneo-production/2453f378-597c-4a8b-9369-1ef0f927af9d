# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 21:34:05 2019

@author: Ramen Pal
"""

import numpy as np
import random2
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
   
    def encode(chrom_num, u_lim, l_lim, ip_data):
        mpool=[]
        lmt=len(ip_data)
        i=0
        while i < chrom_num: 
            c=random2.randrange(l_lim,u_lim,1)
            j=0
            a =np.zeros(c)
            while j < c: 
                r=int(random2.randrange(0,lmt,1))
                a[j]=ip_data[r]
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
    def fitcal(mpool, chrom_num, ip_data):
        i=0
        ln_data=len(ip_data)
        intradis=[]
        interdis=[]
        fitpool=[]
        #print('I am in fitcal')
        ###############################################################
        ###############################################################
        ############ Partition Matrix Generation ######################
        ###############################################################
        ###############################################################        
        while i<chrom_num :
            partition=np.zeros((ln_data),np.uint8)
            chrom=[]
            chrom=mpool[i]
            lmt=len(chrom)
            k=0
            while k< ln_data:
                l=0
                eudis=[]
                while l < lmt:
                    if ip_data[k]>chrom[l]:
                        a1=ip_data[k]-chrom[l]
                    else:
                        a1=chrom[l]-ip_data[k]
                    temp=np.sqrt((a1**2))
                    eudis.append(temp)
                    me=eudis.index(min(eudis))
                    l = l+1
                        #print(l)
                if eudis != []:  
                    me=eudis.index(min(eudis))
                    partition[k]=me

                eudis=0
                k +=1
            
            ###############################################################
            ###############################################################
            ############ Intra Cluster Distence ###########################
            ###############################################################
            ############################################################### 
            count=0
            peudis=[]
            #print('I am in Intra Cluster')
            while count < len(chrom):
                #print(count)
                cluster=np.where(partition == count) 
                #print(cluster)
                #print(len(cluster[0]))
                if len(cluster[0])>0:
                    #print('Cluster=',cluster)
                    #print('chrom=',chrom)
                    #print(cluster)
                    #if len(cluster[1]>1):
                    ncl=len(cluster)
                    #print(ncl)
                    #print(ncl)
                    nc=0
                    r=[]
                    while nc < ncl:
                        m1=cluster[0][nc]
                        r.append(ip_data[m1])
                        nc +=1
                         
                        #print(r)
                    edist=np.std(r)
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
            count1=0
            peudis=[]
            lmtc=len(chrom)
            while count < lmtc:
                if count1<lmtc:
                    count1+=count
                else:
                    count1=count
                while count1 < lmtc:
                    p1=chrom[count1]
                    p2=chrom[count]
                    if p1>p2:
                        a1=p1-p2
                        a1=a1**2
                    else:
                        a1=p2-p1
                        a1=a1**2

                    temp=np.sqrt(a1)
                    peudis.append(temp)
                    count1 +=2
                count +=2
            if len(chrom)>1:
                tinter=sum(peudis)/len(peudis)
            else:
                tinter=0.1
            #print('tinter',tinter)
            interdis.append(tinter)
            i +=1
            #print(intradis)
        lmt=len(intradis)
        temp1=np.zeros(lmt)
        temp2=np.zeros(lmt)
        i=0
        intradis=np.array(intradis)
        interdis=np.array(interdis)
        while i<lmt:
            temp1[i]=intradis[i]
            temp2[i]=interdis[i]
            i+=1
        fitpool.append(temp1)
        fitpool.append(temp2)
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
       dom_count=[]
       front=[]
       while i<chrom_num:
           Np=0
           j=0
           spt=[]
           while j<chrom_num:
               if i !=j and j not in spt:
                   if (fit[0][i] < fit[0][j]) and (fit[1][i] > fit[1][j]):
                       spt.append(j)
                   if (fit[0][i] > fit[0][j]) and (fit[1][i] < fit[1][j]):
                       Np+=1
               j+=1
           dom_count.append(Np)
           Sp.append(spt)
           i+=1  
       
       count=0
       while count<chrom_num:
           i=0
           hfront=[]
           while i<chrom_num:
               if dom_count[i]==0:
                   hfront.append(i)
                   dom_count[i]=99999
                   count+=1
               i+=1
           front.append(hfront)
           j=0
           lmt=len(hfront)
           while j< lmt:
               temp=Sp[hfront[j]]
               if temp!=[]:
                   lmt2=len(temp)
                   k=0
                   while k<lmt2:
                       dom_count[temp[k]]-=1
                       k+=1
               j+=1
       pop_front=dom_count
       i=0
       while i< len(front):
            hf=front[i]
            j=0
            while j< len(hf):
                pop_front[hf[j]]=i
                j+=1
            i+=1
       return( front, pop_front)
    def nondomsort2(chrom_num, fit):
       i=0
       ###############################################################
       ############ Find Domination count(Np) & Sp ###################
       ###############################################################
       Sp=[]
       dom_count=[]
       front=[]
       while i<chrom_num:
           Np=0
           j=0
           spt=[]
           while j<chrom_num:
               if i !=j and j not in spt:
                   if (fit[0][i] < fit[0][j]) and (fit[1][i] > fit[1][j]):
                       spt.append(j)
                   if (fit[0][i] > fit[0][j]) and (fit[1][i] < fit[1][j]):
                       Np+=1
               j+=1
           dom_count.append(Np)
           Sp.append(spt)
           i+=1  
       
       count=0
       while count<chrom_num:
           i=0
           hfront=[]
           while i<chrom_num:
               if dom_count[i]==0:
                   hfront.append(i)
                   dom_count[i]=99999
                   count+=1
               i+=1
           front.append(hfront)
           j=0
           lmt=len(hfront)
           while j< lmt:
               temp=Sp[hfront[j]]
               if temp!=[]:
                   lmt2=len(temp)
                   k=0
                   while k<lmt2:
                       dom_count[temp[k]]-=1
                       k+=1
               j+=1
       pop_front=dom_count
       i=0
       while i< len(front):
            hf=front[i]
            j=0
            while j< len(hf):
                pop_front[hf[j]]=i
                j+=1
            i+=1
       temp=[]
       temp=front[0]
       lmt=len(temp)
       first=[]
       i=0
       j=0
       count=0
       while i<lmt:
           count=0
           j=i+1
           while j<lmt:
               if i!=j:

                   if (fit[0][temp[i]] == fit[0][temp[j]]) or (fit[1][temp[i]] == fit[1][temp[j]]):
                       count+=1
               j+=1
           if count==0:
               first.append(i)
           i+=1
       
       front[0]=first
       return( front, pop_front)
       
    def plot_op(frnt,fitpool):
        #Sp=front1[0]
        #print("Sp")
        #print(Sp)
        #plt.plot(p_front[0], p_front[1])
        j=0
        intra=[]
        inter=[]
        fit1=fitpool[0]
        fit2=fitpool[1]
        #fit1 = [i * -1 for i in fit1]
        #fit2 = [j * -1 for j in fit2]
        temp1=fit1
        temp2=fit2
        lmt=len(frnt)
        while j<lmt:
            intra.append(temp1[frnt[j]])
            inter.append(temp2[frnt[j]])
            j+=1 
        #x_pareto = np.asarray(intra)
        #print(x_pareto)
        #y_pareto = np.asarray(inter)
        #plt.plot(x_pareto, y_pareto, color='r')
        #plt.show()
        #plt.hold(True)
        #function1 = [i for i in intra]
        #function2 = [j for j in inter]
        plt.xlabel('Intra Cluster Distance', fontsize=15)
        plt.ylabel('Inter Cluster Distance', fontsize=15)
        #print(intra)
        #print(inter)
        plt.scatter(intra, inter)
        plt.show()
        #plt.hold(True)
        return plt

    
    def crowd_dist(fitpool, front_full, nobj, fit_type):
        pt=0
        crowd_dist=[0.000]*len(fitpool[0])
        lmt=len(front_full)
        while pt<lmt:
            
            #felem=[]
            front=front_full[pt]
            obj=0
            front_hold=[]
            numfront=len(front)
            i=0
            while i<numfront:
                front_hold.append(front[i])
                i+=1
        
            fitlen=len(fitpool[0])
            INF=float(inf)
            while obj<nobj:
                i=0
                temp2=[]
                temp=[]
                fit_c=0
                while fit_c<fitlen:
                    temp.append(fitpool[obj][fit_c])
                    fit_c+=1
                while i<numfront:
                    temp2.append(temp[front[i]])
                    i+=1
                if fit_type[obj]==0:
                    front=ghelp.insertion_sort_rev(temp2,front_hold)
                else:
                    front=ghelp.insertion_sort(temp2,front_hold)
                
                if max(temp2)!=min(temp2):
                    crowd_dist[front[0]]=INF
                    crowd_dist[front[numfront-1]]=INF

                j=1
                INF=float(inf)
                while j < numfront-1:
                    if crowd_dist[front[j]]!=INF and max(temp2)!=min(temp2):
                        if temp[front[j-1]] > temp[front[j+1]]:
                            diff=temp[front[j-1]]-temp[front[j+1]]
                        else:
                            diff=temp[front[j+1]]-temp[front[j-1]]
                        crowd_dist[front[j]]+= abs((diff/(max(temp2)-min(temp2))))
                    if max(temp2)==min(temp2):
                        crowd_dist[front[j]]+=0
                    j +=1
                obj+=1
            pt+=1
        return crowd_dist, front
    
    def part_matrix(mpool, chrom_num, ip_data):
        ###############################################################
        ###############################################################
        ############ Partition Matrix Generation ######################
        ###############################################################
        ###############################################################  
        ln_data=len(ip_data)
        partition=np.zeros((ln_data),np.uint8)
        chrom=[]
        chrom=mpool
        lmt=len(mpool)
        k=0
        while k< ln_data:
            l=0
            eudis=[]
            while l < lmt:
                if ip_data[k]>chrom[l]:
                    a1=ip_data[k]-chrom[l]
                else:
                    a1=chrom[l]-ip_data[k]
                temp=np.sqrt((a1**2))
                eudis.append(temp)
                me=eudis.index(min(eudis))
                l = l+1
                #print(l)
            if eudis != []:  
                me=eudis.index(min(eudis))
                partition[k]=me
            eudis=0
            k +=1
        return(partition)
        
        
    def dunni(clust_cent1,clust_cent2,chro,partition,image, band):
        i=0
        j=1
        intradist=[]
        interdist=[]
        lm=int(len(chro)/2)
        while i<lm:
            j+=i
            while j<lm:
                if i!=j:
                    #print(clust_cent1[i],clust_cent2[i],clust_cent1[j],clust_cent2[j])
                    interdist.append(ghelp.inter_clust(chro,clust_cent1[i],clust_cent2[i],clust_cent1[j],clust_cent2[j],partition,image, band))
                j+=1
            i+=1
        i=0
        lm=int(lm/2)
        while i<lm:
            intradist.append(ghelp.intra_clust(i,partition,image, band))
            i+=1
        
        a=min(interdist)
        b=min(intradist)
        b2=max(intradist)
        c=max(a,b)
        dunn=a/b2
        if a>b:
            s_score=(a-b)/c
        else:
            s_score=(b-a)/c
        return dunn, s_score
    
    def dbi(clust_cent1,clust_cent2,chro,partition,image, band):
        i=0
        j=0
        db=0
        lm=int(len(chro)/2)
        
        while i<lm:
            db_temp=[]
            j+=i
            while j<lm:
                if i!=j:
                    a=ghelp.intra_clust(i,partition,image, band)
                    b=ghelp.intra_clust(j,partition,image, band)
                    #print(clust_cent1[i],clust_cent2[i],clust_cent1[j],clust_cent2[j])
                    #print(i,j)
                    c=ghelp.inter_clust(chro,clust_cent1[i],clust_cent2[i],clust_cent1[j],clust_cent2[j],partition,image, band)
                    #print(c)
                    temp=(a+b)/c
                    db_temp.append(temp)
                j+=1
            #print(db_temp)
            temp=0
            if len(db_temp)!=0:
                temp=max(db_temp)
            db+=temp
            i+=1
        #print(db)
        db=db/lm
        return db
    
    def silhouette_score(clust_cent1, clust_cent2,chro, partition, image, band):
        i=0
        j=1
        intradist=[]
        interdist=[]
        lm=int(len(chro)/2)
        while i<lm:
            j+=i
            while j<lm:
                if i!=j:
                    #print(clust_cent1[i],clust_cent2[i],clust_cent1[j],clust_cent2[j])
                    interdist.append(ghelp.inter_clust(chro,clust_cent1[i],clust_cent2[i],clust_cent1[j],clust_cent2[j],partition,image, band))
                j+=1
            i+=1
        i=0
        lm=int(lm/2)
        while i<lm:
            intradist.append(ghelp.intra_clust(i,partition,image, band))
            i+=1
        a=min(interdist)
        b=min(intradist)
        c=max(a,b)
        if a>b:
            s_score=(a-b)/c
        else:
            s_score=(b-a)/c
        return s_score

    def xbi(clust_cent1,clust_cent2,chro,partition,image, band):
        i=0
        j=0
        db=0
        lm=int(len(chro)/2)
        
        while i<lm:
            db_temp=[]
            j+=i
            while j<lm:
                if i!=j:
                    a=ghelp.intra_clust(i,partition,image, band)
                    b=ghelp.intra_clust(j,partition,image, band)
                    #print(clust_cent1[i],clust_cent2[i],clust_cent1[j],clust_cent2[j])
                    #print(i,j)
                    c=ghelp.inter_clust(chro,clust_cent1[i],clust_cent2[i],clust_cent1[j],clust_cent2[j],partition,image, band)
                    #print(c)
                    temp=(a+b)/c
                    db_temp.append(temp)
                j+=1
            #print(db_temp)
            temp=0
            if len(db_temp)!=0:
                temp=max(db_temp)
            db+=temp
            i+=1
        #print(db)
        xb=db/lm
        return xb
   