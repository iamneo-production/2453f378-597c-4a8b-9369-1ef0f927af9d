import random as rand
import numpy as np
import cv2
from GA_pack import GA_pack as ga
from GA_help import GA_help as gahelp
from secrets import randbelow


# GA helper Package 
class GA_Main: 
   
    def combine_offspring(nc,ip_data,nobj,fit_type, mpool, umpool, fit_p):
        
        newpool=[]
        newpop=[]
        newcrowd_d=[]
        lmt=nc
        lmt2=nc
        fitpool1=[]
        chrom_num=lmt+lmt2
        fitpool1=ga.fitcal(umpool, lmt2, ip_data)
        fit1=fitpool1[0]
        fit2=fitpool1[1]
        fit1=fit1.tolist()
        fit2=fit2.tolist()
        fitpool1=[]
        fit1_h=[]
        fit2_h=[]
        i=0
        while i<lmt2:
            fit1_h.append(fit1[i])
            fit2_h.append(fit2[i])
            i+=1
        i=0
        while i<lmt:
            #child.append(parent(i))
            umpool.append(mpool[i])
            fit1.append(fit_p[0][i])
            fit2.append(fit_p[1][i])
            fit1_h.append(fit_p[0][i])
            fit2_h.append(fit_p[1][i])
            i+=1   
        #print('umpool legth',len(umpool))
        i=0
        fitpool1=[]
        fitpool1.append(fit1)
        fitpool1.append(fit2)
        fit_handle=[]
        fit_handle.append(fit1_h)
        fit_handle.append(fit2_h)
        front=[]
        front,Sp=ga.nondomsort(chrom_num,fitpool1)
        
        i=0
        j=0
        #INF=float(inf)
        new_fit=[]
        crowd_d,abc=ga.crowd_dist(fitpool1,front,nobj,fit_type)
        fit1=[]
        fit2=[]
        #print('Fit=', len(front[j]))
        while len(front[j])+len(newpop)<=lmt:
            temp=front[j]
            k=0
            j+=1
            lmt3=len(temp)
            #print('lmt3=',lmt3)
            while k<lmt3:
                newpool.append(umpool[temp[k]])
                newpop.append(temp[k])
                fit1.append(fit_handle[0][temp[k]])
                fit2.append(fit_handle[1][temp[k]])
                #newpop_front.append(j-1)
                newcrowd_d.append(crowd_d[temp[k]])
                k+=1
        temp=front[j]
        
        count=0
        handle_d=[]
        handle_c=[]
        while count<len(temp):
            handle_d.append(crowd_d[temp[count]])
            handle_c.append(temp[count])
            count+=1
        handle_c=np.array(handle_c)
        handle_d=np.array(handle_d)
        newfront=gahelp.insertion_sort_rev(handle_d,handle_c)
        l=0
        while len(newpool) < nc:
            newpool.append(umpool[newfront[l]])
            fit1.append(fit_handle[0][newfront[l]])
            fit2.append(fit_handle[1][newfront[l]])
            newpop.append(newfront[l])
            #newpop_front.append(j)
            newcrowd_d.append(crowd_d[newfront[l]])
            l+=1
        fit1=np.array(fit1)
        fit2=np.array(fit2)
        new_fit.append(fit1)
        new_fit.append(fit2)
        return new_fit,newpop, newcrowd_d, newpool
    
    def CrowdedbinTournamentSelec(pop,pop_front,crowd_d,nchrom):
        i=0
        low=0
        up=len(pop)-1
        npool=[]
        while i<nchrom:
            c1Pos=randbelow(up)
            c2Pos=randbelow(up)
            if pop_front[c1Pos] < pop_front[c2Pos]:
                npool.append(pop[c1Pos])
            elif pop_front[c1Pos] > pop_front[c2Pos]:
                npool.append(pop[c2Pos])
            else:
                if crowd_d[c1Pos]>crowd_d[c2Pos]:
                    npool.append(pop[c1Pos])
                else:
                    npool.append(pop[c2Pos])
            i+=1
        return npool
    def binTournamentSelec(pop,pop_front,nchrom):
        #mPool=[]
        i=0
        low=0
        up=len(pop)-1
        npool=[]
        while i<nchrom:
            c1Pos=randbelow(up)
            c2Pos=randbelow(up)
            if pop_front[c1Pos] < pop_front[c2Pos]:
                npool.append(pop[c1Pos])
            elif pop_front[c1Pos] > pop_front[c2Pos]:
                npool.append(pop[c2Pos])
            else:
                    temp=rand.randint(1,2)
                    if temp==1:
                        npool.append(pop[c1Pos])
                    else:
                        npool.append(pop[c2Pos])

            #print('c1Pos=',c1Pos)
            #print('c2Pos=',c2Pos)
            i+=1
        return npool   
    def crossover(n1pool):
        npool=n1pool
        low=0
        up=len(npool)-1
        i=0
        lt=int(len(npool)/2)+1
        while i<lt:
            crossprob=rand.random()
            if crossprob<0.8:
                ch1=rand.randint(low,up)
                ch2=rand.randint(low,up)
                up1=len(npool[ch1])-1
                up2=len(npool[ch2])-1
                if len(npool[ch1])>2:
                    low1=int(len(npool[ch1])/2)
                else:
                    low1=up1
                if len(npool[ch2])>2:
                    low2=int(len(npool[ch2])/2)
                else: 
                    low2=up2
                if up1!=low1 and up1> low1: 
                    crosspoint1=rand.randrange(low1,up1)
                else:
                    crosspoint1=up1+1
                if up2!=low2 and up2>low2:
                    crosspoint2=rand.randrange(low2,up2)
                else:
                    crosspoint2=up2+1
                if crosspoint1 %2 !=0:
                    crosspoint1+=1
                if crosspoint2 %2 !=0:
                    crosspoint2+=1
                j=crosspoint1
                temp1=npool[ch1]
                temp1=temp1.tolist()
                temp3=[]
                temp4=[]
                if j<up1:
                    while j<=up1:
                        temp3.append(temp1.pop(j))
                        up1-=1
                k=crosspoint2
                #print(k)
                up2=len(npool[ch2])-1
                #print(up2)
                temp2=npool[ch2]
                temp2=temp2.tolist()
                #print(temp2)
               
                if k<up2:
                    while k<=up2:
                        #print(k)
                        temp4.append(temp2.pop(k))
                        #print(temp4)
                        up2-=1
                if len(temp3)>1:
                    j=0
                    while j<len(temp3):
                        temp2.append(temp3[j])
                        j+=1
                if len(temp4)>1:
                    k=0
                    while k<len(temp4):
                        temp1.append(temp4[k])
                        k+=1
                #choice=rand.randint(1,10)
                npool[ch1]=np.zeros(len(temp1))
                temp1=np.array(temp1)
                npool[ch2]=np.zeros(len(temp2))
                temp2=np.array(temp2)
                npool[ch1]=temp1
                npool[ch2]=temp2

            #print(i)
            i+=1
        return(npool)
          
    def mutation(n1pool,row,ip_data):
        npool=n1pool
        low1=0
        up1=len(npool)-1
        i=0
        lt=int(len(npool)/2)+1
        while i<lt:
            ch1=rand.randrange(low1,up1)
            #print(ch1)
            temp1=[]
            temp1=npool[ch1]
            #temp1=temp1.tolist()
            #j=0
            #while j<len(npool[ch1]):
                #temp1.pop(j)
                #j+=1
            #print(npool[ch1])
            low=min(npool[ch1])
            #print(low)
            up=max(npool[ch1])
            #print(up)
            j=0
            while j<len(temp1):
                muteprob=rand.random()
                if muteprob<0.1:
                    r=rand.randrange(0,row-1)
                    temp1[j]=ip_data[r]
                j+=1
            npool[ch1]=np.zeros(len(temp1))
            temp1=np.array(temp1)
            npool[ch1]=temp1
            i+=1
        return(npool)  
        
    def segmentation(pmat,row,col,noclust):
        nc=0
        pmatr = np.zeros((row,col),np.uint8)
        pmatg=np.zeros((row,col),np.uint8)
        pmatb=np.zeros((row,col),np.uint8)
        cval=int(255/noclust)
        color=cval
        #colors = np.array([x for x in 'bgrcmykbgrcmykbgrcmykbgrcmyk'])
        #colors = np.hstack([colors] * 20)
        temp=0
        while nc<noclust:  
            i=0
            while i<row:
                j=0
                while j<col:
                    if pmat[i][j]==nc:
                        pmatr[i][j]=color
                        #pmatg[i][j]=color
                        #pmatb[i][j]=color
                    j+=1
                i+=1
            #temp+=3
            color+=cval
            #print(nc)
            nc+=1
        #pmat3 = np.zeros((row,col),np.uint8)
        #pmat3[:,:,0]=pmatr
        #pmat3[:,:,1]=pmatg
        #pmat3[:,:,2]=pmatb
        #img=hsv2rgb(pmat)
        cv2.imwrite('Segmented_Image.tif',pmatr)
        return pmatr
                    