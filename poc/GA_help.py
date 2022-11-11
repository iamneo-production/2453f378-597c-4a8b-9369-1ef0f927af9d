import numpy as np
import rasterio
from numpy import inf
import math

# GA helper Package 
class GA_help: 
    
    def front(Sp, chrom_num):  
       j=0
       hfront=[]
       #print(chrom_num)
       #print(len(Sp))
       while j< chrom_num:
           if Sp[j*2]==0:
               hfront.append(j)
               Sp[j*2]=None
           j +=1      
       return hfront, Sp
    
    def front_help(Sp, hfront):  
       i=0
       temp=0
       handle=[]
       hfront1=[]
       while i< len(hfront):
           temp=(hfront[i]*2)+1
           handle=Sp[temp]
           j=0
           while j< len(handle):
               temp2=(handle[j]*2)
               if temp2 not in hfront1 and Sp[temp2]!=None:
                   Sp[temp2]=Sp[temp2]-1
                   print(Sp[temp2])
               if Sp[temp2]==0:
                   hfront1.append(int(temp2/2))
                   Sp[temp2]=None
               j +=1
           i +=1      
       return hfront1, Sp
   
    def insertion_sort(nums, front):
            # Start on the second element as we assume the first element is sorted
            nums=np.array(nums)
            front=np.array(front)
            i=0
            while i < len(nums):
                item_to_insert = nums[i]
                item=front[i]
                # And keep a reference of the index of the previous element
                j = i - 1
                # Move all items of the sorted segment forward if they are larger than
                # the item to insert
                while j >= 0 and nums[j] > item_to_insert:
                    nums[j + 1] = nums[j]
                    front[j+1]=front[j]
                    j -= 1
                # Insert the item
                nums[j + 1] = item_to_insert
                front[j+1]=item
                i+=1
            front=front.tolist()
            return front
    

    def insertion_sort_rev(nums, front):

            # Start on the second element as we assume the first element is sorted
            i=0
            nums=np.array(nums)
            front=np.array(front)
            while i < len(nums):
                item_to_insert = nums[i]
                item=front[i]
                # And keep a reference of the index of the previous element
                j = i - 1
                # Move all items of the sorted segment forward if they are larger than
                # the item to insert
                while j >= 0 and nums[j] > item_to_insert:
                    nums[j + 1] = nums[j]
                    front[j+1]=front[j]
                    j -= 1
                # Insert the item
                nums[j + 1] = item_to_insert
                front[j+1]=item
                i+=1
            nums=nums.tolist()
            nums.reverse()
            front=front.tolist()
            front.reverse()
            #newf.reverse()
            return front
        

    def duplicate_elim(mpool,nc):
        ch=0
        while ch<nc:
            chrom=mpool[ch]
            chrom=chrom.tolist()
            i=0
            n=len(chrom)
            INF=float(inf)
            if n>2:
                while i<n:
                    j=i+1
                    while j< len(chrom):
                        if chrom[i]==chrom[j]:
                            chrom[i]=INF
                        j+=1
                    i+=1
                i=0
            chrom2=[z for z in chrom if z != INF]
            mpool[ch]=np.zeros(len(chrom2))
            chrom2=np.array(chrom2)
            mpool[ch]=chrom2
            ch+=1
        return (mpool)
    
    def intra_clust(val, partition, img, band):
        
        red=rasterio.open(img).read(1)
        green=rasterio.open(img).read(2)
        blue=rasterio.open(img).read(3)
        if band ==4:
            nir=rasterio.open(img).read(4)
        #red=red[400:,400:]
        #green=green[400:,400:]
        #blue=blue[400:,400:]
        #if band ==4:
         #   nir=nir[400:,400:]
        peudis=[]
        cluster=np.where(partition == val) 
        print(val)
        #print('cluster:',cluster)
        ncl=len(cluster[1])
        nc=0
        nc1=0
        peudis1=[]
        #dup=0
        while nc <ncl:
            #print('nc=',nc)
            nc1+=nc
            while nc1 <ncl:
                #print('nc1=',nc1)
                #print('ncL=',ncl)
                
                if nc1 != nc:
                    p1=cluster[0][nc]
                    p2=cluster[1][nc]
                    p3=cluster[0][nc1]
                    p4=cluster[1][nc1]
                    if red[p1][p2]>red[p3][p4]:
                        a1=red[p1][p2]-red[p3][p4]
                    else:
                        a1=red[p3][p4]-red[p1][p2]
                    if green[p1][p2]>green[p3][p4]:
                        a2=green[p1][p2]-green[p3][p4]
                    else:
                        a2=green[p3][p4]-green[p1][p2]
                    if blue[p1][p2]>blue[p3][p4]:
                        a3=blue[p1][p2]-blue[p3][p4]
                    else:
                        a3=blue[p3][p4]-blue[p1][p2]
                    if band ==4:
                        if nir[p1][p2]>nir[p3][p4]:
                            a4=nir[p1][p2]-nir[p3][p4]
                        else:
                            a4=nir[p3][p4]-nir[p1][p2]

                        temp=np.sqrt((a1**2)+(a2**2)+(a3**2)+(a4**2))
                    else:
                        temp=np.sqrt((a1**2)+(a2**2)+(a3**2))
                    peudis1.append(temp)
                #else:
                    #dup=1
                nc1+=1
            if len(peudis1)!=0:
                edist=sum(peudis1)/len(peudis1)
                peudis.append(edist)
            nc+=1
        tintra=max(peudis)
        #print('Intra:',tintra)
        #print('Intra=',tintra)
        return tintra
            
    def inter_clust(chrom, p1,p2,p3,p4, partition, img, band):
        
        red=rasterio.open(img).read(1)
        green=rasterio.open(img).read(2)
        blue=rasterio.open(img).read(3)
        if band ==4:
            nir=rasterio.open(img).read(4)
        #red=red[400:,400:]
        #green=green[400:,400:]
        #blue=blue[400:,400:]
        #if band ==4:
            #nir=nir[400:,400:]
        #print('Inter Cluster')
        #cluster=np.where(partition == val) 
        if red[p1][p2]>red[p3][p4]:
            a1=red[p1][p2]-red[p3][p4]
        else:
            a1=red[p3][p4]-red[p1][p2]
        if green[p1][p2]>green[p3][p4]:
            a2=green[p1][p2]-green[p3][p4]
        else:
            a2=green[p3][p4]-green[p1][p2]
        if blue[p1][p2]>blue[p3][p4]:
            a3=blue[p1][p2]-blue[p3][p4]
        else:
            a3=blue[p3][p4]-blue[p1][p2]
        if band ==4:
            if nir[p1][p2]>nir[p3][p4]:
                a4=nir[p1][p2]-nir[p3][p4]
            else:
                a4=nir[p3][p4]-nir[p1][p2]
            tinter=math.sqrt((a1**2)+(a2**2)+(a3**2)+(a4**2))
        else: 
            tinter=math.sqrt((a1**2)+(a2**2)+(a3**2))
        #print('Inter=',tinter)
        return tinter
    
    
    #def fuzzy_membership(chrom, p1,p2,p3,p4, partition, img, band):
     #   red=rasterio.open(img).read(1)
      #  green=rasterio.open(img).read(2)
       # blue=rasterio.open(img).read(3)
        #if band ==4:
         #   nir=rasterio.open(img).read(4)
        #red=red[400:,400:]
        #green=green[400:,400:]
        #blue=blue[400:,400:]
        #if band ==4:
         #   nir=nir[400:,400:]
        #peudis=[]
        #cluster=np.where(partition == val)        
        