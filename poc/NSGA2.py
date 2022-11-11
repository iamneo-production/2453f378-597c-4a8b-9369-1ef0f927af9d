#  from osgeo import gdal, osr
#from SSIM_PIL import compare_ssim

import matplotlib.pyplot as plot
from GA_pack import GA_pack as ga1
from GA_Main import GA_Main as ga2
from GA_help import GA_help as ghelp

class nsga:

    def nsga2(nc,ulim,llim,gen,df4):    
        mpool=[]
        newpool=[]

        nobj=2
        fit_type=[0,1]

        mpool=ga1.encode(nc, ulim, llim, df4)
        mpool=ghelp.duplicate_elim(mpool, nc)
        #print('Encoding is done')
        i=0
        newpool=[]
        front_best=[]
        newpop=[]
        newpop_front=[]
        newcrowd_d=[] 
        newpool=[]
        front=[]
        Sp=[]
        crowd_dist=[]
        fitpool=[]
        fit_repository=[]

        while i<gen:
            if i!=0:
                fitpool,newpop, newcrowd_d, mpool=ga2.combine_offspring(nc, df4,nobj,fit_type,mpool,newpool,fitpool)
                #fitpool=ga1.fitcal(newpool,nc,band,ln_row,ln_col,image)
                #print('fit cal is done')
                #a=len(mpool)
                #print('length=',a)
                front,pop_front=ga1.nondomsort(nc,fitpool)
                #ga1.plot_op(front[0], fitpool)
                #print('sorting is done')
                #ga1.plot_op(front[0], fitpool)
                front_best.append(front[0])
                #front,pop_front=ga1.nondomsort(nc,fitpool)
                #fit_repository.append(fitpool)
                newpool=ga2.CrowdedbinTournamentSelec(mpool,pop_front,newcrowd_d,nc)
                #print('selection is done')
                #i+=1
            else:
                fitpool=ga1.fitcal(mpool,nc,df4)
                #fit_repository.append(fitpool)
                front,pop_front=ga1.nondomsort(nc,fitpool)
                newpool=ga2.binTournamentSelec(mpool,pop_front,nc)
                #front,pop_front=ga1.nondomsort2(nc,fitpool)
                #ga1.plot_op(front[0], fitpool)
                #print('sorting is done')
                front_best.append(front[0])
                #mpool=ga2.binTournamentSelec(mpool,newpop_front,newcrowd_d,newpool,nc)
                #print('selection is done')
            newpool=ga2.crossover(newpool)
            newpool=ghelp.duplicate_elim(newpool, nc)
            #print('Crossover is done')
            newpool=ga2.mutation(newpool,len(df4),df4)
            #newpool=ghelp.duplicate_elim(newpool, nc)
            #print('Mutation is done')
            #ga1.plot_op(front_best)
            print('End of generation:',(i+1))
            i+=1
        return newpool
    def nsga2_1(nc,ulim,llim,gen,df4):    
        mpool=[]
        newpool=[]

        nobj=2
        fit_type=[0,1]

        mpool=ga1.encode(nc, ulim, llim, df4)
        mpool=ghelp.duplicate_elim(mpool, nc)
        #print('Encoding is done')
        i=0
        newpool=[]
        front_best=[]
        newpop=[]
        newpop_front=[]
        newcrowd_d=[] 
        newpool=[]
        front=[]
        Sp=[]
        crowd_dist=[]
        fitpool=[]
        fit_repository=[]

        while i<gen:
            if i!=0:
                fitpool,newpop, newcrowd_d, mpool=ga2.combine_offspring(nc, df4,nobj,fit_type,mpool,newpool,fitpool)
                #fitpool=ga1.fitcal(newpool,nc,band,ln_row,ln_col,image)
                #print('fit cal is done')
                #a=len(mpool)
                #print('length=',a)
                front,pop_front=ga1.nondomsort(nc,fitpool)
                #ga1.plot_op(front[0], fitpool)
                #print('sorting is done')
                ga1.plot_op(front[0], fitpool)
                front_best.append(front[0])
                #front,pop_front=ga1.nondomsort(nc,fitpool)
                #fit_repository.append(fitpool)
                newpool=ga2.CrowdedbinTournamentSelec(mpool,pop_front,newcrowd_d,nc)
                #print('selection is done')
                #i+=1
            else:
                fitpool=ga1.fitcal(mpool,nc,df4)
                #fit_repository.append(fitpool)
                front,pop_front=ga1.nondomsort(nc,fitpool)
                newpool=ga2.binTournamentSelec(mpool,pop_front,nc)
                #front,pop_front=ga1.nondomsort2(nc,fitpool)
                ga1.plot_op(front[0], fitpool)
                #print('sorting is done')
                front_best.append(front[0])
                #mpool=ga2.binTournamentSelec(mpool,newpop_front,newcrowd_d,newpool,nc)
                #print('selection is done')
            newpool=ga2.crossover(newpool)
            newpool=ghelp.duplicate_elim(newpool, nc)
            #print('Crossover is done')
            newpool=ga2.mutation(newpool,len(df4),df4)
            #newpool=ghelp.duplicate_elim(newpool, nc)
            #print('Mutation is done')
            #ga1.plot_op(front_best)
            print('End of generation:',(i+1))
            i+=1
        return newpool