from GA_pack import GA_pack as ga1
import pandas as pd
from NSGA2 import nsga
import numpy as np


nc=int(input("Enter the number of Chromosomes: "))
ulim=int(input("The upper limit of the size of a different length chromosome: "))
llim=int(input("The lower limit of the size of a different length chromosome: "))
gen=int(input("Enter the number of generation: "))

df = pd.read_excel('Portfolio_Clustering_Final.xlsx', sheet_name=['Sheet1'])
df1=list(df.items())
df2 = [item for t in df1 for item in t]
df3=df2[1]
df3=df3.values.tolist()
row=len(df3[0])
col=len(df3)

#i=0
#while i<col:
#    j=1
#    while j<row:
 #       df4.append((df3[i][j]+0.5))
  #      j+=1
   # i+=1
#print(df)
date_wise_cluster=[]
cluster_gen=[]   
c_round=0
while c_round<col:
    print('NSGA2 execution starts. Round ',(c_round+1))
    df4=[]
    newpool=[]
    i=1
    while i<row: 
        df4.append((df3[c_round][i]+0.5))
        i+=1
    newpool=nsga.nsga2(nc,ulim,llim,gen,df4)
    fitpool=ga1.fitcal(newpool,nc,df4)
    front,pop_front=ga1.nondomsort(nc,fitpool)
    k=0
    temp=front[0]
    interc=[]
    while k< len(temp):
        pos=temp[k]
        interc.append(fitpool[1][pos])
        k+=1
    opt=interc.index(max(interc))
    soln=newpool[front[0][opt]]
    date_wise_cluster.append(soln)
    ls=0
    while ls<len(soln):
        cluster_gen.append(soln[ls])
        ls+=1
    
    c_round+=1
                        
print('Final clustering starts...')
newpool=nsga.nsga2(nc,ulim,llim,gen,cluster_gen)                       
fitpool=ga1.fitcal(newpool,nc,cluster_gen)
front,pop_front=ga1.nondomsort(nc,fitpool)
ga1.plot_op(front[0], fitpool)
k=0
temp=front[0]
interc=[]
while k< len(temp):
    pos=temp[k]
    interc.append(fitpool[1][pos])
    k+=1
opt=interc.index(max(interc))
soln=newpool[front[0][opt]]     
i=0
df4=[]
while i<col:
    j=1
    while j<row:
        df4.append((df3[i][j]+0.5))
        j+=1
    i+=1
#print(df)                
partition_matrix=ga1.part_matrix(soln,1,df4)
part_matrix=np.zeros((col,row),np.uint8)
part_matrix=part_matrix.tolist()
k=0
i=0
while i<col:
    j=0
    #part_matrix[i][0]=df3[i][0]
    while j<row-1:
        part_matrix[i][j]=partition_matrix[k]
        k+=1
        j+=1
    i+=1


i=0
portfolio=[]
portfolio.append(soln-0.5)
while i<len(front[0]):
    portfolio.append(newpool[front[0][i]]-0.5)
    i+=1
z=0
writer = pd.ExcelWriter('partition_matrix.xlsx', engine = 'xlsxwriter')
while z<len(portfolio):
    soln=portfolio[z]
    partition_matrix=ga1.part_matrix(soln,1,df4)
    part_matrix=np.zeros((col,row),np.uint8)
    part_matrix=part_matrix.tolist()
    k=0
    i=0
    while i<col:
        j=0
        #part_matrix[i][0]=df3[i][0]
        while j<row-1:
            part_matrix[i][j]=partition_matrix[k]
            k+=1
            j+=1
        i+=1
    df2 = pd.DataFrame(part_matrix)
    temp=str(z)
    df2.to_excel(writer, sheet_name = 'Partition_mat for Portfolio'+temp)
    z+=1
writer.save()
writer.close() 
date_wise_cluster2=np.zeros((col,row),np.uint8)
date_wise_cluster2=date_wise_cluster2.tolist()
i=0
while i<col:
    j=1
    date_wise_cluster2[i][0]=df3[i][0]
    while j<len(date_wise_cluster[i])+1:
        date_wise_cluster2[i][j]=(date_wise_cluster[i][j-1]-0.5)
        j+=1
    i+=1

#df2 = pd.DataFrame(part_matrix)
df5 = pd.DataFrame(date_wise_cluster2)
df1 = pd.DataFrame(portfolio)
writer = pd.ExcelWriter('Portfolio_details.xlsx', engine = 'xlsxwriter')
df1.to_excel(writer, sheet_name = 'Portfolios')
#df2.to_excel(writer, sheet_name = 'Partition Matrix')
df5.to_excel(writer, sheet_name = 'Date wise cluster')
writer.save()
writer.close()

 