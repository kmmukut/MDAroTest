#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import networkx as nx


# In[2]:


df = pd.io.parsers.read_csv(
    filepath_or_buffer='points.csv',
    header=None,
    sep=',',
    )
x=np.array(df)
n=df.shape[0]
d=np.zeros(n)

for c in range(n):
    d[c]=np.sqrt(x[c,0]*x[c,0]+x[c,1]*x[c,1]+x[c,2]*x[c,2]) 


# In[3]:


main = np.hstack((d.reshape(n,1), x.reshape(n,3))).real
main_sorted_arg = np.argsort(main[:, 0])
main_sort = np.asarray(main[main_sorted_arg])
main_sort=np.asarray(main_sort)


# In[4]:


#bond distance
bond_distance = 2/ np.sqrt(3) + 0.001


# In[5]:


# marker = np.arange(len(main_sort))
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.scatter(
#     main_sort[:,1],
#     main_sort[:,2],
#     cmap='rainbow',
#     alpha=0.7,
#     edgecolors='c',
# )
# plt.xlim(-0.5,  14)
# plt.ylim(-0.5, 14)
# for x in marker:
#     plt.text(main_sort[x][1], main_sort[x][2], '({})'.format(x))


# In[6]:


def checkDuplicate(neddleCycle, cycles):
    neddleCycle = np.asarray(neddleCycle)
    cycles = np.asarray(cycles)
    
    for cycle in cycles:
        if len(cycle) != neddleCycle.shape[0]:
            continue
        if np.all(np.isin(neddleCycle, cycle)):
            return True
    return False


main_sort=np.asarray(main_sort)
graph = distance_matrix(main_sort[:, 1:], main_sort[:, 1:]) <= bond_distance
np.fill_diagonal(graph, False)


G = nx.DiGraph()
G = nx.from_numpy_matrix(graph, create_using=nx.DiGraph)
all_cycles = list(nx.simple_cycles(G))
cycles = []

for cycle in all_cycles:
    if len(cycle) > 2 and not checkDuplicate(cycle, cycles):
        cycles.append(cycle)
#cycles.append([21,22,23,44,55,66,77,88,99,100,1201])
# cycles


# In[7]:


boro_array = []
choto_array = []

def isSubSet(current, toCompare):
    
    return np.in1d(toCompare,current).all()

i = 0
for r in cycles:
    if(len(r)) < 8:
        choto_array.append(r)
        i += 1
        continue
    if i > 0:
        prev_item = cycles[i-1]
        
        if isSubSet(r, prev_item):
            #remove r from cycle or do some other thing
            #print("Subset paise agertay")
            #print(r)
            #print(prev_item)
            boro_array.append(r)
            i+=1
            continue
        
        
    if i < len(cycles)-1:
        next_item = cycles[i+1]
        
        if isSubSet(r, next_item):
            #remove r from cycle or do some other thing
            #print("Subset paise porertay")
            #print(r)
            #print(next_item)
            boro_array.append(r)
            i+=1
            continue
        
    
    choto_array.append(r)
    
    i += 1
    
# print("Boro array",boro_array)
# print("Choto array",choto_array)


a = {}
b = {}


for r in choto_array:
    if len(r) in a:
        a[len(r)] = a[len(r)] + 1
    else:
        a[len(r)] = 1

for l in boro_array:
    if len(l) in b:
        b[len(l)] = b[len(l)] + 1
    else:
        b[len(l)] = 1

    


# In[8]:


ARnumber=np.unique(np.hstack(cycles)).shape[0]
ALnumber=df.shape[0]-ARnumber
PercentageAromatic=ARnumber/(ARnumber+ALnumber)
PercentageAlephatic=ALnumber/(ARnumber+ALnumber)
print("Existing Rings\t",a)
print("Molecular Rings\t",b)
print("Number of Aromatic Carbon Atom\t:",ARnumber)
print("Number of Alephatic Carbon Atom\t:",ALnumber)
print("Percentage of Aromatic components\t:",PercentageAromatic)
print("Percentage of Alephatic components\t:",PercentageAlephatic)


# In[9]:


# graph


# In[10]:


# dp = pd.DataFrame(np.asarray(graph, dtype=np.uint8))#[11:22, 11:22])
# dp


# In[11]:


# np.asarray(graph, dtype=np.uint8)


# In[ ]:




