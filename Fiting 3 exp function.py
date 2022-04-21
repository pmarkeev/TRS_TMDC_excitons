# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 11:17:25 2022

@author: MarkeevP
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import glob
from sklearn import preprocessing

#READ DATA

base_folder = 'Traces graphs/MoS2 - 3/' 
folder= base_folder + '*.csv'
folder_save = base_folder + 'results/'
data_names = glob.glob(folder)
# Choose files from the folder w(with pent) n(no pentacene) TMD, p(just pentacene, no TMD)
w = 2
n = 0
p = 4

print('Sanity check!: with=  ' + data_names[w])
print('without=  ' + data_names[n])
print('pent=  ' + data_names[p])

#!!! Check names manually! 

#%% Assign the names (name1 - with pent, name2 - without, name3 - pentacene, BUT ALWAYS DOUBLE CHECK THE NAMING IN OTHER PROGRAM )

title= 'Trace 605-630 nm'
name1 = 'MoS\u2082-Pentacene'
name2 = 'MoS\u2082'
name3 = 'Pentacene'


#%% Read data

data=[]

for i in range(len(data_names)):
    data.append(i)
    #Read data and transpose from horizontal to vertical
    data[i]=pd.read_csv(data_names[i], index_col=0, header=None).T
    data[i].columns = ['0','1']
    
    #make all data positive time
    if data[i]['1'][6]<0:
        data[i]['1']=data[i]['1']*-1
    
    # #subset only positive time points (could be changes for better fit +-1)
    # data[i]=data[i].iloc[5:,:]
    
    print(str(i) + "---" + data_names[i])

#data[w] = data[w].iloc[2:,:]
#data[n] = data[n].iloc[2:,:]
data[p] = data[p].iloc[5:,:]

#%% NORMALIZE DATA
traceNorm=preprocessing.normalize([data[w].iloc[:,1], data[n].iloc[:,1]])
data[w]=pd.DataFrame({'0':data[w].iloc[:,0].tolist(), '1':traceNorm[0,:]})
data[n]=pd.DataFrame({'0':data[n].iloc[:,0].tolist(), '1':traceNorm[1,:]})
# data[p]=pd.DataFrame({'time':data[p].iloc[:,0].tolist(), 'signal':traceNorm[2,:]})
data[p]['1']=data[p]['1']*100


#%% Plot all data before fitting
x1 = data[n].iloc[:,0]
x2 = data[w].iloc[:,0]

plt.figure(dpi=600)
plt.scatter(x1,data[w].iloc[:,1], label = name1, s=4, color='red')
plt.scatter(x2,data[n].iloc[:,1], label = name2, s=4, color='blue')
plt.scatter(data[p].iloc[:,0], data[p].iloc[:,1], label = name3, s=4, color='green')
plt.xscale('log')
#plt.xlim(-1.5,0.5)
#plt.title(title)
plt.xlabel('Pump-Probe Delay (ps)')
plt.ylabel(u'ΔReflectivity')
#plt.legend(fontsize=9, scatterpoints=3)

# !!! save here without fitting line
    #plt.savefig(folder_save + title + ' w='+str(w) + ' n='+str(n) + ' p='+str(p) + '.png', transparent=True, dpi=300)


#%% fitting 3 exp on no pent

# objective function 3 down (for traces from TAS data)
#name_add = ''
def objective_n(t, a1, a3,   t1, t3):
    return a1*np.exp((-t/t1))+a3*np.exp(-t/t3)
# 6 unit guess
unit_guess=[0.5, 0.5, 0.6, 1000]
bound=([0, 0, 0, 0], [1, 1, 50, 2000])
index = ['a1', 'a3',  't1', 't3']

#fit curve
param1, cov1 = curve_fit(objective_n, x1[:], data[n].iloc[:,1], p0=unit_guess, bounds=bound, maxfev=2000)

# calculate the output for the range
y_line1 = objective_n(x1[:], *param1)

# create a line plot for the mapping function
plt.plot(x1[:], y_line1, '--', color='blue', label='Fit of '+ name2, alpha=0.4)

#put a values of a's in pecentage form
param1_new = np.empty(2)
for i in range(2):
    param1_new[i] = param1[i]/param1[0:2].sum()*100
param1[0:2]=list(param1_new)
    
results_n=pd.DataFrame({str(name2):param1, 'error':np.sqrt(np.diag(cov1))}, index = index)
results_n = results_n.apply(lambda x: round(x, 3))
print(results_n)

#%% fitting 4 exp on with pent keeping t1 and t3 the same

#t1 = param1[-3]
#t4 = param1[-1]

def objective_w(t, a1, a3,   t1, t3):
    return a1*np.exp((-t/t1))+a3*np.exp(-t/t3)
# 6 unit guess
unit_guess=[0.5, 0.5, 0.6, 1000]
bound=([0, 0, 0, 0], [1, 1, 50, 2000])
index = ['a1', 'a3',  't1', 't3']

#fit curve
param2, cov2 = curve_fit(objective_w, x2[:], data[w].iloc[:,1], p0=unit_guess, bounds=bound, maxfev=2000)

# calculate the output for the range
y_line2 = objective_w(x2[:], *param2)

# create a line plot for the mapping function
plt.plot(x2[:], y_line2, '--', color='red', label='Fit of '+ name1, alpha=0.4)
plt.legend(fontsize=9, scatterpoints=3, frameon=False)

#put a values of a's in pecentage form
param2_new = np.empty(2)
for i in range(2):
    param2_new[i] = param2[i]/param2[0:2].sum()*100
param2[0:2]=list(param2_new)

results_w=pd.DataFrame({str(name1):param2, 'error':np.sqrt(np.diag(cov2))}, index = index)
results_w = results_w.apply(lambda x: round(x, 3))
print(results_w)

# # #SAVING GRAPH
#plt.savefig(folder_save + title + ' w='+str(w) + ' n='+str(n) + ' p='+str(p) + name_add + '.png', transparent=True, dpi=300)
# #SAVING TEXT 
#results.to_csv(folder_save + 'results of fitting ' + title + name_add + ' w='+ str(w) + ' n='+ str(n) + ' p='+ str(p) +'.csv')
# pd.DataFrame([x1, y_line1]).to_csv(folder_save + 'fit line ' + name_add + name1 + ' w='+str(w) +'.csv')
# pd.DataFrame([x2, y_line2]).to_csv(folder_save + 'fit line ' + name_add + name2 + ' n='+str(n) +'.csv')

plt.show()
