# -*- coding: utf-8 -*-
"""
Created on Tue May 11 15:26:15 2021

@author: MarkeevP
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing
import glob

#!!! DOUBLE CHECK WHAT IS GETTING SAVED 

#%% READ DATA 

folder='WSe2 Au new (for article)/double_check'
data_names = glob.glob(folder+'/*.csv')
file_number=0
# for save destination use name to save at the same folder or save_to in a specific folder
#name=data_names[file_number][:-4]

# from what range get the trace
t1=675
t2=695
# min max on 3d graph of z axis
vmin = -0.001
vmax = 0.001

#save_name = data_names[file_number][-21:-6]
save_name = ' WSe2 ave '+ str(file_number) 

# #read data file excel
data1 = pd.read_excel('C:/Users/MarkeevP/Desktop/ELI-Alps/Transient absorption/python/Pentacene/TR_Results_Pentacene_only_1-36time_points_04May2021_P1.xlsx', sheet_name='Results_Sample_xxYY2020', index_col=2)
data1 = data1.drop(data1.columns[[0,1,2,-1,-2,-3,-4,-5,-6]], axis=1)
data1 = data1.drop(data1.index[[0,1,2,3]])
### drop nan last row in some cases
# data1 = data1.drop(data1.index[-1])

# # read csv with dots
# # data1=pd.read_csv(data_names[file_number], index_col=0)
# data1=pd.read_csv(
#     'C:/Users/MarkeevP/Desktop/ELI-Alps/Transient absorption/python/MoS2 double check again 18.08.2021/Pent ave.csv', index_col=0)
# # read csv with dots by name (not glob)

#data1 = pd.read_csv(folder+'/WSe2 pent.csv', index_col=0)
name = folder+'/Pent'

#subset data from 550-750 nm wavelength
data1=data1.iloc[:,863:1350]

# Make arrays for columns(wavelength) and rows(time) data
col=data1.columns.astype(float) #col is array of all column names
rows=data1.index.astype(float) #col is array of all column names
rows=[round(x, 2) for x in rows]


#!!! check the parameters file in the folder to set the same limits for subtraction!!!

#go thru col values (wavelength values) and find number  of column for needed wavelength
#change g2 g3 and g4 to set where subtraction will happen (where you know that data is 0 and have only noise)  
for i1 in range(len(col)):
    if int(col[i1])==565:
        g2=i1
    if int(col[i1])==610:
        g3=i1
    if int(col[i1])==715:
        g4=i1
    if int(col[i1])==690:
        gflip=i1

#%% DATA FLIP (there could be a random sign change in raw data)

def flip(data):
    dataRes0=pd.DataFrame(1, columns=data.columns, index=data.index)
    for i in range(data.shape[0]):
        if data.iloc[i,gflip]>0:
            dataRes0.iloc[i,:]=data.iloc[i,:]*-1
        else:
            dataRes0.iloc[i,:]=data.iloc[i,:]
    return dataRes0

datares0 = flip(data1)       

#find index of the interested Wavelength with lambda, not for loop. this index is used to plot a line on the center of exciton
col1=list(col)
takeClosest = min(col1,key=lambda x : abs(x-611))

#takeClosest = lambda num,collection : min(collection,key=lambda x : abs(x-num))
gline1 = list(col).index( min(col1,key=lambda x : abs(x-550)))
gline2 = list(col).index( min(col1,key=lambda x : abs(x-674)))

#Plot flipped data (vmin and vmax - z-axis limits)
plt.imshow(datares0, aspect='auto',  vmin = -0.004, vmax = 0.004, interpolation='gaussian')

#define appropriate axes to the plot
x = col.astype(int) # the grid to which your data corresponds
nx = x.shape[0] #length of x
no_labels = 11 # how many labels to see on axis x
step_x = int(nx / (no_labels - 1)) # step between consecutive labels
x_positions = np.arange(0,nx,step_x) # pixel count at label position
x_labels = x[::step_x] # labels you want to see
plt.xticks(x_positions, x_labels) #set defined x positions and labels to the graph 

# .apply(lambda x: '%.2e' % x) - apply format of rounding numbers for 2 decimals

# for y axis use the same procedure as for x
y = datares0.index
ny = y.shape[0]
step_y =  int(ny/10)
y_positions = np.arange(0, ny, step_y)
y_labels = y[::step_y] 
y_labels=[round(x, 2) for x in y_labels]
plt.yticks(y_positions, y_labels)

plt.colorbar(label=u'Δ mOD', format='%.1e') #add color bar
plt.xlabel('Probe Wavelength (nm)')
plt.ylabel('Pump-Probe Delay (ps)') 
plt.tight_layout()
# plt.savefig(name + '_original_flipped.png', dpi=150) #save 

# # add vertical lines that show where trace was taken
# plt.axvline(x=gline1, ymin=0, ymax=1, c='r', ls='--') #plot a line of the center of exciton
# plt.axvline(x=gline2, ymin=0, ymax=1, c='r', ls='--') #plot a line of the center of exciton

plt.show()

#%% FLATTEN THE DATA

#process the data: by each line subtract the value at certain posicion (g2, g3) and along certain fitted slopes
def process_data(data):
    
    #create empty (all ones) dataframe for processed data
    dataRes=pd.DataFrame(1, columns=data.columns, index=data.index)
    
    xfit1=np.array([g2,g3]).reshape(-1,1)
    xfit2=np.array([g3,g4]).reshape(-1,1)  
                
    for i in range(data.shape[0]):
        
        # 2 linear fits along the pekas (could be only one used)
        yfit1=np.array([data.iloc[i,g2],data.iloc[i,g3]]).reshape(-1,1)
        fit1 = LinearRegression().fit(xfit1, yfit1) 
        a1 = fit1.coef_
        b1 = fit1.intercept_
        
        yfit2=np.array([data.iloc[i,g3],data.iloc[i,g4]]).reshape(-1,1)
        fit2 = LinearRegression().fit(xfit2, yfit2) 
        a2 = fit2.coef_
        b2 = fit2.intercept_
        
        # loop to subtract fitted lines from the data at each point
        for j in range (len(col[g2:g3])):
            dataRes.iloc[i,g2+j] = data.iloc[i,g2+j]-float(a1*(g2+j)+b1)
        
        for j in range (len(col[g3:g4])):
            dataRes.iloc[i,g3+j] = data.iloc[i,g3+j]-float(a2*(g3+j)+b2)
        
        # subtract same value at horizontal (stable) parts of the graphs 
        dataRes.iloc[i,:g2]=data.iloc[i,:g2]-data.iloc[i,g2]
        dataRes.iloc[i,g4:]=data.iloc[i,g4:]-data.iloc[i,g4]
                    
    return dataRes

#datares1 = process_data(data1)
datares1 = process_data(datares0)

#%% PLOT EVERYTHING


#Plot original data (vmin and vmax - z-axis limits)
plt.imshow(data1, aspect='auto',  vmin = -0.004, vmax = 0.004, interpolation='gaussian', cmap='coolwarm')

#define appropriate axes to the plot
x = col.astype(int) # the grid to which your data corresponds
nx = x.shape[0] #length of x
no_labels = 11 # how many labels to see on axis x
step_x = int(nx / (no_labels - 1)) # step between consecutive labels
x_positions = np.arange(0,nx,step_x) # pixel count at label position
x_labels = x[::step_x] # labels you want to see
plt.xticks(x_positions, x_labels) #set defined x positions and labels to the graph 

# for y axis use the same procedure as for x
y = data1.index
ny = y.shape[0]
step_y =  int(ny/10)
y_positions = np.arange(0, ny, step_y)
y_labels = y[::step_y] 
y_labels=[round(x, 2) for x in y_labels]
plt.yticks(y_positions, y_labels)

plt.colorbar(label=u'Δ mOD', format='%.1e') #add color bar
plt.xlabel('Probe Wavelength (nm)')
plt.ylabel('Pump-Probe Delay (ps)') 
plt.tight_layout()
# plt.savefig(name +  '_original.png', dpi=150) #save 
plt.show()


#Plot processed 1 data (vmin and vmax - z-axis limits)
plt.figure(dpi=1200)
plt.imshow(datares1, aspect='auto',  vmin = -0.004, vmax = 0.004, 
           interpolation='gaussian', cmap='coolwarm')

#define appropriate axes to the plot
x = col.astype(int) # the grid to which your data corresponds
nx = x.shape[0] #length of x
no_labels = 11 # how many labels to see on axis x
step_x = int(nx / (no_labels - 1)) # step between consecutive labels
x_positions = np.arange(0,nx,step_x) # pixel count at label position
# x_labels = x[::step_x] # labels you want to see
plt.xticks(x_positions , np.arange(550,751,20)) # hard coded labels to the x axis
# plt.xticks(x_positions, x_labels) #set defined x positions and labels to the graph 

# for y axis use the same procedure as for x
y = data1.index
ny = y.shape[0]
step_y =  int(ny/10)
y_positions = np.arange(0, ny, step_y)
y_labels = y[::step_y] 
y_labels=[round(x, 2) for x in y_labels]
plt.yticks(y_positions, y_labels)

plt.colorbar(label=u'ΔReflectance', format='%.1e') #add color bar
plt.xlabel('Probe Wavelength (nm)')
plt.ylabel('Pump-Probe Delay (ps)') 
plt.tight_layout()

#add lines
# plt.axvline(x=gline1, ymin=0, ymax=1, c='r', ls='--', alpha=0.6) #plot a line of the center of exciton
# plt.axvline(x=gline2, ymin=0, ymax=1, c='r', ls='--', alpha=0.6) #plot a line of the center of exciton

# save processed image
#plt.savefig(name + save_name + '_param(' + str(int(col[g2]))+'-'+str(int(col[g3]))+'-' + str(int(col[g4]))+')_processed.png', dpi=150)
plt.show()


#%% To plot processed data with propped Y axis
datares2 = datares1.iloc[1: , :].apply(lambda x: x*1000)
plt.figure(dpi=1200)
plt.pcolormesh( datares2.columns.to_numpy(), datares2.index.to_numpy(),
               datares2.to_numpy(),  vmin = -2, vmax = 2, shading ='gouraud', cmap='coolwarm')
plt.yscale('log')
#plt.gca().invert_yaxis()

plt.ylim(1750, 0.1)

plt.xticks(x_positions , np.arange(550,751,20))
plt.colorbar(label=u'ΔReflectance', format='%.1f') #add color bar
plt.xlabel('Probe Wavelength (nm)')
plt.ylabel('Pump-Probe Delay (ps)') 
plt.tight_layout()
%% TRACES

# GET TRACE
# set wavelength at which you want to get traces
# t1=605                                 # Now is changed at the start of the file
# t2=630
g1=[]
#go thru col values (wavelength values) and find number  of column for needed wavelength
for i in range(len(col)):
    if col[i]>t1 and col[i]<t2:
        g1.append(i)
#function to get traces 
def get_trace(data):
    trace=[]
    for i in range(data.shape[0]):
        trace.append(data.iloc[i,g1].mean())
    trace=pd.DataFrame([data.index,trace])
    return trace
#get trace
trace1=get_trace(datares1)
trace1=trace1.iloc[:,5:]

# #!!! delete one element in specific data case !!! 
# trace1.drop(40, axis=1, inplace=True)

#plot 2 traces (with and without pentacene)
plt.figure(dpi=1200)
plt.scatter(trace1.iloc[0,:], trace1.iloc[1,:], marker='x', color='blue')
#plt.title('Trace '+str(t1)+'-'+str(t2)+' nm')
plt.xlabel('Pump-Probe Delay (ps)')
plt.ylabel(u'Δ mOD')
plt.xscale('log')

# #!!! manual apply y axis
# plt.yticks([-0.0016,-0.0014,-0.0012, -0.0010, -0.0008], [-1.6,-1.4,-1.2,-1.0,-0.8])


#plt.xlim(-1,6)
#plt.tight_layout()
plt.legend()
# SAVING plot

# plt.savefig(name + save_name + '_TRACE at '+str(t1)+'--'+str(t2)+'.png', dpi=300) 

plt.show()
#save trace data
#np.savetxt(name + save_name + '_TRACE at '+str(t1)+'--'+str(t2)+' .csv', trace1, delimiter=",")



%% Spectras

# GET SPECTRUM
#function to get spectrum at s (any time)
s=0.53
s_index = list(rows).index( min(rows,key=lambda x : abs(x-s)))
def get_spectrum(data):
    Spectrum=[]
    for i in range(data.shape[1]):
        Spectrum.append(data.iloc[s_index,i])
    Spectrum=pd.DataFrame([data.columns,Spectrum])
    return Spectrum

#get spectrum for resulted data sets
spectra=get_spectrum(datares1)
#col1=list(col)
cut = list(col).index( min(col,key=lambda x : abs(x-550))) #cut spectrum at 550
spectra=spectra.iloc[:,cut:]

plt.figure(dpi=1200)
plt.plot(spectra.iloc[0,:], spectra.iloc[1,:])

#!!! manual Y and X 
plt.xticks(x_positions , np.arange(550,751,20))
plt.yticks([-0.0035,-0.003,-0.0025, -0.002,-0.0015,-0.0010, -0.0005, -0.0000, 0.0005],
            [-3.5,-3.0, -2.5,-2.0,-1.5,-1.0,-0.5,0, 0.5])

plt.xlabel('Probe Wavelength (nm)')
plt.ylabel(u'Δ mOD')

# SAVING

# plt.savefig(name + save_name + '_SPECTRUM at '+str(s)+' delay.png', dpi=150)

plt.show()

# np.savetxt(name +  save_name + '_SPECTRUM at '+str(data1.index[s_index])+' delay.csv', spectra, delimiter=",", fmt='%s') 



