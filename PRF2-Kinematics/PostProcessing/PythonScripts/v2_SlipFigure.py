#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 14:58:21 2019

@author: thomas
"""

#Slip Figure
#We will be showing the velocity of each moving part as well as the body.
#There will be regions where the velocity of all 3 move in the same direction
#This is what we are calling slip
#There will be slip regions before and after expansion and compression
#The goal is to that the oscillation (in the lab frame) is not as simple as it seems
#We want to highlight the regions of slip for both SSL and LSL swimmers
#There is a difference and we want to make sure that difference is apparent in this figure
#First, store posData
#Subset 3 periods worth of data starting at SSTIME
#Calculate velocities of SS, LS, and CM
#Plot vS, vL, and vCM for both SSL and LSL

#MODULES
import os,sys
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)
from matplotlib.ticker import FormatStrFormatter
from scipy.signal import savgol_filter

#CONSTANTS
cwd_PYTHON = os.getcwd()
PERIOD = 0.1
FREQ = 10.0
DENS = 2.0
MAXAMP = 0.3
RSMALL = 0.15
RLARGE = 0.3
RSL_DEFAULT = 0.75
dt = 1.0e-3

FIGNUM = 0

mpl.rcParams['axes.linewidth'] = 1.0 #set the value globally

print('BODYLENGTH = ',RSMALL + RLARGE + RSL_DEFAULT*1.3)

def StoreData(cwd_DATA):
    #Reset position data every Re
    allData = []
    #Load position data
    allData = pd.read_csv(cwd_DATA+'/pd.txt',delimiter = ' ')
    return allData

def CalcVel(data):
    #Calculate the velocity of each individual sphere
    #Caculate the components and calculate the magnitude
    #Variables to use
    nTime = int(len(data.time))
    #Velocity in the x-direction and y-direction
    data['vL'], data['vS'], data['vCM'] = 0.0, 0.0, 0.0
    for idxTime in range(nTime):
        if((idxTime >= 1) and (idxTime < nTime-1)):
            #Finite difference: midpoint method
            d2t = data.time[idxTime+1] - data.time[idxTime-1]
            #Large Sphere
            data.loc[idxTime,'vL'] = (data.loc[idxTime+1,'Ly'] - data.loc[idxTime-1,'Ly'])/d2t
            #Small Sphere
            data.loc[idxTime,'vS'] = (data.loc[idxTime+1,'Sy'] - data.loc[idxTime-1,'Sy'])/d2t
            #Center of Mass
            data.loc[idxTime,'vCM'] = (data.loc[idxTime+1,'yCM'] - data.loc[idxTime-1,'yCM'])/d2t
    #Velocity at initial time == vel(1)
    data.loc[0,'vL'] = data.loc[1,'vL']
    data.loc[0,'vS'] = data.loc[1,'vS']
    data.loc[0,'vCM'] = data.loc[1,'vCM']
    #Velocity at final time == vel(ntime-1)
    data.loc[nTime-1,'vL'] = data.loc[nTime-2,'vL']
    data.loc[nTime-1,'vS'] = data.loc[nTime-2,'vS']
    data.loc[nTime-1,'vCM'] = data.loc[nTime-2,'vCM']
    #Smooth the Velocity Data!
    data.vL = SmoothVelocity(data,'vL')/(RSMALL*FREQ)
    data.vS = SmoothVelocity(data,'vS')/(RSMALL*FREQ)
    data.vCM = SmoothVelocity(data,'vCM')/(RSMALL*FREQ)
    
    return data

def SmoothVelocity(data,varName):
    #After plotting the raw data for forces, the term mdu/dt was very noisy
    #We will apply a smoothing algorithm to all velocities in the y-direction
        
    #Use Savitsky-Golay filter from scipy
    yval = data[varName].tolist()
    varValue = savgol_filter(yval, 51, 3)
    
    return varValue

def PlotVel(data,name):
    #Plot the position to make sure the appending was done properly
    fig = plt.figure(num=1,figsize=(6,4),dpi=200)
    ax = fig.add_subplot(111)
    csfont = {'fontname':'Times New Roman'}
    #ax.set_title(r'Center of Mass',fontsize=20)
    ax.set_xlabel(r'time $\tau$',fontsize=14,**csfont)
    ax.set_ylabel(r'$\hat v_{CM}$, $\hat v_{R}$',fontsize=14,**csfont)
    #Horizontal Lines
    ax.plot([-0.25,1.25],[0.0,0.0] ,color='k',lw=2)
    ax.plot(data.tau,data.vL,color=(255/255,127/255,0),lw=2)
    ax.plot(data.tau,data.vCM,color=(152/255,78/255,163/255),lw=2)
    ax2 = ax.twinx()
    ax2.set_ylabel(r'$\hat v_{r}$',fontsize=14,**csfont,color=(55/255,126/255,184/255))
    ax2.plot(data.tau,data.vS,color=(55/255,126/255,184/255),lw=2)
    
    #Set axes limits
    minY, maxY = int(np.amin(data.vL))-1, int(np.amax(data.vL))+1
    ax.axis([-0.1,1.1,-3.5,3.5])
    ax2.axis([-0.1,1.1,-7.5,7.5])
    
    #Plot Vertical Lines
    PlotVertLines(ax,0.0,-3.5,3.5)
    
    #Change Axes Parameters
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    ax.set_axisbelow(False)
    ax.set_xticks(np.arange(0.0,1.1,step=0.25))
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(1.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.75))
    #Ax2
    ax2.tick_params(which='major',axis='y',direction='in',length=6,width=1,labelcolor='tab:blue')
    ax2.tick_params(which='minor',axis='y',direction='in',length=4,width=0.75,labelcolor='tab:blue')
    ax2.yaxis.set_major_locator(MultipleLocator(3.0))
    ax2.yaxis.set_minor_locator(MultipleLocator(1.5))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'/../Version2/PaperFigures/Images/v2_Slip_'+str(name)+'.svg')
    fig.clf()
    plt.close()
    
    #sys.exit(0)
    return

def PlotSlipRegions(data,name):
    #There will be 2 slip regions
    #One from tau = [-0.05,0.05]
    #Another from [0.45,0.55]
    #Tau = [0.95, 1.05]
    fig = plt.figure(num=2,figsize=(2.5,2),dpi=200)
    ax = fig.add_subplot(111)
    csfont = {'fontname':'Times New Roman'}
    #ax.set_title(r'Center of Mass',fontsize=20)
    #ax.set_xlabel(r'time $\tau$',fontsize=14,**csfont)
    #ax.set_ylabel(r'$v_{CM}$ and $v_{R}$ (m/s)',fontsize=14,**csfont)
    #Horizontal Lines
    ax.plot([-0.25,1.25],[0.0,0.0] ,color='k',lw=2)
    ax.plot(data.tau,data.vL,color='tab:red',lw=2)
    ax.plot(data.tau,data.vCM,color='tab:green',lw=2)
    ax2 = ax.twinx()
    #ax2.set_ylabel(r'$v_{r}$ (m/s)',fontsize=14,**csfont,color='tab:blue')
    ax2.plot(data.tau,data.vS,color='tab:blue',lw=2)
    
    #Set axes limits
    minY, maxY = int(np.amin(data.vL))-1, int(np.amax(data.vL))+1
    ax.axis([0.95,1.05,-1.1,1.1])
    ax2.axis([0.95,1.05,-2.6,2.6])
    
    #Plot Vertical Lines
    PlotVertLines(ax,0.0,-5,5)
    
    #Change Axes Parameters
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    ax.set_axisbelow(False)
    ax.set_xticks(np.arange(0.95,1.06,step=0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.01))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))
    #Ax2
    ax2.tick_params(which='major',axis='y',direction='in',length=6,width=1,labelcolor='tab:blue')
    ax2.tick_params(which='minor',axis='y',direction='in',length=4,width=0.75,labelcolor='tab:blue')
    ax2.yaxis.set_major_locator(MultipleLocator(1.0))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.5))

    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'/../PaperFigures/Results/Slip/Slip1_'+str(name)+'.png')
    #Do the same but for #2
    #Tau = [0.45,0.55]
    #Set axes limits
    ax.axis([0.45,0.55,-1.1,1.1])
    ax2.axis([0.45,0.55,-2.6,2.6])
    #Change Axes Parameters
    ax.set_xticks(np.arange(0.45,0.56,step=0.05))
    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'/../PaperFigures/Results/Slip/Slip2_'+str(name)+'.png')
    
    fig.clf()
    plt.close()
    
    #sys.exit(0)
    return

def PlotVertLines(ax,xMin,yMin,yMax):
    #Vertical Lines Indicating Expansion/Compression
    ax.plot([xMin,xMin],[yMin,yMax],color='k',ls='--',zorder=0)
    ax.plot([xMin+0.25,xMin+0.25],[yMin,yMax],color='gray',ls='--',zorder=0)
    ax.plot([xMin+0.5,xMin+0.5],[yMin,yMax],color='k',ls='--',zorder=0)
    ax.plot([xMin+0.75,xMin+0.75],[yMin,yMax],color='gray',ls='--',zorder=0)
    ax.plot([xMin+1.0,xMin+1.0],[yMin,yMax],color='k',ls='--',zorder=0)
    return

if __name__ == '__main__':
    #This is where the main part of the code will be conducted
    #Functions will be called here
    
    '''SSL and LSL FIGURES'''
    nameList = ['Re2.5','Re70.0']
    #nameList = ['Re70.0']
    ssList = [3.8,6.6]
    for idx in range(len(nameList)):
        name = nameList[idx]
        SSTIME = ssList[idx]
        #Obtain SSL/LSL Data directory
        cwd_DATA = cwd_PYTHON+'/../HydroForces/ForceData/'+name+'/'
        #Store Data in database
        sslData = StoreData(cwd_DATA)
        #keep data of importance
        dataDict = {'time':sslData.time,'Ly':sslData.yUp,'Sy':sslData.yLow}
        data = pd.DataFrame(data=dataDict)
        #Subset data to 3 periods after Steady State
        data = data[data.time >= SSTIME].copy()
        data = data[data.time < SSTIME+3.0*PERIOD].copy()
        data = data.reset_index(drop=True)
        print('Min Time = ',np.amin(data.time))
        print('Max Time = ',np.amax(data.time))
        #Create new important variables
        data['tau'] = data.time/PERIOD
        data.tau = data.tau - 10.0*(SSTIME+0.075)
        data['yCM'] = 0.8*data.Ly + 0.2*data.Sy
        #Calculate velocities of CM and Small Sphere
        data = CalcVel(data)
        #vL, vS, and vCM calculated
        #Subset data to 1.5 periods (tau = 0.0 is when time = 1.075s)
        plotData = data[500:2000].copy()
        plotData = plotData.reset_index(drop=True)
        print(str(name)+' Data:')
        print(len(plotData.time))
        print('Answer should be 1500')
        print('Delta Time = ',plotData.loc[len(plotData.time)-1,'time'] - plotData.loc[0,'time'])
        print('tau start = ',data.loc[0,'tau'])
        #Plot Velocities
        PlotVel(plotData,name)
        #Plot Slip
        #PlotSlipRegions(plotData,name+'_2')
    