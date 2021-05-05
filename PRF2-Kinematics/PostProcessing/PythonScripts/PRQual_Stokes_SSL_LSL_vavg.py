#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:00:41 2020

@author: thomas
"""

#Power/Recovery: Qualitative Figure
#We will be showing the power and recovery locations for swimmers of different Re
#We can either show it using v_CM or by showing y_CM. What will be used is tbd
#Power/Recovery will be determined by the small sphere's velocity, v_r
#SSL: When v_r < 0 = Power; else Recovery
#LSL: When v_r > 0 = Power; else Recovery
#We will plot X rows by 1 column, where X is # of Re values
#Re will increase as we go down rows and Stokes will be included
#X-axis: tau, Y-axis: y_CM or v_CM
#Tau = [-0.25,1.25]
#We will need to append data ranging from [-0.25,0.0) and (1.0,1.25]
#Data will be sorted by tau, and the indices will be reset
#We plot v_CM vs Tau
#Next, we need to see when power/recovery starts/ends.
#This will be done using a root finder on v_r
#v_r will be calculated using pos data from posRSL#.#Shifted.csv
#RSL will initially be 1.3
#We will change the color of the curve based on if it is P or R
#We will also also vertical dashed lines where P and R start/end

#MODULES
import os,sys
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)
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

ReData = [None]*20
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
    data.vL = SmoothVelocity(data,'vL')
    data.vS = SmoothVelocity(data,'vS')
    data.vCM = SmoothVelocity(data,'vCM')
    
    return data

def SmoothVelocity(data,varName):
    #After plotting the raw data for forces, the term mdu/dt was very noisy
    #We will apply a smoothing algorithm to all velocities in the y-direction
        
    #Use Savitsky-Golay filter from scipy
    yval = data[varName].tolist()
    varValue = savgol_filter(yval, 51, 3)
    
    return varValue

def AppendData(data):
    #'amp','Re','time','tau','COM','lowY'
    #Add tau = [-0.25,0) to the database
    #Sort data in descending order of tau
    data = data.sort_values(by=['tau'],ascending=False)
    data = data.reset_index(drop=True)
    #Parameters for appending loop
    numRows = int(0.25/dt)
    minIndex = data.index.max()+1
    maxIndex = minIndex + numRows
    posIdx = 1
    for idxRow in range(minIndex,maxIndex):
        data.loc[idxRow,'amp'] = 0.18
        data.loc[idxRow,'Re'] = data.loc[0,'Re']
        data.loc[idxRow,'time'] = data.loc[idxRow - 1,'time'] - dt
        data.loc[idxRow,'tau'] = data.loc[idxRow - 1,'tau'] - dt
        data.loc[idxRow,'COM'] = (data.loc[idxRow-1,'COM'] - (data.loc[posIdx-1,'COM'] - data.loc[posIdx,'COM']))
        data.loc[idxRow,'lowY'] = (data.loc[idxRow-1,'lowY'] - (data.loc[posIdx-1,'lowY'] - data.loc[posIdx,'lowY']))
        posIdx += 1
        
    #Add tau = (1.0,1.25] to the database
    #Sort data in ascending order of tau
    data = data.sort_values(by=['tau'],ascending=True)
    data = data.reset_index(drop=True)
    #Parmaeters for appending loop
    minIndex = data.index.max()+1
    maxIndex = minIndex + numRows
    for idxRow in range(minIndex,maxIndex):
        data.loc[idxRow,'amp'] = 0.18
        data.loc[idxRow,'Re'] = data.loc[0,'Re']
        data.loc[idxRow,'time'] = data.loc[idxRow - 1,'time'] + dt
        data.loc[idxRow,'tau'] = data.loc[idxRow - 1,'tau'] + dt
        data.loc[idxRow,'COM'] = (data.loc[idxRow-1,'COM'] + (data.loc[posIdx,'COM'] - data.loc[posIdx-1,'COM']))
        data.loc[idxRow,'lowY'] = (data.loc[idxRow-1,'lowY'] + (data.loc[posIdx,'lowY'] - data.loc[posIdx-1,'lowY']))
        posIdx += 1
    
    #data = data.sort_by(['tau'],ascending=True)
    #data = data.reset_index(drop=True)
    
    print('Length of data = ',len(data.tau))
    return data.copy()

def PlotPos(data):
    #Plot the position to make sure the appending was done properly
    fig = plt.figure(num=1,figsize=(4,6),dpi=200)
    ax1 = fig.add_subplot(211)
    ax1.set_title(r'Center of Mass',fontsize=20)
    #ax1.set_xlabel(r'time $\tau$',fontsize=14)
    ax1.set_ylabel(r'$y_{CM}$ m',fontsize=14)
    ax1.plot(data.tau,data.COM)
    ax1.set_xlim(-0.25,1.25)
    ax2 = fig.add_subplot(212)
    ax2.set_title(r'Small Sphere',fontsize=20)
    ax2.set_xlabel(r'time $\tau$',fontsize=14)
    ax2.set_ylabel(r'$y_{r}$ m',fontsize=14)
    ax2.plot(data.tau,data.lowY)
    ax2.set_xlim(-0.25,1.25)
    fig.tight_layout()
    fig.clf()
    plt.close()
    
    #sys.exit(0)
    return

def PlotVel(data):
    #Plot the position to make sure the appending was done properly
    fig = plt.figure(num=2,figsize=(4,6),dpi=200)
    ax1 = fig.add_subplot(211)
    ax1.set_title(r'Center of Mass',fontsize=20)
    #ax1.set_xlabel(r'time $\tau$',fontsize=14)
    ax1.set_ylabel(r'$v_{CM}$ m/s',fontsize=14)
    ax1.plot(data.tau,data.VelSmooth)
    ax1.set_xlim(-0.25,1.25)
    ax2 = fig.add_subplot(212)
    ax2.set_title(r'Small Sphere',fontsize=20)
    ax2.set_xlabel(r'time $\tau$',fontsize=14)
    ax2.set_ylabel(r'$v_{r}$ m/s',fontsize=14)
    ax2.plot(data.tau,data.lowVySmooth)
    ax2.set_xlim(-0.25,1.25)
    fig.tight_layout()
    fig.clf()
    plt.close()
    
    #sys.exit(0)
    return

def FindRoots(data):
    #Identify 2 locations
    #1) where v goes from + to -
    #2) where v goes from - to +
    #+ to - v will always happen first
    nTime = len(data)
    idxTime = 1
    rootList = []
    
    #Find where v either
    #goes from + to -
    #goes from - to +
    #Save each root in rootList
    while(idxTime < nTime):
        if(data[idxTime]*data[idxTime-1] <= 0.0 and 
           data[idxTime] >= 0.0):
            #- to + root has been found
            print('- to +: b4Vel = %.3f\tvel = %.3f\ta4vel = %.3f'%(data[idxTime-1],
                                                            data[idxTime],
                                                            data[idxTime+1]))
            #Save + to - index value
            rootList.append(idxTime)
        elif(data[idxTime]*data[idxTime-1] <= 0.0 and 
               data[idxTime] <= 0.0):
            #+ to - root has been found
            print('+ to -: b4Vel = %.3f\tvel = %.3f\ta4vel = %.3f'%(data[idxTime-1],
                                                                data[idxTime],
                                                                data[idxTime+1]))
            #Save- to + index value
            rootList.append(idxTime)
        idxTime += 1
        
    return rootList

def PlotPowRec(data,rootList,name):
    global FIGNUM,FREQ
    #Plot the position to make sure the appending was done properly
    fig = plt.figure(num=4,figsize=(6,2),dpi=200)
    #Velocity
    ax = fig.add_subplot(111)
    csfont = {'fontname':'Times New Roman'}
    #ax.set_title(r'Re = %.2f'%ReValue,fontsize=20,**csfont)
    ax.set_xlabel(r'time $\tau$',fontsize=14,**csfont)
    #ax.set_ylabel(r'Re = %.1f'%ReValue,fontsize=14,**csfont)
    ax.set_ylabel(r'$v_{CM}$',fontsize=14,**csfont)
    
    ax2 = ax.twinx() # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel(r'$v_{r}$',fontsize=14,**csfont,color='tab:blue')
    #Important Parameters
    NetDisp = data.loc[999,'yCM'] - data.loc[0,'yCM']
    NetVel = NetDisp*FREQ
    print('NetDisp = ',NetDisp)
    print('AvgVel = ',NetDisp*10.0)
    if(NetDisp <= 0.0):
        #SSL
        #Power = Green
        #Recovery = Red
        colorList = ['green','red']
    else:
        #LSL
        #Power = Green
        #Recovery = Red
        colorList = ['red','green']
    BODYLENGTH = RSMALL+RLARGE+RSL_DEFAULT*float(rslValue)
    #Horizontal Lines
    ax.plot([-0.25,1.25],[0.0,0.0] ,color='k',lw=2,alpha=0.25)
    ax.plot([-0.25+dt*rootList[0],-0.25+dt*rootList[2]],[0.0,0.0],color='k',lw=2)
    #Data
    ax.plot(data.tau,data.vCM,color='k',lw=2,alpha=0.25)
    ax.plot(data.loc[rootList[0]:rootList[2],'tau'],data.loc[rootList[0]:rootList[2],'vCM'],color='k',lw=2)
    #Small Sphere
    #ax2.plot(data.tau,data.lowVySmooth/BODYLENGTH,lw=2,ls='-',color='tab:blue',zorder=0)
    ax2.plot(data.tau,data.vS*NetVel,color='tab:blue',lw=2,alpha=0.25)
    ax2.plot(data.loc[rootList[0]:rootList[2],'tau'],data.loc[rootList[0]:rootList[2],'vS']*NetVel,color='tab:blue',lw=2,ls='-')

    #Add Dotted Black Vertical lines where velocity switches direction
    minY, maxY = -1.4, 1.4
    ax.plot([data.loc[rootList[0],'tau'],data.loc[rootList[0],'tau']],[0.0,data.loc[rootList[0],'vCM']],color='k',lw=2,zorder=2)
    ax.plot([data.loc[rootList[1],'tau'],data.loc[rootList[1],'tau']],[0.0,data.loc[rootList[1],'vCM']],color='k',lw=2,zorder=2)
    ax.plot([data.loc[rootList[2],'tau'],data.loc[rootList[2],'tau']],[0.0,data.loc[rootList[2],'vCM']],color='k',lw=2,zorder=2)
    ax.plot([-0.25+dt*rootList[0],-0.25+dt*rootList[2]],[minY,minY] ,color='k',lw=2)
    ax.plot([-0.25+dt*rootList[0],-0.25+dt*rootList[2]],[maxY,maxY] ,color='k',lw=2)
    
    #Let's fill between the curves
    #Test on the velocity curve
    ax = FillBetween(ax,rootList[0],rootList[1],data.vCM,colorList[1])
    ax = FillBetween(ax,rootList[1],rootList[2],data.vCM,colorList[0])
  
    #Set axes limits
    ax.axis([-0.25,1.25,minY,maxY])
    ax2.set_ylim(-4.5,4.5)
    
    #Plot Vertical Lines
    PlotVertLines(ax,0.0,minY,maxY)
    
    #Change Axes Parameters
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    ax.set_axisbelow(False)
    ax.set_xticks(np.arange(-0.25,1.26,step=0.25))
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    #Ax2
    ax2.tick_params(which='major',axis='y',direction='in',length=6,width=1,labelcolor='tab:blue')
    ax2.tick_params(which='minor',axis='y',direction='in',length=4,width=0.75,labelcolor='tab:blue')
    ax2.yaxis.set_major_locator(MultipleLocator(3.0))
    ax2.yaxis.set_minor_locator(MultipleLocator(1.5))
    
    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'/../PaperFigures/Results/PR/'+str(name)+'/PRQual_'+str(name)+'_vavg.svg')
    fig.clf()
    plt.close()
    FIGNUM += 1
    
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

def FillBetween(ax,start,end,vel,c):
    #Create x coords between start and end
    xval = np.linspace(-0.25+dt*start,-0.25+dt*end,end-start+1)
    y1 = vel[start:end+1]
    y2 = np.zeros(len(xval))
    ax.fill_between(xval,y1,y2,color=c,alpha=0.5)
    
    return ax

if __name__ == '__main__':
    #This is where the main part of the code will be conducted
    #Functions will be called here
    
    '''STOKES FIGURE'''
    #Obtain Stokes Data directory
    cwd_DATA = cwd_PYTHON+'/../Stokes/Data'
    #Store Data in database
    stokesData = StoreData(cwd_DATA)
    #keep data of importance
    dataDict = {'time':stokesData.time,'Ly':stokesData.yUp,'Sy':stokesData.yLow}
    data = pd.DataFrame(data=dataDict)
    #Create new important variables
    SSTIME = 1.0
    data['tau'] = data.time/PERIOD
    data.tau = data.tau - 10.0*(SSTIME + 0.075)
    data['yCM'] = 0.8*data.Ly + 0.2*data.Sy
    #Calculate velocities of CM and Small Sphere
    data = CalcVel(data)
    #vL, vS, and vCM calculated
    #Subset data to 1.5 periods (tau = 0.0 is when time = 1.075s)
    plotData = data[10500:12000].copy()
    plotData = plotData.reset_index(drop=True)
    print('Stokes Data:')
    print(len(plotData.time))
    print('Answer should be 1500')
    print('Delta Time = ',plotData.loc[len(plotData.time)-1,'time'] - plotData.loc[0,'time'])
    #Use Rootfinder to locate where v_r switches directions
    rootList = FindRoots(plotData.vS)
    print('List of v_r roots')
    print(rootList)
    #Use the root list to color the vel CM curves
    PlotPowRec(plotData,rootList,'Stokes')
    
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
        #Use Rootfinder to locate where v_r switches directions
        rootList = FindRoots(plotData.vS)
        print('List of v_r roots')
        print(rootList)
        #Use the root list to color the vel CM curves
        PlotPowRec(plotData,rootList,name)
    