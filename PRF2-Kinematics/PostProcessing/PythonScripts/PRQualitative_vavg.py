#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 11:31:12 2020

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

def StoreData(cwd_RSL,rsl):
    #Reset position data every Re
    allData = []
    #Load position data
    allData = pd.read_csv(cwd_RSL+'/pos'+rsl+'Shifted.csv',delimiter = ',')
    return allData

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

def CalcVel(data):
    DT = 1.0e-4
    data['VelCOM'], data['lowVy'] = 0.0, 0.0
    nTime = len(data.time)
    for idxTime in range(nTime):
        if((idxTime >= 1) and (idxTime < nTime-1)):
            #Finite difference: midpoint method
            data.loc[idxTime,'VelCOM'] = (data.loc[idxTime+1,'COM'] - 
                         data.loc[idxTime-1,'COM'])/(2.0*DT)
            data.loc[idxTime,'lowVy'] = (data.loc[idxTime+1,'lowY'] - 
                         data.loc[idxTime - 1,'lowY'])/(2.0*DT)
    #Calculate end point values
    data.loc[0,'VelCOM'] = data.loc[1,'VelCOM']
    data.loc[0,'lowVy'] = data.loc[1,'lowVy']
    data.loc[nTime - 1, 'VelCOM'] = data.loc[nTime - 2,'VelCOM']
    data.loc[nTime - 1, 'lowVy'] = data.loc[nTime - 2,'lowVy']
    return data.copy()

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

def PlotPowRec(data,rootList,ReValue,rsl,rslValue):
    global FIGNUM,FREQ
    #Plot the position to make sure the appending was done properly
    fig = plt.figure(num=4,figsize=(6,2),dpi=200)
    #Velocity
    ax = fig.add_subplot(111)
    csfont = {'fontname':'Times New Roman'}
    #ax.set_title(r'Re = %.1f'%ReValue,fontsize=12,**csfont)
    ax.set_xlabel(r'time $\tau$',fontsize=14,**csfont)
    #ax.set_ylabel(r'Re = %.1f'%ReValue,fontsize=14,**csfont)
    ax.set_ylabel(r'$v_{CM}$',fontsize=14,**csfont)
    
    ax2 = ax.twinx() # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel(r'$v_{r}$',fontsize=14,**csfont,color='tab:blue')
    #Important Parameters
    NetDisp = data.loc[999,'COM'] - data.loc[0,'COM']
    NetVel = NetDisp*FREQ
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
    #ax.plot([-0.25,1.25],[0.0,0.0] ,color='k',lw=2,zorder=2)
    ax.plot([-0.25,1.25],[0.0,0.0] ,color='k',lw=2,alpha=0.25)
    ax.plot([-0.25+dt*rootList[0],-0.25+dt*rootList[2]],[0.0,0.0],color='k',lw=2)
    #Data
    #ax.plot(data.tau,data.VelSmooth/BODYLENGTH,color='k',lw=2,zorder=2)
    ax.plot(data.tau,data.VelSmooth,color='k',lw=2,alpha=0.25)
    ax.plot(data.loc[rootList[0]:rootList[2],'tau'],data.loc[rootList[0]:rootList[2],'VelSmooth'],color='k',lw=2)
    #Small Sphere
    #ax2.plot(data.tau,data.lowVySmooth/BODYLENGTH,lw=2,ls='-',color='tab:blue',zorder=0)
    ax2.plot(data.tau,data.lowVySmooth*NetVel,color='tab:blue',lw=2,alpha=0.25)
    ax2.plot(data.loc[rootList[0]:rootList[2],'tau'],data.loc[rootList[0]:rootList[2],'lowVySmooth']*NetVel,color='tab:blue',lw=2,ls='-')

    #ax2.plot(data.tau,data.lowVySmooth,color='b',lw=2)
    #ax2.plot(data.loc[:rootList[0]+1,'tau']             ,data.loc[:rootList[0]+1,'VelSmooth']             ,color='k',lw=2)
    #ax2.plot(data.loc[rootList[0]+1:rootList[1]+1,'tau'],data.loc[rootList[0]+1:rootList[1]+1,'VelSmooth'],color='k',lw=2)
    #ax2.plot(data.loc[rootList[1]+1:rootList[2]+1,'tau'],data.loc[rootList[1]+1:rootList[2]+1,'VelSmooth'],color='k',lw=2)
    #ax2.plot(data.loc[rootList[2]:,'tau']               ,data.loc[rootList[2]:,'VelSmooth']               ,color='k',lw=2)
    #Add Dotted Black Vertical lines where velocity switches direction
    minY, maxY = -1.4, 1.4
    ax.plot([data.loc[rootList[0],'tau'],data.loc[rootList[0],'tau']],[0.0,data.loc[rootList[0],'VelSmooth']],color='k',lw=2,zorder=2)
    ax.plot([data.loc[rootList[1],'tau'],data.loc[rootList[1],'tau']],[0.0,data.loc[rootList[1],'VelSmooth']],color='k',lw=2,zorder=2)
    ax.plot([data.loc[rootList[2],'tau'],data.loc[rootList[2],'tau']],[0.0,data.loc[rootList[2],'VelSmooth']],color='k',lw=2,zorder=2)
    #ax.plot([data.loc[rootList[0],'tau'],data.loc[rootList[0],'tau']],[minY,maxY],color='k',lw=2)
    #ax.plot([data.loc[rootList[1],'tau'],data.loc[rootList[1],'tau']],[minY,maxY],color='k',lw=2)
    #ax.plot([data.loc[rootList[2],'tau'],data.loc[rootList[2],'tau']],[minY,maxY],color='k',lw=2)
    ax.plot([-0.25+dt*rootList[0],-0.25+dt*rootList[2]],[minY,minY] ,color='k',lw=2)
    ax.plot([-0.25+dt*rootList[0],-0.25+dt*rootList[2]],[maxY,maxY] ,color='k',lw=2)
    
    #Let's fill between the curves
    #Test on the velocity curve
    #ax2 = FillBetween(ax2,0          ,rootList[0],data.VelSmooth,colorList[0])
    ax = FillBetween(ax,rootList[0],rootList[1],data.VelSmooth,colorList[1])
    ax = FillBetween(ax,rootList[1],rootList[2],data.VelSmooth,colorList[0])
    #ax2 = FillBetween(ax2,rootList[2],1500       ,data.VelSmooth,colorList[1])
    
    ax.axis([-0.25,1.25,minY,maxY])
    ax2.set_ylim(-4.5,4.5)
    #ax2.set_ylim(minY,maxY)
    
    #Plot Vertical Lines
    PlotVertLines(ax,0.0,minY,maxY)
    
    #Change Axes Parameters
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    ax.set_axisbelow(False)
    ax.set_xticks(np.arange(-0.25,1.26,step=0.25))
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    ax2.tick_params(which='major',axis='y',direction='in',length=6,width=1,labelcolor='tab:blue')
    ax2.tick_params(which='minor',axis='y',direction='in',length=4,width=0.75,labelcolor='tab:blue')
    ax2.yaxis.set_major_locator(MultipleLocator(3.0))
    ax2.yaxis.set_minor_locator(MultipleLocator(1.5))
    
    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'/../PaperFigures/Results/PR/'+rsl+'/V_vavg/PRQual_V_vavg'+str(FIGNUM)+'.svg')
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
    #y3 = np.linspace(-1.25,1.25,len(xval))
    #print('y1.shape = ',y1.shape)
    #print('y2.shape = ',y2.shape)
    ax.fill_between(xval,y1,y2,color=c,alpha=0.5)
    
    return ax

if __name__ == '__main__':
    #This is where the main part of the code will be conducted
    #Functions will be called here
    cwd_DATA = cwd_PYTHON+'/../CoMData'
    #dirsRSL = [d for d in os.listdir(cwd_DATA) if not d.startswith('.')]
    dirsRSL = ['RSL1.3']
    for rsl in dirsRSL:
        FIGNUM = 0
        rslString, rslValue = tuple(re.split('(\d.*)',rsl)[:2])
        #if(rsl == 'RSL1.3'):
        cwd_RSL = cwd_DATA+'/'+rsl
        allData = StoreData(cwd_RSL,rsl)
        #keep data of importance
        dataDict = {'amp':allData.amp,'Re':allData.Re,'time':allData.time,
                    'tau':allData.tau,'COM':allData.COM,'lowY':allData.lowCOMy}
        data = pd.DataFrame(data=dataDict)
        #Split database based upon amplitude
        #Only save A_f = 0.6 => A = 0.18m
        ampValue = 0.18
        ampData = data[data.amp == ampValue].copy()
        #Store in a temporary database
        tempAmpData = ampData.copy()
        tempAmpData = tempAmpData.sort_values(by=['Re','tau'])
        tempAmpData = tempAmpData.reset_index(drop=True)
        #Get List of Re Values in tempData
        ReList = tempAmpData.Re.values.tolist()
        ReList = sorted(list(set(ReList)))
        print(ReList)
        print(len(ReList))
        #Loop over Re values
        for idxRe in range(len(ReList)):
            #Save Re filtered data
            ReData[idxRe] = tempAmpData[tempAmpData.Re == ReList[idxRe]].copy()
            #Store in a temporary database
            tempReData = ReData[idxRe].copy()
            tempReData = tempReData.reset_index(drop=True)
            #Now that we have a single simulation's position data, we can begin
            #Append 1/4 oscillation before and after
            print('b4 append: Length = ',len(tempReData.tau))
            posReData = AppendData(tempReData)
            #DATA APPENDED CORRECTLY! WOOH
            #PlotPos(posReData)
            #Add Velocities to the database
            allReData = CalcVel(posReData)
            #Smooth Velocity Data
            #COM
            #Use Savitsky-Golay filter from scipy
            yval = allReData.VelCOM.tolist()
            allReData['VelSmooth'] = savgol_filter(yval, 51, 3)
            #lowY
            #Use Savitsky-Golay filter from scipy
            yval = allReData.lowVy.tolist()
            allReData['lowVySmooth'] = savgol_filter(yval, 51, 3)
            #See if the data looks like it should
            #PlotVel(allReData)
            #It does! Woohoo
            #Use Rootfinder to locate where v_r switches directions
            rootList = FindRoots(allReData.lowVySmooth)
            print('List of v_r roots')
            print(rootList)
            #Use the root list to color the pos/vel CM curves
            PlotPowRec(allReData,rootList,ReList[idxRe],rsl,rslValue)
            
            #sys.exit(0)
    