#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:34:36 2019

@author: thomas
"""

#Data is coming from posRSL#.#.csv located in ../CoMData/RSL#.#
#We will look try to plot the large sphere's displacement and velocity over time
#Different figures for RSL. Different subplots for Amplitude

#1) Obtain database from directory ../CoMData/RSL#.#/
#2) Loop over rsl directories
#3) Split database up based on the amplitude
#4) Use the plotting structure from bigplots

#MODULES
import os, sys
import numpy as np
import pandas as pd
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)

#CONSTANTS
PERIOD = 0.1
FREQ = 10.0
DENS = 2.0
MAXAMP = 0.3
RSMALL = 0.15
RLARGE = 0.3
RSL_DEFAULT = 0.75
dt = 1.0e-4
fignum=0
allPosData, allVelData = [],[]
posData, velData = [],[]
posNameList = ['Ly','Sy','COM']
velNameList = ['VelLarge','VelSmall','VelCOM']
listRSL = ['1.0','1.1','1.2','1.3','1.4','1.5','1.6','1.8','2.0']

figList, axList = [None]*6, [None]*6

data, ampData, ReData = [None], [None]*5, [None]*20
ReListAll = [None]*5

def MakeDirectory(directory,seed):
    if not os.path.exists(directory+'/'+str(seed)):
        os.makedirs(directory+'/'+str(seed))
    return

def CalcVel(data,posName,velName):
    data[velName] = 0.0
    nTime = len(data.time)
    for idxTime in range(nTime):
        if((idxTime >= 1) and (idxTime < nTime-1)):
            #Finite difference: midpoint method
            data.loc[idxTime,velName] = (data.loc[idxTime+1,posName] - 
                         data.loc[idxTime-1,posName])/(data.loc[idxTime+1,'time'] - data.loc[idxTime-1,'time'])
    #Calculate end point values
    data.loc[0,velName] = data.loc[1,velName]
    data.loc[nTime - 1,velName] = data.loc[nTime - 2,velName]
    return data.copy()

def PlotDataAll(ax, data, plotParams, plotType, idxRSL,posName):
    global fignum,ampData,posNameList, listRSL
    #Change directory to BigPlots directory
    cwd_CURRENT = os.getcwd()
    os.chdir(cwd_CURRENT+'/../../Figures/Movement/BigPlots/')
    cwd_BP = os.getcwd()
    varName = plotParams.varName
    print('varName = ',varName)
    MakeDirectory(cwd_BP,varName)
    cwd_DATA = cwd_BP+'/'+varName+'/'
    os.chdir(cwd_DATA)
    print('cwd_DATA = ',cwd_DATA)
    #Color Palette for plotting
    R1=[255/255,255/255,153/255,153/255,204/255]
    G1=[153/255,204/255,255/255,204/255,153/255]
    B1=[204/255,153/255,153/255,255/255,255/255]
    idxLC = 0
    #3)Split database up based on amplitude
    #ampData[0] => A0.12 data
    #ampData[1] => A0.15 data
    #ampData[2] => A0.18 data
    #ampData[3] => A0.21 data
    #ampData[4] => A0.24 data
    for idxAmp in range(len(ampData)):
        #Set up subplot
        idxRow = idxRSL
        #print('idxRow = ',idxRow)
        idxCol = idxAmp
        if(varName in posNameList):
            if(idxCol == 0):
                ax[idxRow,idxCol].set_ylabel(r'$\hat{d}_0$ = %.1f: $\Delta y$ (m)'%(float(listRSL[idxRSL])*RSL_DEFAULT/RSMALL),fontsize=20)
            if(idxRow == 8):
                ax[idxRow,idxCol].set_xlabel(r'time $\tau$',fontsize=20)
        else:
            if(idxCol == 0):
                ax[idxRow,idxCol].set_ylabel(r'$\hat{d}_0$ = %.1f: v (m/s)'%(float(listRSL[idxRSL])*RSL_DEFAULT/RSMALL),fontsize=20)
            if(idxRow == 8):
                ax[idxRow,idxCol].set_xlabel(r'time $\tau$',fontsize=20)
        ShadeValue = 1 #used to scale line color so there is a gradient as rsl changes
        #Calculate the amplitude used in the subplot
        ampValue = 0.12 + 0.03*idxAmp
        ampData[idxAmp] = data[data.amp == ampValue]
        #Store in a temporary database
        tempAmpData = ampData[idxAmp]
        tempAmpData = tempAmpData.sort_values(by=['Re','tau'])
        tempAmpData = tempAmpData.reset_index(drop=True)
        #Get List of Re Values in tempData
        ReList = tempAmpData.Re.values.tolist()
        ReList = sorted(list(set(ReList)), reverse = True) #ascending or descending order???
        #Plot Horizontal axis
        ax[idxRow,idxCol].plot([0.0,1.0],[0.0,0.0],'k',linewidth=2)
        #Loop over Re values
        for idxRe in range(1,len(ReList)): #5 distinct Re values
            #Select RGB Color
            R=R1[idxLC]*ShadeValue
            G=G1[idxLC]*ShadeValue
            B=B1[idxLC]*ShadeValue
            #Save Re filtered data
            ReData[idxRe] = tempAmpData[tempAmpData.Re == ReList[idxRe]]
            #Store in a temporary database
            tempReData = ReData[idxRe]
            tempReData = tempReData.reset_index(drop=True)
            #Plot disp (y) or vel (v) vs time
            if(len(tempReData.tau) < 1010):
                if(varName in posNameList):
                    ax[idxRow,idxCol].plot(tempReData.tau,tempReData[varName] - tempReData.loc[0,varName], 
                            color = (R,G,B), linewidth = 2,label = round(ReList[idxRe],2))
                else:
                    tempReData = CalcVel(tempReData,posName,varName)
                    tempReData[varName] = SmoothCurve(tempReData,varName)
                    ax[idxRow,idxCol].plot(tempReData.tau,tempReData[varName], 
                            color = (R,G,B), linewidth = 2, label = round(ReList[idxRe],2))
            #Decrease ShadeValue
            ShadeValue -= 0.025
        #Add Vertical Dashed Lines every quarter period
        PlotVertLines(ax[idxRow,idxCol],0.0,plotParams['min'],plotParams['max'])
        #Increment line color to match amplitude
        idxLC += 1
    
        #Axes Parameters
        ax[idxRow,idxCol].tick_params(which='major',axis='both',direction='in',length=8,width=2,labelsize=14)
        ax[idxRow,idxCol].tick_params(which='minor',axis='both',direction='in',length=6,width=1.5)
        ax[idxRow,idxCol].set_xticks(np.arange(0.0,1.1,step=0.2))
        ax[idxRow,idxCol].xaxis.set_minor_locator(MultipleLocator(0.1))
        #Change Y-axis numbers major and minor
        #Ly
        if(varName == 'Ly'):
            ax[idxRow,idxCol].set_yticks(np.arange(0.0,0.25,step=0.1))
            ax[idxRow,idxCol].yaxis.set_minor_locator(MultipleLocator(0.05))
        #Sy
        elif(varName == 'Sy'):
            ax[idxRow,idxCol].set_yticks(np.arange(-0.6,0.25,step=0.2))
            ax[idxRow,idxCol].yaxis.set_minor_locator(MultipleLocator(0.1))
        #COM
        elif(varName == 'COM'):
            ax[idxRow,idxCol].set_yticks(np.arange(-0.08,0.13,step=0.04))
            ax[idxRow,idxCol].yaxis.set_minor_locator(MultipleLocator(0.02))
        #vLy
        elif(varName == 'VelLarge'):
            ax[idxRow,idxCol].set_yticks(np.arange(-6.0,6.0,step=2.0))
            ax[idxRow,idxCol].yaxis.set_minor_locator(MultipleLocator(1.0))
        #vSy
        elif(varName == 'VelSmall'):
            ax[idxRow,idxCol].set_yticks(np.arange(-20.0,20.0,step=5.0))
            ax[idxRow,idxCol].yaxis.set_minor_locator(MultipleLocator(2.5))
        #vCOM
        else:
            ax[idxRow,idxCol].set_yticks(np.arange(-2.0,2.0,step=1.0))
            ax[idxRow,idxCol].yaxis.set_minor_locator(MultipleLocator(0.5))
        #ax[idxRow,idxCol].yaxis.set_minor_locator(MultipleLocator(0.01))
        #ax[idxRow,idxCol].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax[idxRow,idxCol].set_ylim(plotParams['min'],plotParams['max'])
        ax[idxRow,idxCol].set_xlim(0.0,1.0)
        
    os.chdir(cwd_CURRENT)
    print('cwd = ',os.getcwd())
    return cwd_BP

def PlotVertLines(ax,xMin,yMin,yMax):
    #Vertical Lines Indicating Expansion/Compression
    ax.plot([xMin+0.25,xMin+0.25],[yMin,yMax],color='gray',ls='--',zorder=0)
    ax.plot([xMin+0.5,xMin+0.5],[yMin,yMax],color='k',ls='--',zorder=0)
    ax.plot([xMin+0.75,xMin+0.75],[yMin,yMax],color='gray',ls='--',zorder=0)
    return

def SmoothCurve(data,varName):
    #After plotting the raw data for forces, the term mdu/dt was very noisy
    #We will apply a smoothing algorithm to all velocities in the y-direction
        
    #Use Savitsky-Golay filter from scipy
    yval = data[varName].tolist()
    varValue = savgol_filter(yval, 51, 3)
    
    return varValue

if __name__ == '__main__':
    #Obtain current directory
    cwd_PYTHON = os.getcwd()
    #Change to data directories
    os.chdir(cwd_PYTHON+'/../CoMData')
    #1) Save data directory
    cwd_DATA = os.getcwd()
    #Specify RSL Directories and values to be analyzed
    #All directories
    dirsRSL = ['RSL1.0','RSL1.1','RSL1.2','RSL1.3','RSL1.4','RSL1.5','RSL1.6','RSL1.8','RSL2.0']
    #1 directory
    #dirsRSL = ['RSL1.0','RSL1.1']
    
    #Create Figures for Plotting Position and Velocity Data
    #Large Sphere, Small Sphere, COM
    #In this order: Ld,Sd,COMd,Lv,Sv,COMv
    nRows = 9
    nCols = 5
    #Ld
    figLd, axLd = plt.subplots(nrows=nRows, ncols=nCols, num=0,figsize=(20,27),dpi=200)
    #figLd.suptitle('Large Sphere Displacement', fontsize=20,x=0.5,y=1.1)
    #bigAxLd = figLd.add_subplot(111)
    #bigAxLd.set_title('Large Sphere Displacement')
    #bigAxLd.tick_params(labelcolor=(1.,1.,1.,0.0),top='off',bottom='off',left='off',right='off')
    #bigAxLd._frameon = False
    #Sd
    figSd, axSd = plt.subplots(nrows=nRows, ncols=nCols, num=1,figsize=(20,27),dpi=200)
    #figSd.suptitle('Small Sphere Displacement', fontsize=20)
    #COMd
    figCd, axCd = plt.subplots(nrows=nRows, ncols=nCols, num=2,figsize=(20,27),dpi=200)
    #figCd.suptitle('COM Displacement', fontsize=20)
    #Lv
    figLv, axLv = plt.subplots(nrows=nRows, ncols=nCols, num=3,figsize=(20,27),dpi=200)
    #figLv.suptitle('Large Sphere Velocity', fontsize=20)
    #Sv
    figSv, axSv = plt.subplots(nrows=nRows, ncols=nCols, num=4,figsize=(20,27),dpi=200)
    #figSv.suptitle('Small Sphere Velocity', fontsize=20)
    #COMv
    figCv, axCv = plt.subplots(nrows=nRows, ncols=nCols, num=5,figsize=(20,27),dpi=200)
    #figCv.suptitle('COM Velocity', fontsize=20)
    
    #Save figures in a list for future use
    figList[0] = figLd
    figList[1] = figSd
    figList[2] = figCd
    figList[3] = figLv
    figList[4] = figSv
    figList[5] = figCv
    
    #Create Subplot titles for all plots
    for idx in range(len(ampData)):
        ampValue = 0.12 + 0.03*idx
        axLd[0,idx].set_title(r'$\Delta y_R$: $\epsilon$ = %.1f'%(ampValue/RSMALL),fontsize=20)
        axSd[0,idx].set_title(r'$\Delta y_r$: $\epsilon$ = %.1f'%(ampValue/RSMALL),fontsize=20)
        axCd[0,idx].set_title(r'$\Delta y_{CM}$: $\epsilon$ = %.1f'%(ampValue/RSMALL),fontsize=20)
        axLv[0,idx].set_title(r'$v_R$: $\epsilon$ = %.1f'%(ampValue/RSMALL),fontsize=20)
        axSv[0,idx].set_title(r'$v_r$: $\epsilon$ = %.1f'%(ampValue/RSMALL),fontsize=20)
        axCv[0,idx].set_title(r'$v_{CM}$: $\epsilon$ = %.1f'%(ampValue/RSMALL),fontsize=20)
    
    '''#Add Row Annotation
    pad = 5 # in points
    rows = ['RSL={}'.format(row) for row in ['1.0','1.1','1.2','1.3','1.4',
            '1.5','1.6','1.8','2.0']]
    for idx in range(9):
        for row in rows:
            axLd[idx,0].annotate(row, xy=(0,0.5), xytext=(0,pad),
                                 textcoords='offset points', size='large', ha='right',va='center')'''
    
    #Save axes in a list for future use
    axList[0] = axLd
    axList[1] = axSd
    axList[2] = axCd
    axList[3] = axLv
    axList[4] = axSv
    axList[5] = axCv        
    
    #2) Loop over RSL directories
    for idxRSL in range(len(dirsRSL)):        
        #Change to RSL#.# directory
        os.chdir(cwd_DATA+'/'+dirsRSL[idxRSL])
        #Save directory of velocity data
        cwd_VEL = os.getcwd()
        #Store data from velRSL#.#.csv file
        allPosData = pd.read_csv('pos'+dirsRSL[idxRSL]+'Shifted.csv',delimiter = ',')
        allVelData = pd.read_csv('vel'+dirsRSL[idxRSL]+'Shifted.csv',delimiter = ',')
        #keep data of importance
        posDataDict = {'amp':allPosData.amp,'Re':allPosData.Re,'time':allPosData.time,
                       'tau':allPosData.tau,'Ly':allPosData.upCOMy,'Sy':allPosData.lowCOMy,
                       'COM':allPosData.COM}
        posData = pd.DataFrame(data=posDataDict)
        velDataDict = {'amp':allVelData.amp,'Re':allVelData.Re,'tau':allVelData.tau,
                       'VelLarge':allVelData.VelLarge,'VelSmall':allVelData.VelSmall,
                       'VelCOM':allVelData.VelCOM}
        velData = pd.DataFrame(data=velDataDict)
        '''#Check End points of velocity
        print(velData.head())
        print(velData.tail())
        #Smooth the Velocity Data!
        velData['VelLarge'] = SmoothCurve(velData,'VelLarge')
        velData['VelSmall'] = SmoothCurve(velData,'VelSmall')
        velData['VelCOM'] = SmoothCurve(velData,'VelCOM')
        print(velData.head())
        print(velData.tail())
        plt.close()
        sys.exit(0)'''
        #Plot the displacement over time for each part of the spherobot
        #Make a dataframe of plot parameters
        plotName = ['Ly','Sy','COM','VelLarge','VelSmall','VelCOM']
        #plotMin = [-0.1,-0.25,-0.05,-5.0,-20.0,-1.5]
        #plotMax = [0.125,0.35,0.12,5.0,20.0,1.5]
        plotMin = [-0.045,-0.45,-0.05,-5.0,-17.5,-1.75]
        plotMax = [0.175,0.15,0.1,5.0,17.5,1.75]
        plotDict = {'varName':plotName,'min':plotMin,'max':plotMax}
        plotParams = pd.DataFrame(data=plotDict)
        
        #Plot all RSL on same plot
        plotType = 'same_all'
        for idx in range(len(plotParams.varName)):
            print('varName index = ',idx)
            if(idx < 3):
                cwd_BP = PlotDataAll(axList[idx],posData,plotParams.loc[idx,:],plotType,idxRSL,None)
            else:
                cwd_BP = PlotDataAll(axList[idx],posData,plotParams.loc[idx,:],plotType,idxRSL,plotName[idx-3])
    #Loop over variable names
    #Change to appropriate directory
    #Tight Layout figure
    #Save figure in that directory
    for idx in range(len(plotParams.varName)):
        varName = plotParams.loc[idx,'varName']
        os.chdir(cwd_BP+'/'+varName+'/')
        figList[idx].tight_layout()
        figList[idx].savefig(varName+'_'+plotType+'_Final.png')
        figList[idx].clf()
