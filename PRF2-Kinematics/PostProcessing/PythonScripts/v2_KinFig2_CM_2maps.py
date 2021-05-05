#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 12:06:40 2020

@author: thomas
"""

#Data is coming from posRSL#.#.csv located in ../CoMData/RSL#.#
#We will look try to plot the large sphere's displacement and velocity over time
#Different figures for RSL. Different subplots for Amplitude

#1) Obtain database from directory ../CoMData/RSL1.3/
#3) Split database up based on the amplitude
#4) Plot 2 rows (top = LS, bot  = SS). Subplot size based off of maximum/minimum
#   displacement. 

#MODULES
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)
from scipy.signal import savgol_filter

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
listRSL = ['1.0','1.1','1.2','1.3','1.4','1.5','1.6','1.8','2.0']

figList, axList = [None]*2, [None]*2

data, posReData, velReData = [None], [None]*20, [None]*20
ReListAll = [None]*5

mpl.rcParams['axes.linewidth'] = 1.5 #set the value globally

def MakeDirectory(directory,seed):
    if not os.path.exists(directory+'/'+str(seed)):
        os.makedirs(directory+'/'+str(seed))
    return

def PlotDataAll(ax, posData, velData, plotParams):
    global fignum, posReData, velReData, listRSL
    #Change directory to BigPlots directory
    cwd_CURRENT = os.getcwd()
    os.chdir(cwd_CURRENT+'/../../PaperFigures/Results/')
    cwd_FIG = os.getcwd()
    #Color Palette for plotting
    R1=[255/255,255/255,153/255,153/255,204/255]
    G1=[153/255,204/255,255/255,204/255,153/255]
    B1=[204/255,153/255,153/255,255/255,255/255]
    idxLC = 2
    #3)Split database up based on amplitude
    #Axes Labels
    csfont = {'fontname':'Times New Roman'}
    ax[0].set_ylabel(r'$\Delta y_{CM}$ (m)',fontsize=16,**csfont)
    ax[1].set_ylabel(r'$v_{CM}$ (m/s)',fontsize=16,**csfont)
    ax[0].set_xlabel(r'time',fontsize=16,**csfont)
    ax[1].set_xlabel(r'time',fontsize=16,**csfont)
    ShadeValue = 1.0 #used to scale line color so there is a gradient as Re changes
    ampValue = 0.18
    #Obtain temporary position data
    tempPosData = posData[posData.amp == ampValue]
    tempPosData = tempPosData.sort_values(by=['Re','tau'])
    tempPosData = tempPosData.reset_index(drop=True)
    #Obtain temporary velocity data
    tempVelData = velData[velData.amp == ampValue]
    tempVelData = tempVelData.sort_values(by=['Re','tau'])
    tempVelData = tempVelData.reset_index(drop=True)
    #Get List of Re Values in tempData
    ReList = tempPosData.Re.values.tolist()
    ReList = sorted(list(set(ReList)), reverse = True) #descending order
    print(ReList)
    #Loop over Re Values
    RedShadeValue = 1.0
    GreenShadeValue = 1.0
    #Colormap colors
    fracRe = [a/150.0 for a in ReList]
    #ax[0].set_prop_cycle('color',[plt.cm.plasma(i) for i in fracRe])
    #ax[1].set_prop_cycle('color',[plt.cm.plasma(i) for i in fracRe])
    #ax[0].set_prop_cycle('color',[plt.cm.spring(i) for i in np.linspace(0.15,1.0,12)])
    #ax[0].set_prop_cycle('color',[plt.cm.winter(i) for i in np.linspace(0.0,1.0,7)])
    #ax[1].set_prop_cycle('color',[plt.cm.plasma(i) for i in np.linspace(1,0.15,19)])
    idxSSL, idxLSL = 12, 0
    for idxRe in range(1,len(ReList)): #5 distinct Re values
        #Save Re filtered data for pos
        posReData[idxRe] = tempPosData[tempPosData.Re == ReList[idxRe]]
        #Store in a temporary database
        tempPosReData = posReData[idxRe]
        tempPosReData = tempPosReData.reset_index(drop=True)
        #Save Re filtered data for vel
        velReData[idxRe] = tempVelData[tempVelData.Re == ReList[idxRe]]
        #Store in a temporary database
        tempVelReData = velReData[idxRe]
        tempVelReData = tempVelReData.reset_index(drop=True)
        #Smooth the Velocity Data!
        tempVelReData.VelCOM = SmoothVelocity(tempVelReData,'VelCOM')
        #Select LineStyle
        lineStyle = '-'
        #Select RGB Color
        if(tempPosReData.loc[1000,'COM'] <= tempPosReData.loc[0,'COM']):
            plotColor = plt.cm.spring(idxSSL/14.0)
            idxSSL -= 1
        else:
            plotColor = plt.cm.winter(idxLSL/7.0)
            idxLSL += 1
        #Plot 1) horizontal axis, 2) disp or vel v time
        #print('R = %.2f\tG = %.2f\tB = %.2f'%(R,G,B))
        print('plotColor = ',plotColor)
        #Displacement
        ax[0].plot(tempPosReData.tau,tempPosReData.COM - tempPosReData.loc[0,'COM'],
          linewidth=2,linestyle=lineStyle,color=plotColor)
        ax[1].plot(tempVelReData.tau,tempVelReData.VelCOM,
          linewidth=2,linestyle=lineStyle,color=plotColor)
        #Decrease ShadeValue
        #ShadeValue -= 0.035
     
    PlotLowReData(ax)
    print('were good')
        
    #Loop over Re Values
    #Dashed Overplot
    ShadeValue = 1.0
    #ax[0].set_prop_cycle('color',[plt.cm.plasma(i) for i in fracRe])
    #ax[1].set_prop_cycle('color',[plt.cm.plasma(i) for i in fracRe])
    #ax[0].set_prop_cycle('color',[plt.cm.plasma(i) for i in np.linspace(1,0.15,19)])
    #ax[1].set_prop_cycle('color',[plt.cm.plasma(i) for i in np.linspace(1,0.15,19)])
    idxLSL = 0
    for idxRe in range(1,len(ReList)): 
        #Save Re filtered data for pos
        posReData[idxRe] = tempPosData[tempPosData.Re == ReList[idxRe]]
        #Store in a temporary database
        tempPosReData = posReData[idxRe]
        tempPosReData = tempPosReData.reset_index(drop=True)
        #Save Re filtered data for vel
        velReData[idxRe] = tempVelData[tempVelData.Re == ReList[idxRe]]
        #Store in a temporary database
        tempVelReData = velReData[idxRe]
        tempVelReData = tempVelReData.reset_index(drop=True)
        #Smooth the Velocity Data!
        tempVelReData.VelCOM = SmoothVelocity(tempVelReData,'VelCOM')
        #Select LineStyle
        lineStyle = ':'
        #Select RGB Color
        idxLC = 2
        R=R1[idxLC]*ShadeValue
        G=G1[idxLC]*ShadeValue
        B=B1[idxLC]*ShadeValue
        print('R = %.2f\tG = %.2f\tB = %.2f'%(R,G,B))
        #Plot 1) horizontal axis, 2) disp or vel v time
        if(tempPosReData.loc[1000,'COM'] > tempPosReData.loc[0,'COM']):
            #Displacement
            plotColor = plt.cm.winter(idxLSL/7.0)
            ax[0].plot(tempPosReData.tau,tempPosReData.COM - tempPosReData.loc[0,'COM'],
              linewidth=2,linestyle=lineStyle,color=plotColor)
            ax[1].plot(tempVelReData.tau,tempVelReData.VelCOM,
              linewidth=2,linestyle=lineStyle,color=plotColor)
            idxLSL += 1
            #Decrease ShadeValue
            ShadeValue -= 0.125
    
    PlotStokesData(ax)    
    #Plot Horizontal Lines
    ax[0].plot([-100.0,100.0],[0.0,0.0],'k',linewidth=2)
    ax[1].plot([-100.0,100.0],[0.0,0.0],'k',linewidth=2)
    #Add Vertical Dashed Lines every quarter period
    PlotVertLines(ax[0],0.0)
    PlotVertLines(ax[1],0.0)
    
    ax[0].set_xlim(0.0,1.0)
    ax[1].set_xlim(0.0,1.0)
    ax[0].set_ylim(plotParams.loc[0,'Min'],plotParams.loc[0,'Max'])
    ax[1].set_ylim(plotParams.loc[1,'Min'],plotParams.loc[1,'Max'])
    #Increment line color to match amplitude
    #idxLC += 1
    
    #Axes Parameters
    ax[0].tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax[0].tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    ax[0].set_axisbelow(False)
    ax[0].set_xticks(np.arange(0.0,1.1,step=0.2))
    ax[0].xaxis.set_minor_locator(MultipleLocator(0.1))
    ax[0].yaxis.set_minor_locator(MultipleLocator(0.01))
    ax[1].tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax[1].tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    ax[1].set_axisbelow(False)
    ax[1].set_xticks(np.arange(0.0,1.1,step=0.2))
    ax[1].xaxis.set_minor_locator(MultipleLocator(0.1))
    ax[1].yaxis.set_minor_locator(MultipleLocator(0.25))
    
    os.chdir(cwd_CURRENT)
    print('cwd = ',os.getcwd())
    return cwd_FIG

def PlotVertLines(ax,xMin):
    #Vertical Lines Indicating Expansion/Compression
    ax.plot([xMin+0.25,xMin+0.25],[-1000,1000],color='gray',ls='--')
    ax.plot([xMin+0.5,xMin+0.5],[-1000,1000],color='k',ls='--')
    ax.plot([xMin+0.75,xMin+0.75],[-1000,1000],color='gray',ls='--')
    return

def PlotStokesData(ax):
    #Obtain current directory
    cwd_CURRENT = os.getcwd()
    os.chdir(cwd_CURRENT+'/../../Stokes/Data')
    cwd_STOKES = os.getcwd()+'/'
    #Load Data from Re directory
    posData,velData = StoreStokesData(cwd_STOKES)
    
    #Stokes Disp
    posData['COM'] = 0.2*posData.yLow + 0.8*posData.yUp
    ax[0].plot(posData.tau,posData.COM - posData.loc[750,'COM'],color='k',linewidth=2)
    
    #Stokes Vel
    velData['velCOM'] = 0.2*velData.vyLow + 0.8*velData.vyUp
    ax[1].plot(velData.tau,velData.velCOM,color='k',linewidth=2)
    os.chdir(cwd_CURRENT)
    
    return

def PlotLowReData(ax):
    global cwd_PYTHON
    #Obtain current directory
    cwd_CURRENT = os.getcwd()
    os.chdir(cwd_PYTHON)
    cwd_LOW = os.getcwd()+'/'
    '''#Load Data from Re directory
    posData,velData = StoreLowReData(cwd_LOW,'pd_Re0_01.txt')
    #LowRe: Re = 0.01
    posData['COM'] = 0.2*posData.yLow + 0.8*posData.yUp
    ax[0].plot(posData.tau,posData.COM - posData.loc[750,'COM'],color=(0.30,0.18,0.24),linewidth=2)
    #Vel
    velData['velCOM'] = 0.2*velData.vyLow + 0.8*velData.vyUp
    ax[1].plot(velData.tau,velData.velCOM,color=(0.30,0.18,0.24),linewidth=2)'''
    
    #Load Data from Re directory: Re = 0.1
    posData,velData = StoreLowReData(cwd_LOW,'pd_Re0_1.txt')
    #LowRe: Re = 0.1
    posData['COM'] = 0.2*posData.yLow + 0.8*posData.yUp
    ax[0].plot(posData.tau,posData.COM - posData.loc[750,'COM'],color=plt.cm.spring(0.05),linewidth=2)
    #Vel
    velData['velCOM'] = 0.2*velData.vyLow + 0.8*velData.vyUp
    ax[1].plot(velData.tau,velData.velCOM,color=plt.cm.spring(0.05),linewidth=2)
    print('plotColor = ',plt.cm.spring(0.05))

    #Load Data from Re directory: Re = 0.3
    posData,velData = StoreLowReData(cwd_LOW,'pd_Re0_3.txt')
    #LowRe: Re = 0.3
    posData['COM'] = 0.2*posData.yLow + 0.8*posData.yUp
    ax[0].plot(posData.tau,posData.COM - posData.loc[750,'COM'],color=plt.cm.spring(0.10),linewidth=2)
    #Vel
    velData['velCOM'] = 0.2*velData.vyLow + 0.8*velData.vyUp
    ax[1].plot(velData.tau,velData.velCOM,color=plt.cm.spring(0.10),linewidth=2)
    print('plotColor = ',plt.cm.spring(0.10))
    os.chdir(cwd_CURRENT)
    
    return
    
def StoreStokesData(dirName):
    global ssTime, DT
    DT = 1.0e-4
    ssTime = 1.0
    #Load position and force data
    pdData = []
    print('dirName = ',dirName)
    pdData = pd.read_csv(dirName+'pd.txt',delimiter = ' ')
    #Limit the amount of data to 2 periods worth
    #Start and Stop Indices
    idxStart = int(ssTime/DT)
    idxStop = int((ssTime+0.4)/DT)
    print('idxStart = ',idxStart)
    print('idxStop = ',idxStop)
    #Subset Data
    pdData = pdData.iloc[idxStart:idxStop+1]
    #Reset data indices
    pdData = pdData.reset_index(drop=True)
    #Remove unnecessary variables
    posDict = {'xLow':pdData.xLow, 'yLow':pdData.yLow, 'xUp':pdData.xUp, 
               'yUp':pdData.yUp, 'time':pdData.time}
    posData = pd.DataFrame(data=posDict)
    posData['tau'] = posData.time/0.1
    posData.tau = posData.tau - 10.0*(ssTime + 0.075)
    #Calculate the velocity given position
    velData = CalcVelocity(posData)
    #Plot Velocity to make sure it is correct
    velData['tau'] = velData.time/0.1
    velData.tau = velData.tau - 10.0*(ssTime + 0.075)
    
    return (posData, velData)

def StoreLowReData(dirName,fileName):
    global ssTime, DT
    DT = 1.0e-4
    if(fileName == 'pd_Re0_01.txt'):
        ssTime = 0.5
    else:
        ssTime = 1.0
    #Load position and force data
    pdData = []
    print('dirName = ',dirName)
    pdData = pd.read_csv(dirName+fileName,delimiter = ' ')
    #Limit the amount of data to 2 periods worth
    #Start and Stop Indices
    idxStart = int(ssTime/DT)
    if(fileName == 'pd_Re0_01.txt'):
        idxStop = int((ssTime+0.2)/DT)
    else:
        idxStop = int((ssTime+0.4)/DT)
    print('idxStart = ',idxStart)
    print('idxStop = ',idxStop)
    #Subset Data
    pdData = pdData.iloc[idxStart:idxStop+1]
    #Reset data indices
    pdData = pdData.reset_index(drop=True)
    #Remove unnecessary variables
    posDict = {'xLow':pdData.xLow, 'yLow':pdData.yLow, 'xUp':pdData.xUp, 
               'yUp':pdData.yUp, 'time':pdData.time}
    posData = pd.DataFrame(data=posDict)
    posData['tau'] = posData.time/0.1
    posData.tau = posData.tau - 10.0*(ssTime + 0.075)
    #Calculate the velocity given position
    velData = CalcVelocity(posData)
    #Plot Velocity to make sure it is correct
    velData['tau'] = velData.time/0.1
    velData.tau = velData.tau - 10.0*(ssTime + 0.075)
    
    return (posData, velData)


def CalcVelocity(posData):
    #Calculate the velocity of each individual sphere
    #Caculate the components and calculate the magnitude
    #Variables to use
    nTime = int(len(posData.time))
    #Velocity in the x-direction and y-direction
    posData['vxUp'] = 0.0
    posData['vyUp'] = 0.0
    posData['vxLow'] = 0.0
    posData['vyLow'] = 0.0
    for idxTime in range(nTime):
        if((idxTime >= 1) and (idxTime < nTime-1)):
            #Finite difference: midpoint method
            d2t = posData.time[idxTime+1] - posData.time[idxTime-1]
            #Large Sphere
            posData.at[idxTime,'vxUp'] = (posData.at[idxTime+1,'xUp'] - posData.at[idxTime-1,'xUp'])/d2t
            posData.at[idxTime,'vyUp'] = (posData.at[idxTime+1,'yUp'] - posData.at[idxTime-1,'yUp'])/d2t
            #Small Sphere
            posData.at[idxTime,'vxLow'] = (posData.at[idxTime+1,'xLow'] - posData.at[idxTime-1,'xLow'])/d2t
            posData.at[idxTime,'vyLow'] = (posData.at[idxTime+1,'yLow'] - posData.at[idxTime-1,'yLow'])/d2t
    #Velocity at initial time == vel(1)
    posData.at[0,'vxUp'] = posData.at[1,'vxUp']
    posData.at[0,'vyUp'] = posData.at[1,'vyUp']
    posData.at[0,'vxLow'] = posData.at[1,'vxLow']
    posData.at[0,'vyLow'] = posData.at[1,'vyLow']
    #Velocity at final time == vel(ntime-1)
    posData.at[nTime,'vxUp'] = posData.at[nTime-1,'vxUp']
    posData.at[nTime,'vyUp'] = posData.at[nTime-1,'vyUp']
    posData.at[nTime,'vxLow'] = posData.at[nTime-1,'vxLow']
    posData.at[nTime,'vyLow'] = posData.at[nTime-1,'vyLow']
    #Smooth the Velocity Data!
    posData.vyUp = SmoothVelocity(posData,'vyUp')
    posData.vyLow = SmoothVelocity(posData,'vyLow')
    #Calculate the Magnitude
    posData['vMagUp'] = np.sqrt(posData.vxUp**2.0 + posData.vyUp**2.0)
    posData['vMagLow'] = np.sqrt(posData.vxLow**2.0 + posData.vyLow**2.0)
    
    #Create new Velocity Database
    velDict = {'vxUp':posData.vxUp,'vyUp':posData.vyUp,'vxLow':posData.vxLow,
               'vyLow':posData.vyLow,'vMagUp':posData.vMagUp,
               'vMagLow':posData.vMagLow,'time':posData.time}
    velData = pd.DataFrame(data=velDict)
    
    return velData

def SmoothVelocity(data,varName):
    #After plotting the raw data for forces, the term mdu/dt was very noisy
    #We will apply a smoothing algorithm to all velocities in the y-direction
        
    #Use Savitsky-Golay filter from scipy
    yval = data[varName].tolist()
    varValue = savgol_filter(yval, 51, 3)
    
    return varValue

def main():
    global allPosData, allVelData, posData, velData
    #Obtain current directory
    cwd_PYTHON = os.getcwd()
    #Change to data directories
    os.chdir(cwd_PYTHON+'/../CoMData')
    #1) Save data directory
    cwd_DATA = os.getcwd()
    
    #Create Figures for Plotting Position and Velocity Data
    #Displacement and Velocity
    #In this order: D, V
    nRows = 1
    nCols = 2
    #DandV
    fig, ax = plt.subplots(nrows=nRows, ncols=nCols, num=0, 
                             figsize=(10,4),dpi=200)
    
    #Create Subplot titles for all plots
    #ax[0].set_title('Displacement: A = 0.18m',fontsize=20)
    #ax[1].set_title('Velocity: A = 0.18m',fontsize=20)
        
    #Change to RSL#.# directory
    os.chdir(cwd_DATA+'/RSL1.3')
    #Save directory of velocity data
    cwd_VEL = os.getcwd()
    #Store data from posRSL1.3Shifted.csv file
    #Store data from velRSL1.3Shifted.csv file
    allPosData = pd.read_csv('posRSL1.3Shifted.csv',delimiter = ',')
    allVelData = pd.read_csv('velRSL1.3Shifted.csv',delimiter = ',')
    #keep data of importance
    posDataDict = {'amp':allPosData.amp,'Re':allPosData.Re,'time':allPosData.time,
                   'tau':allPosData.tau,'Ly':allPosData.upCOMy,'Sy':allPosData.lowCOMy,
                   'COM':allPosData.COM}
    posData = pd.DataFrame(data=posDataDict)
    velDataDict = {'amp':allVelData.amp,'Re':allVelData.Re,'tau':allVelData.tau,
                   'VelLarge':allVelData.VelLarge,'VelSmall':allVelData.VelSmall,
                   'VelCOM':allVelData.VelCOM}
    velData = pd.DataFrame(data=velDataDict)
    #Plot the displacement over time for each part of the spherobot
    #Make a dataframe of plot parameters
    plotName = 'COM'
    plotMin = [-0.025,-1.875]
    plotMax = [0.07,1.875]
    plotDict = {'Min':plotMin, 'Max':plotMax}
    plotParams = pd.DataFrame(data=plotDict)
    
    #Plot all RSL on same plot
    plotType = 'figure'
    cwd_PLOT = PlotDataAll(ax,posData,velData,plotParams)
        
    #Change to appropriate directory
    #Tight Layout figure
    #Save figure in that directory
    os.chdir(cwd_PLOT+'/')
    fig.tight_layout()
    #fig.savefig(plotName+'_'+plotType+'GR-Smooth.png')
    fig.savefig('v2_COM_SW.png')
    fig.clf()
    plt.close()
    
    return

#---------------------__END_MAIN-------------------------------------
main()
