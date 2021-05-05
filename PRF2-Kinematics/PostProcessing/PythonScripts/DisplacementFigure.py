#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 5 1:36:09 2019

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
listRSL = ['1.0','1.1','1.2','1.3','1.4','1.5','1.6','1.8','2.0']
csfont = {'fontname':'Times New Roman'}

figList, axList = [None]*2, [None]*2

data, ampData, ReData = [None], [None]*5, [None]*20
ReListAll = [None]*5

def MakeDirectory(directory,seed):
    if not os.path.exists(directory+'/'+str(seed)):
        os.makedirs(directory+'/'+str(seed))
    return

def PlotDataAll(ax, data, plotParams, plotType, idxRSL):
    global fignum,ampData, listRSL
    #Change directory to BigPlots directory
    cwd_CURRENT = os.getcwd()
    os.chdir(cwd_CURRENT+'/../../PaperFigures/Results/')
    cwd_FIG = os.getcwd()
    varName = plotParams.varName
    print('varName = ',varName)
    #MakeDirectory(cwd_FIG,varName)
    #cwd_DATA = cwd_FIG+'/'+varName+'/'
    #os.chdir(cwd_DATA)
    #print('cwd_DATA = ',cwd_DATA)
    #Color Palette for plotting
    R1=[255/255,255/255,153/255,153/255,204/255]
    G1=[153/255,204/255,255/255,204/255,153/255]
    B1=[204/255,153/255,153/255,255/255,255/255]
    idxLC = 0
    idxCol = 0
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
        #idxCol = 0
        if(varName == 'Disp'):
            ax[0,0].set_ylabel(r'$\Delta y_R$ (m)',fontsize=24,**csfont)
            ax[1,0].set_ylabel(r'$\Delta y_r$ (m)',fontsize=24,**csfont)
        if(varName == 'Vel'):
            ax[0,0].set_ylabel(r'$v_R$ (m/s)',fontsize=24,**csfont)
            ax[1,0].set_ylabel(r'$v_r$ (m/s)',fontsize=24,**csfont)
        ax[1,idxCol].set_xlabel(r'time $\tau$',fontsize=24,**csfont)
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
        #Loop over Re values
        ReCount = 0
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
            #Plot 1) horizontal axis, 2) disp or vel v time
            ax[0,idxCol].plot([0.0,1.0],[0.0,0.0],'k',linewidth=2)
            ax[1,idxCol].plot([0.0,1.0],[0.0,0.0],'k',linewidth=2)
            if(varName == 'Disp'):
                '''ax[0,idxCol].plot(tempReData.tau,tempReData.Ly - tempReData.loc[0,'Ly'] - (tempReData.COM - tempReData.loc[0,'COM']),
                  color=(R,G,B), linewidth = 1)
                ax[1,idxCol].plot(tempReData.tau,tempReData.Sy - tempReData.loc[0,'Sy']- (tempReData.COM - tempReData.loc[0,'COM']),
                  color=(R,G,B), linewidth = 1)'''
                ax[0,idxCol].plot(tempReData.tau,tempReData.Ly - tempReData.loc[0,'Ly'],
                  color=(R,G,B), linewidth = 2)
                ax[1,idxCol].plot(tempReData.tau,tempReData.Sy - tempReData.loc[0,'Sy'],
                  color=(R,G,B), linewidth = 2)
            if(varName == 'Vel'):
                '''ax[0,idxCol].plot(tempReData.tau,tempReData.VelLarge - tempReData.VelCOM,
                  color=(R,G,B), linewidth = 1)
                ax[1,idxCol].plot(tempReData.tau,tempReData.VelSmall - tempReData.VelCOM,
                  color=(R,G,B), linewidth = 1)'''
                ax[0,idxCol].plot(tempReData.tau,tempReData.VelLarge,
                  color=(R,G,B), linewidth = 2)
                ax[1,idxCol].plot(tempReData.tau,tempReData.VelSmall,
                  color=(R,G,B), linewidth = 2)
            #Add Vertical Dashed Lines every quarter period
            for x in range(1,4):
                xvalue = 0.25*x
                ax[0,idxCol].plot([xvalue,xvalue],[plotParams['RMin'],plotParams['RMax']],color='k',linestyle='--')
                ax[1,idxCol].plot([xvalue,xvalue],[plotParams['rMin'],plotParams['rMax']],color='k',linestyle='--')
            #Decrease ShadeValue
            #ShadeValue -= 0.1
            ShadeValue -= 0.025
            ReCount += 1
        print('ReCount = ',ReCount)
        #Axes parameters
        #Ly
        ax[0,idxCol].tick_params(which='major',axis='both',direction='in',length=8,width=2,labelsize=18)
        ax[0,idxCol].tick_params(which='minor',axis='both',direction='in',length=6,width=1.5)
        ax[0,idxCol].set_xticks(np.arange(0.0,1.1,step=0.2))
        ax[0,idxCol].xaxis.set_minor_locator(MultipleLocator(0.1))
        ax[0,idxCol].set_yticks(np.arange(0.0,0.25,step=0.1))
        ax[0,idxCol].yaxis.set_minor_locator(MultipleLocator(0.05))
        #Sy
        ax[1,idxCol].tick_params(which='major',axis='both',direction='in',length=8,width=2,labelsize=18)
        ax[1,idxCol].tick_params(which='minor',axis='both',direction='in',length=6,width=1.5)
        ax[1,idxCol].set_xticks(np.arange(0.0,1.1,step=0.2))
        ax[1,idxCol].xaxis.set_minor_locator(MultipleLocator(0.1))
        ax[1,idxCol].set_yticks(np.arange(-0.6,0.25,step=0.2))
        ax[1,idxCol].yaxis.set_minor_locator(MultipleLocator(0.1))
        
        ax[0,idxCol].set_xlim(0.0,1.0)
        ax[1,idxCol].set_xlim(0.0,1.0)
        ax[0,idxCol].set_ylim(plotParams['RMin'],plotParams['RMax'])
        ax[1,idxCol].set_ylim(plotParams['rMin'],plotParams['rMax'])
        #Increment line color to match amplitude
        idxLC += 1
        idxCol += 1
    os.chdir(cwd_CURRENT)
    print('cwd = ',os.getcwd())
    return cwd_FIG

def main():
    global allPosData, allVelData, posData, velData
    #Obtain current directory
    cwd_PYTHON = os.getcwd()
    #Change to data directories
    os.chdir(cwd_PYTHON+'/../CoMData')
    #1) Save data directory
    cwd_DATA = os.getcwd()
    #Specify RSL Directories and values to be analyzed
    #All directories
    #dirsRSL = ['RSL1.0','RSL1.1','RSL1.2','RSL1.3','RSL1.4','RSL1.5','RSL1.6','RSL1.8','RSL2.0']
    #1 directory
    dirsRSL = ['RSL1.3']
    
    #Create Figures for Plotting Position and Velocity Data
    #Displacement and Velocity
    #In this order: D, V
    nRows = 2
    nCols = 5
    #D
    figD, axD = plt.subplots(nrows=nRows, ncols=nCols, num=0, 
                             gridspec_kw = {'height_ratios': [1.0,2.5]}, 
                             figsize=(20,10.5),dpi=200)
    #V
    figV, axV = plt.subplots(nrows=nRows, ncols=nCols, num=1,
                               gridspec_kw = {'height_ratios': [1.0,2.0]},
                               figsize=(20,9),dpi=200)
    
    #Save figures in a list for future use
    figList[0] = figD
    figList[1] = figV
    
    #Create Subplot titles for all plots
    for idx in range(len(ampData)):
        ampValue = 0.12 + 0.03*idx
        axD[0,idx].set_title('$\epsilon$ = %.1f'%(ampValue/RSMALL),fontsize=24,**csfont)
        #axD[1,idx].set_title('r Disp: Amp = %.2fm'%ampValue,fontsize=18)
        axV[0,idx].set_title('$/epsilon$ = %.1f'%(ampValue/RSMALL),fontsize=24,**csfont)
        #axV[1,idx].set_title('r Vel: Amp = %.2fm'%ampValue,fontsize=18)
    
    '''#Add Row Annotation
    pad = 5 # in points
    rows = ['RSL={}'.format(row) for row in ['1.0','1.1','1.2','1.3','1.4',
            '1.5','1.6','1.8','2.0']]
    for idx in range(9):
        for row in rows:
            axLd[idx,0].annotate(row, xy=(0,0.5), xytext=(0,pad),
                                 textcoords='offset points', size='large', ha='right',va='center')'''
    
    #Save axes in a list for future use
    axList[0] = axD
    axList[1] = axV       
    
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
        #Plot the displacement over time for each part of the spherobot
        #Make a dataframe of plot parameters
        #plotMin = [-0.05,-0.5,-0.25,-5.0,-20.0,-1.5]
        #plotMax = [0.2,0.15,0.25,5.0,20.0,1.5]
        #plotDict = {'varName':plotName,'min':plotMin,'max':plotMax}
        plotName = ['Disp','Vel']
        plotRMin = [-0.045,-5.0]
        plotrMin = [-0.45,-15.0]
        plotRMax = [0.175,5.0]
        plotrMax = [0.15,15.0]
        plotDict = {'varName':plotName,'RMin':plotRMin, 'RMax':plotRMax,
                    'rMin':plotrMin,'rMax':plotrMax}
        plotParams = pd.DataFrame(data=plotDict)
        
        #Plot all RSL on same plot
        plotType = 'Final_SI'
        for idx in range(len(plotParams.varName)):
            print('varName index = ',idx)
            if(idx < 1):
                cwd_PLOT = PlotDataAll(axList[idx],posData,plotParams.loc[idx,:],plotType,idxRSL)
            else:
                cwd_PLOT = PlotDataAll(axList[idx],velData,plotParams.loc[idx,:],plotType,idxRSL)
    #Loop over variable names
    #Change to appropriate directory
    #Tight Layout figure
    #Save figure in that directory
    for idx in range(len(plotParams.varName)):
        varName = plotParams.loc[idx,'varName']
        os.chdir(cwd_PLOT+'/')
        figList[idx].tight_layout()
        figList[idx].savefig(varName+'_'+plotType+'.png')
        figList[idx].clf()
    return

#---------------------__END_MAIN-------------------------------------
main()
