#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 14:45:42 2018

@author: thomas
"""

import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats

#CONSTANTS
radiusLarge = 0.30 #Radius of large sphere
radiusSmall = 0.15 #Radius of small sphere
RSL = 0.75 #Default resting spring length
FREQ = 10.0 #frequency of oscillation
AMPFRAC = 0.8 #Fraction of max amplitude
MAXAMP = 0.30 #maximum amplitude of oscillation
DENSITY = 2.0 #Fluid density (kg/m-3)
L = 0.3 #Length Scale (Currently just radius of large sphere)
ssFrac = 2.0
ALPHA = radiusLarge/radiusSmall

#LISTS
rslFracList = [1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.0]
ampFracList = [0.4,0.5,0.6,0.7,0.8]
ampSelectList = [0.5,0.7]
rslSelectList = [1.2,1.6]

def PlotFigures(data,xlabel,ylabel,Name):
    #This function solely plots V vs Nu for all data. It plots the following figures
    #1) V/fR vs Re Cumulative (both reg and inset semi-log-x)
    #2) V/fR vs Re select 2 amplitudes (both reg and inset semi-log-x)  (4x4)
    #2) V/fR vs Re select 2 RSL (both reg and inset semi-log-x)  (4x4)
    
    #Color Palette for plotting
    R1=[255/255,255/255,153/255,153/255,204/255]
    G1=[153/255,204/255,255/255,204/255,153/255]
    B1=[204/255,153/255,153/255,255/255,255/255]
    
    #Figure Number
    figNum = 1
    
    #1) V/fR vs Re Cumulative
    fig = plt.figure(num=figNum, figsize=(8,8),dpi=200)
    #axReg = fig.add_subplot(111)
    axReg = fig.add_subplot(211)
    #axReg.set_xlabel(r'$\frac{A_{r}^2}{\delta^2}$',fontsize=24)
    axReg.set_ylabel(r'$\frac{V}{f L}$',fontsize=24)
    #inset location
    #left, bottom, width, height = [0.45,0.17,0.5,0.25] #WholeSweep
    #left, bottom, width, height = [0.45,0.2,0.5,0.25] #Re_s
    #axSemiLog = fig.add_axes([left,bottom,width,height])
    #Semi-log plot
    axSemiLog = fig.add_subplot(212)
    axSemiLog.set_xlabel(r'$\frac{A_{r}^2}{\delta^2}$',fontsize=24)
    axSemiLog.set_ylabel(r'$\frac{V}{f L}$',fontsize=24)
    
    j = 0 #used to select the line color in the loop
    #Loop over every amplitude fraction
    for amp in ampFracList:
        i = 1 #used to scale line color so there is a gradient as rsl changes
        #Loop over every rsl fraction
        for rsl in rslFracList:
            #Select RGB color
            R=R1[j]*i
            G=G1[j]*i
            B=B1[j]*i
            #Plot each (ampFrac,rslFrac) line
            axReg.plot(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                                data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                                color=(R,G,B),linewidth=2)
            axReg.scatter(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                               data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                               color=(R,G,B),s=18)
            axSemiLog.semilogx(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                                data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                                color=(R,G,B),linewidth=2)
            axSemiLog.scatter(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                               data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                               color=(R,G,B),s=18)
            i -= 0.66/(len(rslFracList)+1)
        j += 1
    #label='amp='+str(amp)+' rsl='+str(rsl),color=(R,G,B)
    #lgd = axReg.legend(loc=2, bbox_to_anchor=(1.05,1),borderaxespad=0,ncol=2,fontsize='x-small')
    #Finish up plot and save
    
    #Axes for inset
    axReg.plot([np.amin(data.xval),np.amax(data.xval)],[0,0],'k',linewidth=3) #plot y=0
    #axReg.axis([0.0,np.amax(data.xval)-50,np.amin(data.yval)-0.1,np.amax(data.yval)]) #WholeSweep
    #axReg.axis([0.0,np.amax(data.xval),np.amin(data.yval)-0.1,np.amax(data.yval)])
    axReg.axis([0.1,150.0,np.amin(data.yval)-0.01,np.amax(data.yval)+0.01]) #Re_s
    axSemiLog.semilogx([np.amin(data.xval),np.amax(data.xval)],[0,0],'k',linewidth=3) #plot y=0
    #axSemiLog.axis([np.amin(data.xval),150.0,np.amin(data.yval)-0.01,0.05]) #Re_s
    #axSemiLog.axis([np.amin(data.xval),np.amax(data.xval)-50,np.amin(data.yval)-0.05,0.05])    
    
    #Axes for 2Plot
    #axReg.axis([0.0,np.amax(data.xval)-50,np.amin(data.yval)-0.05,np.amax(data.yval)]) #WholeSweep
    axReg.axis([0.0,150.0,np.amin(data.yval)-0.01,np.amax(data.yval)+0.01]) #Re_s
    #axReg.axis([0.0,150.0,np.amin(data.yval)-0.01,np.amax(data.yval)+0.01])
    #axSemiLog.axis([np.amin(data.xval),np.amax(data.xval)-50,np.amin(data.yval)-0.05,0.05]) #WholeSweep
    axSemiLog.axis([np.amin(data.xval),150.0,np.amin(data.yval)-0.01,np.amax(data.yval)+0.01]) #Re_s
    #axSemiLog.axis([np.amin(data.xval),150.0,np.amin(data.yval)-0.01,0.05])
    
    axReg.tick_params(labelsize=16)
    axSemiLog.tick_params(labelsize=16)
    matplotlib.style.use('default')
    fig.tight_layout()
    fig.savefig('CollapseTransition'+str(Name)+'.svg')#,bbox_extra_artists=(lgd,),bbox_inches='tight')
    fig.clf()
    figNum += 1
    
    '''#1) V/fR vs Re Cumulative (Top = Regime I: Bottom = Regime II)
    fig = plt.figure(num=figNum, figsize=(8,8),dpi=120)
    #axReg = fig.add_subplot(111)
    axRegI = fig.add_subplot(211)
    axRegI.set_title('Regime I',fontsize=18)
    #axRegI.set_xlabel(r'$\frac{\alpha A r}{\delta^2}$',fontsize=16)
    axRegI.set_ylabel(r'$\frac{V}{f \alpha r}$',fontsize=16)
    #inset location
    #left, bottom, width, height = [0.4,0.15,0.5,0.25]
    #axSemiLog = fig.add_axes([left,bottom,width,height])
    #Semi-log plot
    axRegII = fig.add_subplot(212)
    axRegII.set_title('Regime II',fontsize=18)
    axRegII.set_xlabel(r'$\frac{\alpha A r}{\delta^2}$',fontsize=16)
    axRegII.set_ylabel(r'$\frac{V}{f \alpha r}$',fontsize=16)
    
    j = 0 #used to select the line color in the loop
    #Loop over every amplitude fraction
    for amp in ampFracList:
        i = 1 #used to scale line color so there is a gradient as rsl changes
        #Loop over every rsl fraction
        for rsl in rslFracList:
            #Select RGB color
            R=R1[j]*i
            G=G1[j]*i
            B=B1[j]*i
            #Plot each (ampFrac,rslFrac) line
            axRegI.plot(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                                data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                                color=(R,G,B),linewidth=1)
            axRegI.scatter(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                               data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                               color=(R,G,B),s=9)
            axRegII.plot(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                                data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                                color=(R,G,B),linewidth=1)
            axRegII.scatter(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                               data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                               color=(R,G,B),s=9)
            i -= 0.66/(len(rslFracList)+1)
        j += 1
    #label='amp='+str(amp)+' rsl='+str(rsl),color=(R,G,B)
    #lgd = axReg.legend(loc=2, bbox_to_anchor=(1.05,1),borderaxespad=0,ncol=2,fontsize='x-small')
    #Finish up plot and save
    #Axes for 2Regimes
    axRegI.plot([0.0,np.amax(data.xval)],[0,0],'k',linewidth=3) #plot y=0
    axRegII.plot([0.0,np.amax(data.xval)],[0,0],'k',linewidth=3) #plot y=0
    axRegI.axis([0.0,50.0,np.amin(data.yval)-0.05,0.01])
    axRegII.axis([50.0,np.amax(data.xval)-50,-0.01,np.amax(data.yval)])
    matplotlib.style.use('default')
    fig.tight_layout()
    fig.savefig('CollapseTransition'+str(Name)+'.png')#,bbox_extra_artists=(lgd,),bbox_inches='tight')
    fig.clf()
    figNum += 1'''
    
    '''#2) V/fR vs Re select 2 amplitude
    fig = plt.figure(num=figNum, figsize =(16,16),dpi=120)
    #amp = 0.5
    j = 1 #used to select the line color in the loop
    #Loop over every amplitude fraction
    idx = 0
    for amp in ampSelectList:
        i = 1 #used to scale line color so there is a gradient as rsl changes
        #Subplot for Changing RSL Constant Amplitude
        axAmp = fig.add_subplot(2,2,idx+1)
        axAmp.set_title('amp = '+str(amp),fontsize=24)
        #axAmp.set_xlabel(r'$\frac{\alpha A r}{\delta^2}$',fontsize=24)
        axAmp.set_xlabel('',fontsize=24)
        axAmp.set_ylabel(r'$\frac{V}{f L}$',fontsize=24)
        #Loop over every rsl fraction
        for rsl in rslFracList:
            #Subplot for Changing Amplitude Constant RSL
            #Select RGB color
            R=R1[j]*i
            G=G1[j]*i
            B=B1[j]*i
            #Plot each (ampFrac,rslFrac) line
            axAmp.plot(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                            data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                            label='rsl='+str(rsl),color=(R,G,B))
            axAmp.scatter(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                            data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                            label=None,color=(R,G,B))
            axAmp.legend(loc='upper left', fontsize=16)
            axAmp.plot([min(data.xval),max(data.xval)],[0,0],'k') #plot y=0
            #axAmp.axis([0.0,np.amax(data.xval)-50,np.amin(data.yval)-0.05,np.amax(data.yval)]) #WholeSweep
            axAmp.axis([0.0,150.0,np.amin(data.yval),np.amax(data.yval)+0.01])
            i -= 0.66/(len(rslFracList)+1)
        #amp = 0.7
        j += 2
        idx += 1
        axAmp.tick_params(labelsize=16)
    
    #3) V vs Nu each rslFrac
    #fig = plt.figure(num=figNum, figsize = (45,8),dpi=160)
    j = 0 #used to select the line color in the loop
    #Loop over every amplitude fraction
    for amp in ampFracList:
        #rslFrac = 1.2
        i = 1 -  1.32/(len(rslFracList)+1)#used to scale line color so there is a gradient as rsl changes
        #Loop over every rsl fraction
        idx = 0 #used to keep track of which subplot the data will be located
        for rsl in rslSelectList:
            #Create Subplots
            axRSL = fig.add_subplot(2,2,idx+3)
            axRSL.set_title('rsl = '+str(rsl),fontsize=24)
            axRSL.set_xlabel(r'$\frac{A_{r}^2}{\delta^2}$',fontsize=24)
            axRSL.set_ylabel(r'$\frac{V}{f L}$',fontsize=24)
            #Select RGB color
            R=R1[j]*i
            G=G1[j]*i
            B=B1[j]*i
            #Plot each (ampFrac,rslFrac) line
            axRSL.plot(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                            data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                            label='amp='+str(amp),color=(R,G,B))
            axRSL.scatter(data[(data.ampFrac == amp) & (data.rslFrac == rsl)].xval,
                            data[(data.ampFrac == amp) & (data.rslFrac == rsl)].yval,
                            label=None,color=(R,G,B))
            axRSL.legend(loc='upper left', fontsize=16)
            axRSL.plot([min(data.xval),max(data.xval)],[0,0],'k') #plot y=0
            #axRSL.axis([0.0,np.amax(data.xval)-50,np.amin(data.yval)-0.05,np.amax(data.yval)]) #WholeSweep
            axRSL.axis([0.0,150.0,np.amin(data.yval),np.amax(data.yval)+0.01])
            #rslFrac = 1.6
            i -= 0.66*4/(len(rslFracList)+1)
            idx += 1
            axRSL.tick_params(labelsize=16)
        j += 1
    fig.tight_layout()
    fig.savefig('collapseSelect'+str(Name)+'.png')
    fig.clf()
    figNum += 1'''
    
    return

def main():
    
    #Read and import WholeSweep data file
    #Columns used by data.'NameofColumnHere'
    data = pd.read_csv('../FinalSweepData.txt',delimiter = ' ')
    data = data.sort_values(by=['visc','ampFrac','rslFrac'])
    
    #Remove data (Would go here below)
    #BAD DATA POINTS
    
    #RSL1.0: amp=0.8, visc <= 0.05
    #RSL>=1.1: amp=0.8, visc <= 0.02
    #RSL2.0: amp=0.7, visc = 0.01

    #Tom Data Removal
    data = data.ix[~((data.rslFrac == 1.0) & (data.ampFrac == 0.8) & (data.visc <= 0.05))]
    data = data.ix[~((data.rslFrac >= 1.1) & (data.ampFrac == 0.8) & (data.visc <= 0.02))]
    data = data.ix[~((data.ampFrac == 0.7) & (data.rslFrac >= 1.6) & (data.visc <= 0.01))]
    
    #Add new column to data by the following: data['NewColumnName'] = Expression
    data['rsl'] = data.rslFrac*RSL
    data['amp'] = data.ampFrac*MAXAMP
    data['visc'] = data.visc/DENSITY #Converts mu to nu 
    data['cd'] = data.rsl - data.amp - radiusLarge - radiusSmall #Closest distance the 2 spheres get to each other
    
    #Create dimensionless parameters
    data['rsldivamp'] = data.rsl/data.amp
    data['delta'] = np.sqrt(data.visc/(2.0*np.pi*FREQ))
    data['VdivFrsl'] = data.vel/(FREQ*data.rslFrac)
    data['VdivFamp'] = data.vel/(FREQ*data.ampFrac)
    data['VdivFdelta'] = data.vel/(FREQ*data.delta)
    data['VdivF'] = data.vel/FREQ
    data['ampdivdelta'] = data.ampFrac/data.delta
    data['rsldivdelta'] = data.rslFrac/data.delta #Prob not gonna be used
    data['surf'] = data.rsl - radiusSmall - radiusLarge
    data['length'] = data.rsl + radiusSmall + radiusLarge 
    data['Re'] = data.amp*0.8*radiusSmall/data.delta**2
    data['Re_s'] = (data.amp*0.8)**2/data.delta**2
    
    #Viscosity Constraints
    #data = data[(data.visc >= 0.01) & (data.visc < 2000.08)]
    
    #Velocity Constraints
    #data = data[data.vel < 0.0]
    
    #Re Constraints
    #maxRe = 250 #WholeSweep
    maxRe = 250
    #data = data[data.Re < maxRe]
    
    #Amp Constraints
    data = data[data.ampFrac >= 0.4]
    
    #RSL Constraints
    #data = data[data.rslFrac <= 1.6]
    
    print(data)
    print('minRe = %.4f\tmaxRe = %.4f'%(np.amin(data.Re),np.amax(data.Re)))
    
    #The next two lines are the x and y axis values for plotting 
    #(be sure to include "data" in data.variable when making changes)
    
    #Axes to use
    #x ='data.Re'
    x = 'data.Re_s'
    y ='data.vel/(FREQ*data.length)'
    
    #Create new dataframe columns with the evaulation of x and y values
    data['xval'] = eval(x)
    data['yval'] = eval(y)

    PlotName = '2Plots-Re_s'

    PlotFigures(data,x,y,PlotName)


#------------------__END MAIN__-----------------------------------
main()