#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 16:21:57 2018

@author: thomas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 12:39:40 2018

@author: thomas
"""

import pandas as pd
import numpy as np
import os
import re
from pathlib import Path
from shutil import copyfile
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

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

#LISTS
rslFracList = [1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.0]
ampFracList = [0.4,0.5,0.6,0.7,0.8]

def func(x, x0, a, b, c, d):
    Re_crit = 50.0
    
    #REGIME I
    #return -1.0*(x/Re_crit)**a*(1-x/Re_crit)**b
    #return -1.0*(x/a)**b*(1.0-(x/Re_crit))**c
    #return -1.0*(x/a)**b*(1.0 - x/Re_crit)**c  
    return -1.0*a*x**b*(1.0 - x/Re_crit)**c  
    #REGIME II
    #return a*np.log10(x)+b  
    #return a + b*(x-x0) + c*(x-x0)**2 + d*(x-x0)**3
    #WHOLE
    #return -1.0*a*x**b*(1 - x/Re_crit)**(1)

def main():
    
    #Color Palette for plotting
    R1=[255/255,255/255,153/255,153/255,204/255]
    G1=[153/255,204/255,255/255,204/255,153/255]
    B1=[204/255,153/255,153/255,255/255,255/255]
    
    #Read and import WholeSweep data file
    #Columns used by data.'NameofColumnHere'
    data = pd.read_csv('WholeSweepData.txt',delimiter = ' ')
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
    
    #Regime II
    #data = data[(data.visc > 0.01) & (data.visc >= 0.02) & (data.visc <=0.06)]
    #Regime I
    
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
    data['Re'] = data.amp*radiusLarge/data.delta**2
    
    #Re Constraint
    maxRe = 125
    data = data[data.Re < maxRe]
    
    #Velocity Constraint
    data = data[(data.vel < 0.00)]# & (data.vel <= 0.05)]
    
    #Specific RSL and Amp Fractions being Used
    RSLFRAC = 1.3
    AMPFRAC = 0.8
    
    #Axes for Plots
    x ='data.delta**-2*data.amp**1*radiusLarge'
    y = 'data.vel/(FREQ*radiusLarge)'
    
    #Create Figure
    fig = plt.figure(num=0, figsize = (8,8),dpi=120)
    ax = fig.add_subplot(111)
    fig1 = plt.figure(num=1, figsize = (8,8),dpi=120)
    ax1 = fig1.add_subplot(111)
    
    #Whole
    data = data[(data.visc > 0.01)]# & (data.visc < 10.08)]# & (data.visc < 10.0)]
    dataRaw = data #Save raw data for future use
    
    #Calculate Function fit for specific Amplitude
    data = data[data.ampFrac == AMPFRAC]
    rslData = data
    j = 0
    i = 1
    for rsl in rslFracList:
        print('beginning of loop')
        
        #Select RGB color
        R=R1[j]*i
        G=G1[j]*i
        B=B1[j]*i
        
        data = rslData
        print(data)
        data = data[data.rslFrac == rsl]
        print(data)
        
        #The next two lines are the x and y axis values for plotting 
        #(be sure to include "data" in data.variable when making changes)
        
        #Axes to use
        x ='data.delta**-2*data.amp**1*radiusLarge'
        #Create new dataframe columns with the evaulation of x and y values
        data['xval'] = eval(x)
        data['yval'] = eval(y)
        print('\n'+'='*60+'\n')
        print(data.xval)
        print(data.yval)
        print('\n'+'-'*60+'\n')
        
        #Use the function defined in func to find a best fit for 'a' and 'b'
        popt, pcov = curve_fit(func,data.xval,data.yval)#,bounds=([1.0e-3,-100,-100],[1.0e3,100,100]))
        print('\n'+'='*60+'\n')
        print(popt)
        print('\n'+'-'*60+'\n')
        #Plot Data and Fit
        xDataFit = np.linspace(1.0e-5,maxRe,100000)
        #ax.plot(data.xval,data.yval,'g',label='rslFrac = %1.1f'%rsl)
        ax.scatter(data.xval,data.yval,color=(R,G,B),label=None)
        #ax.plot(data.xval,func(data.xval, *popt),'r-',label='fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' %tuple(popt))
        #ax.semilogx(data.xval,func(data.xval, *popt),'r-',label='fit: a=%5.3f, b=%5.3f, c=%5.3f' %tuple(popt))
        #ax.plot(data.xval,func(data.xval, *popt),'r-',label='fit: a=%5.3f, b=%5.3f' %tuple(popt))
        #ax.plot(data.xval,func(data.xval, *popt),'r-',label='fit: x0=%.5f, a=%5.5f, b=%5.5f, c=%.5f, d=%.5f' %tuple(popt))
        ax.plot(xDataFit,func(xDataFit, *popt),color=(R,G,B),label='fit: x0=%.5f, a=%5.5f, b=%5.5f, c=%.5f, d=%.5f' %tuple(popt))
        #ax.scatter(data.xval,func(data.xval, *popt),color='r',label=None)
        lgd = ax.legend(loc=2, bbox_to_anchor=(0,-0.1),borderaxespad=0,ncol=1)
        i -= 0.66/(len(rslFracList)+1)
    #REGIME I
    #ax.set_title(r'amp = 0.24: $\frac{VRSL^{1.5}}{fr}$ vs. $\frac{A r}{\delta^2}$ : $(\frac{x}{a})^b(1-\frac{x}{Re_{crit}})^c$')
    ax.set_title(r'amp = %.1f: $\frac{VRSL^{1.5}}{fr}$ vs. $\frac{A r}{\delta^2}$ : $-ax^b(1-\frac{x}{Re_{crit}})^c$'%AMPFRAC)
    #ax.set_title(r'amp = 0.24: $\frac{V}{fr}$ vs. $\frac{A r}{\delta^2}$ : $(\frac{x}{a})^b(Re_{crit} - x)^c$')
    #REGIME II
    #ax.set_title(r'amp = %.1f: $\frac{V}{fR}$ vs. $\frac{A R}{\delta^2}$ : $alog10(x)+b$'%AMPFRAC)
    ax.set_xlabel(r'log($\frac{A R}{\delta^2}$)')
    ax.set_ylabel(r'$\frac{V}{fr}$')
    fig.tight_layout()
    fig.savefig('FunctionFit/attemptRSL.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
    fig.clf()
    
    #Reset available data
    data = dataRaw
    #Calculate Function Fit for specific RSL
    data = data[data.rslFrac == RSLFRAC]
    ampData = data
    j = 0
    for amp in ampFracList:
        print('beginning of loop')
        data = ampData
        print(data)
        data = data[data.ampFrac == amp]
        print(data)
        
        #Select RGB color
        i=1
        R=R1[j]*i
        G=G1[j]*i
        B=B1[j]*i
        
        #The next two lines are the x and y axis values for plotting 
        #(be sure to include "data" in data.variable when making changes)
        
        #Axes to use
        x ='data.delta**-2*data.amp**1*radiusLarge'
        #Create new dataframe columns with the evaulation of x and y values
        data['xval'] = eval(x)
        data['yval'] = eval(y)
        print('\n'+'='*60+'\n')
        print(data.xval)
        print(data.yval)
        print('\n'+'-'*60+'\n')
        
        #Use the function defined in func to find a best fit for 'a' and 'b'
        popt, pcov = curve_fit(func,data.xval,data.yval)#,bounds=([1.0e-3,-100,-100],[1.0e3,100,100]))
        print('\n'+'='*60+'\n')
        print(popt)
        print('\n'+'-'*60+'\n')
        #Plot Data and Fit
        xDataFit = np.linspace(1.0e-5,maxRe,100000)
        #ax1.plot(data.xval,data.yval,'g',label='rslFrac = %1.1f'%rsl)
        ax1.scatter(data.xval,data.yval,color=(R,G,B),label=None)
        #ax1.plot(data.xval,func(data.xval, *popt),'r-',label='fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' %tuple(popt))
        #ax1.semilogx(data.xval,func(data.xval, *popt),'r-',label='fit: a=%5.3f, b=%5.3f, c=%5.3f' %tuple(popt))
        #ax1.plot(data.xval,func(data.xval, *popt),'r-',label='fit: a=%5.3f, b=%5.3f' %tuple(popt))
        ax1.semilogx(xDataFit,func(xDataFit, *popt),color=(R,G,B),label='fit: x0=%.5f, a=%5.5f, b=%5.5f, c=%.5f, d=%.5f' %tuple(popt))
        #ax1.scatter(data.xval,func(data.xval, *popt),color='r',label=None)
        lgd = ax1.legend(loc=2, bbox_to_anchor=(0,-0.1),borderaxespad=0,ncol=1)
        j+=1
    #REGIME I
    #ax1.set_title(r'amp = 0.24: $\frac{VRSL^{1.5}}{fr}$ vs. $\frac{A r}{\delta^2}$ : $(\frac{x}{a})^b(1-\frac{x}{Re_{crit}})^c$')
    ax1.set_title(r'RSL = %.1f: $\frac{VRSL^{1.5}}{fr}$ vs. $\frac{A r}{\delta^2}$ : $-ax^b(1-\frac{x}{Re_{crit}})^c$'%RSLFRAC)
    #ax1.set_title(r'amp = 0.24: $\frac{V}{fr}$ vs. $\frac{A r}{\delta^2}$ : $(\frac{x}{a})^b(Re_{crit} - x)^c$')
    #REGIME II
    #ax1.set_title(r'rsl = %.1f: $\frac{V}{fR}$ vs. $\frac{A R}{\delta^2}$ : $alog10(x)+b$'%RSLFRAC)
    ax1.set_xlabel(r'log($\frac{A R}{\delta^2}$)')
    ax1.set_ylabel(r'$\frac{V}{fr}$')
    fig1.tight_layout()
    fig1.savefig('FunctionFit/attemptAmp.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
    fig1.clf()
    
    
    return


#------------------__END MAIN__-----------------------------------
main()
    
    
    
    