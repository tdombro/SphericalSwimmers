#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 14:45:42 2018

@author: thomas
"""

import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)

#CONSTANTS
radiusLarge = 0.30 #Radius of large sphere
radiusSmall = 0.15 #Radius of small sphere
RSL = 0.75 #Default resting spring length
FREQ = 10.0 #frequency of oscillation
AMPFRAC = 0.7 #Fraction of max amplitude (Radius Sweep)
MAXAMP = 0.30 #maximum amplitude of oscillation
DENSITY = 2.0 #Fluid density (kg/m-3)
L = 0.3 #Length Scale (Currently just radius of large sphere)
ssFrac = 2.0
ALPHA = radiusLarge/radiusSmall

#LISTS
csfont = {'fontname':'Times New Roman'}
#plt.rcParams["font.family"] = "Times New Roman"
mpl.rcParams['axes.linewidth'] = 1.5 #set the value globally

def Plot3D(data3D,amp,rsl):
    
    #Figure Number
    figNum = 2
    
    #1) V/fR vs Re Cumulative
    fig = plt.figure(num=figNum, figsize=(8,6),dpi=200)
    axReg = fig.add_subplot(111)
    axReg.set_xlabel(r'$Re$',fontsize=24,**csfont)
    #axReg.set_ylabel(r'$\frac{\langle v \rangle}{Rf}$',fontsize=24,**csfont)
    axReg.set_ylabel(r'$\langle v \rangle/Rf$',fontsize=24,**csfont)
    #Inset location
    left, bottom, width, height = [0.25,0.625,0.3,0.3] #Re_s
    axInset = fig.add_axes([left,bottom,width,height])
    
    #data for 3D
    axReg.plot(data3D.xval,data3D.yval,color='k',linewidth=2)
    axReg.scatter(data3D.xval,data3D.yval,color='k',s=32,linewidth=2)
    #Axes for 3D
    axReg.plot([np.amin(data3D.xval),np.amax(data3D.xval)],[0,0],'k',linewidth=2) #plot y=0
    axReg.axis([0,80.0,-0.05,0.105]) #Both/RSLAmp
    axReg.tick_params(which='major',axis='both',direction='in',length=14,width=1,zorder=10)
    axReg.tick_params(which='minor',axis='both',direction='in',length=8,width=0.7)
    axReg.xaxis.set_major_locator(MultipleLocator(10.0))
    axReg.xaxis.set_minor_locator(MultipleLocator(2.5))
    axReg.yaxis.set_major_locator(MultipleLocator(0.02))
    axReg.yaxis.set_minor_locator(MultipleLocator(0.005))
    
    #Inset (Currently Zoomed-In)
    axInset.semilogx(data3D.xval,
                        data3D.yval,
                        color='k',linewidth=2)
    axInset.scatter(data3D.xval,
                           data3D.yval,
                           color='k',s=32,label=None)

    axInset.semilogx([6e-4,100.0],[0,0],'k',linewidth=2) #plot y=0
    axInset.axis([np.amin(data3D.xval),30.0,-0.05,0.0125]) #Zoomed-In Inset
    axInset.tick_params(which='major',axis='both',direction='in',length=10,zorder=10)
    axInset.tick_params(which='minor',axis='both',direction='in',length=5)
    axInset.yaxis.set_major_locator(MultipleLocator(0.02))
    axInset.yaxis.set_minor_locator(MultipleLocator(0.01))
    
    axReg.tick_params(labelsize=16)
    mpl.style.use('default')
    fig.tight_layout()
    fig.savefig('3DSinglSwimmerVelocity.png')#,bbox_extra_artists=(lgd,),bbox_inches='tight')
    fig.clf()
    plt.close()
    
    return

if __name__ == '__main__':
    
    data3D = pd.read_csv('3DVelData.csv',delimiter = ',')
    data3D = data3D.sort_values(by=['Re','amp','rsl'])
    data3D = data3D.reset_index(drop=True)

    data3D = data3D[data3D.Re <= 90]
    print('='*40+'\n')
    print(data3D)
    print('-'*40)
    
    data3D['VdivFR'] = data3D.vel/(8.0*data3D.Rlarge)
    
    #Axes to use
    x3D = 'data3D.Re'
    y3D = 'data3D.VdivFR'
    
    #Create new dataframe columns with the evaulation of x and y values
    data3D['xval'] = eval(x3D)
    data3D['yval'] = eval(y3D)
    
    Plot3D(data3D,0.8,1.0)

