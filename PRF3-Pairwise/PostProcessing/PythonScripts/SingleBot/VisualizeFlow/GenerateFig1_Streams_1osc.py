#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 21:59:46 2020

@author: thomas
"""

import numpy as np
import pandas as pd
import os, sys
import time as t
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
import pathlib
from matplotlib.colors import Normalize
norm = Normalize()

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
RHO = 1000.0
DX = 0.025/256.0
PERIOD = 0.1
maxWin = 0.03
minWin = -1.0*maxWin
#config = sys.argv[1]
Re = sys.argv[1]

# constructs a filepath for the plot omages of Re=$Re, config=$config, and field=$field
def plotName(cwd,Re,field,idx):
    #field = Vort, Pres, Vel, AvgW, AvgP, AvgU
    strDir = cwd+"../Figures/AVG/{0}/{1}/".format(Re,field)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return strDir+"{0}_{1}".format(field,idx)

def AddDiscsToPlot(ax,pos):
    RADIUS_LARGE = 0.002
    RADIUS_SMALL = 0.001
    #Add Discs
    circle1 = Circle((pos.loc[0,'xU'], pos.loc[0,'yU']), RADIUS_LARGE, facecolor=(0.0,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle1)
    circle2 = Circle((pos.loc[0,'xL'], pos.loc[0,'yL']), RADIUS_SMALL, facecolor=(0.0,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle2)
    return

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)
    return ax

def PlotCombinedStream(cwd,time,mx,my,w,Ux,Uy,pos,pars):
    Re = pars[0]
    field = pars[1]
    idx = pars[2]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200,num=10)
    #ax.set_title(r'time = %.2fs'%(time),fontsize=16)
    #Streamlines
    x = np.linspace(-0.05,0.05,1024)
    y = np.linspace(-0.05,0.05,1024)
    #Vorticity Field (Contour)
    if(field == 'W_stream'):
        levels = MaxNLocator(nbins=21).tick_values(-5, 5)
        ax.contourf(mx,my,w,cmap='bwr',levels=levels,extend='both')
        ax.streamplot(x,y,Ux.T,Uy.T,color='k',linewidth=1.5,arrowsize=2.0,density=2.0)
    else:
        magU = np.hypot(Ux,Uy)
        normUx, normUy = Ux/magU, Uy/magU
        levels = MaxNLocator(nbins=21).tick_values(0.0, 0.0075)
        ax.contourf(mx,my,magU,cmap='viridis',levels=levels,extend='both')
        ax.streamplot(x,y,Ux.T,Uy.T,color='w',linewidth=1.5,arrowsize=2.0,density=2.0)
    
    #Add Discs
    AddDiscsToPlot(ax,pos)
    xCM = 0.8*pos.loc[0,'xU'] + 0.2*pos.loc[0,'xL']
    yCM = 0.8*pos.loc[0,'yU'] + 0.2*pos.loc[0,'yL']
    #axis = [pos.loc[0,'xU']-0.025,pos.loc[0,'xU']+0.025,pos.loc[0,'yU']-0.025,pos.loc[0,'yU']+0.025]
    axis = [xCM - 0.025,xCM+0.025,yCM-0.025,yCM+0.025]
    ax.axis(axis)
    ax.set_aspect('equal')
    # Turn off tick labels
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    fig.tight_layout()
    ax = set_size(6,6,ax)
    fig.savefig(plotName(cwd,Re,field,idx)+'.png')
    fig.savefig(plotName(cwd,Re,field,idx)+'.svg')
    #fig.savefig('test_avgComb.png')
    fig.clf()
    plt.close()
    return

# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    #return cwd+"/pd.txt"
    #cwd = cwd_PYTHON
    return cwd+"/pd.txt"

def GetPosData(cwd,time):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    topData = data[data['idx'] == 6].copy()
    botData = data[data['idx'] == 19].copy()
    topData = topData.sort_values(by=['time'])
    botData = botData.sort_values(by=['time'])
    topData = topData.reset_index(drop=True)
    botData = botData.reset_index(drop=True)
    dictPos = {'xU':topData['x'],'yU':topData['y'],'xL':botData['x'],'yL':botData['y'],'time':topData['time']}
    pos = pd.DataFrame(data=dictPos)
    pos = pos[pos['time'] == time]
    pos = pos.reset_index(drop=True)
    return pos

def GetPosDataLength(cwd):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    topData = data[data['idx'] == 6].copy()
    botData = data[data['idx'] == 19].copy()
    topData = topData.sort_values(by=['time'])
    botData = botData.sort_values(by=['time'])
    topData = topData.reset_index(drop=True)
    botData = botData.reset_index(drop=True)
    assert len(topData['time']) == len(botData['time'])
    return len(topData['time'])

def GetAvgFieldData(cwd,idx):
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    #cwd = cwd_PYTHON
    #fieldData = pd.read_csv(cwd+'AVG_Acc_%04d.csv'%idx,delimiter=' ')
    fieldData = pd.read_csv(cwd+'AVG_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    WList  = fieldData['avgW'].values.tolist()
    PList  = fieldData['avgP'].values.tolist()
    UxList = fieldData['avgUx'].values.tolist()
    UyList = fieldData['avgUy'].values.tolist()
    #axList = fieldData['avgax'].values.tolist()
    #ayList = fieldData['avgay'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 1024, 1024
    mxArr = np.array(mxList).reshape((Nx,Ny))
    myArr = np.array(myList).reshape((Nx,Ny))
    WArr  = np.array(WList).reshape((Nx,Ny))
    PArr  = np.array(PList).reshape((Nx,Ny))
    UxArr = np.array(UxList).reshape((Nx,Ny))
    UyArr = np.array(UyList).reshape((Nx,Ny))
    #axArr = np.array(axList).reshape((Nx,Ny))
    #ayArr = np.array(ayList).reshape((Nx,Ny))
    #return (mxArr, myArr, WArr, PArr, UxArr, UyArr, axArr, ayArr)
    return (mxArr, myArr, WArr, PArr, UxArr, UyArr)

if __name__ == '__main__':
    
    #READ ALL AVG FILES IN A SIMULATION DIRECTORY
    #EXTRACT AVERAGE FIELD DATA INTO NUMPY ARRAYS
    #PLOT AVERAGED FIELD DATA
    #Simulation Parameters
    #simList = ['HB','SF','L','V']
    cwd_Re = cwd_PYTHON+'../SweepRe/Re'+Re+'/'
    #Extract Position Data
    cwd_POS = cwd_Re
    #Calculate # Periods
    DUMP_INT = 20.0
    nTime = GetPosDataLength(cwd_POS)
    nPer = int(np.trunc(1.0*nTime/DUMP_INT))
    #nPer = 2
    #Paths to data and plots
    cwd_DATA = cwd_Re+'/VTK/AVG/'
    countPer = 0
    for countPer in range(nPer):
        if(countPer == 100):
            #strDir = cwd_PYTHON+"../Figures/AVG/"+config+"/"+Re+"/avgU/"
            #AVGPlot = pathlib.Path(strDir+'avgU_%i.png'%countPer)
            #AVGPlot = pathlib.Path(cwd_DATA+'AVG_Acc_%04d.csv'%countPer)
            AVGPlot = pathlib.Path(cwd_DATA+'AVG_%04d.csv'%countPer)
            #if not AVGPlot.exists ():
            if AVGPlot.exists ():
                start = t.clock()
                #Get Avg Field Data
                mx,my,avgW,avgP,avgUx,avgUy = GetAvgFieldData(cwd_DATA,countPer)
                #Extract Position and Time Data
                time = np.round(0.05 + countPer*PERIOD,2)
                #print('time = ',time)
                posData = GetPosData(cwd_POS,time)
                #print(posData)
                #Plot Averaged Field Data
                #Vorticity And Streamlines
                pars = [Re,'W_stream',countPer]
                PlotCombinedStream(cwd_PYTHON,time,mx,my,avgW,avgUx,avgUy,posData,pars)
                pars[1] = 'U_stream'
                PlotCombinedStream(cwd_PYTHON,time,mx,my,avgW,avgUx,avgUy,posData,pars)

                stend = t.clock()
                diff = stend - start
                print('Time to run for 1 period = %.5fs'%diff)
                sys.stdout.flush()

