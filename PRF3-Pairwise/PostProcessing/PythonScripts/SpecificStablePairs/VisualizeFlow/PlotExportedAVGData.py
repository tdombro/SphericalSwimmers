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
config = sys.argv[1]
Re = sys.argv[2]

# constructs a filepath for the plot omages of Re=$Re, config=$config, and field=$field
def plotName(cwd,Re,config,field,idx):
    #field = Vort, Pres, Vel, AvgW, AvgP, AvgU
    strDir = cwd+"../Figures/AVG/{0}/{1}/{2}/".format(config,Re,field)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return strDir+"{0}_{1}.png".format(field,idx)

def AddDiscsToPlot(ax,pos):
    RADIUS_LARGE = 0.002
    RADIUS_SMALL = 0.001
    #Add Discs
    circle1 = Circle((pos.loc[0,'aXU'], pos.loc[0,'aYU']), RADIUS_LARGE, facecolor=(0.0,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle1)
    circle2 = Circle((pos.loc[0,'aXL'], pos.loc[0,'aYL']), RADIUS_SMALL, facecolor=(0.0,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle2)
    circle3 = Circle((pos.loc[0,'bXU'], pos.loc[0,'bYU']), RADIUS_LARGE, facecolor=(0.7,)*3,
                             linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle3)
    circle4 = Circle((pos.loc[0,'bXL'], pos.loc[0,'bYL']), RADIUS_SMALL, facecolor=(0.7,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle4)
    return

def GetAxisBounds(pos):
    global maxWin,minWin
    #Also Use centroids to change window zoom if need be
    maxCenterX = max(pos.loc[0,'aXU'],pos.loc[0,'aXL'],pos.loc[0,'bXU'],pos.loc[0,'bXL'])
    minCenterX = min(pos.loc[0,'aXU'],pos.loc[0,'aXL'],pos.loc[0,'bXU'],pos.loc[0,'bXL'])
    maxCenterY = max(pos.loc[0,'aYU'],pos.loc[0,'aYL'],pos.loc[0,'bYU'],pos.loc[0,'bYL'])
    minCenterY = min(pos.loc[0,'aYU'],pos.loc[0,'aYL'],pos.loc[0,'bYU'],pos.loc[0,'bYL'])
    maxSwim = max(abs(maxCenterX),abs(minCenterX),abs(maxCenterY),abs(minCenterY))
    #Check if swimmer is outside of window
    if(maxSwim + 0.01 > maxWin): #Swimmer is leaving window area. Make it larger
        maxWin += 0.002
        minWin = -1.0*maxWin
        print("\n\nmaxWin = %.3f\n\n"%maxWin)
    
    return [minWin,maxWin,minWin,maxWin]
    
def PlotVortField(cwd,time,mx,my,w,pos,pars):
    Re = pars[0]
    config = pars[1]
    field = pars[2]
    idx = pars[3]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200)
    ax.set_title(r'time = %.2fs: $\omega$ field'%(time),fontsize=16)
    #Contour
    #levels = MaxNLocator(nbins=21).tick_values(w.min(), w.max())
    levels = MaxNLocator(nbins=21).tick_values(-5, 5)
    cf = ax.contourf(mx,my,w,cmap='bwr',levels=levels,extend='both')
    fig.colorbar(cf,ax=ax,shrink=0.75)
    #Add Discs
    AddDiscsToPlot(ax,pos)
    axis = GetAxisBounds(pos)
    #axis = [pos.loc[0,'aXU']-0.02,pos.loc[0,'aXU']+0.02,pos.loc[0,'aYU']-0.02,pos.loc[0,'aYU']+0.02]
    ax.axis(axis)
    ax.set_aspect('equal')
    fig.tight_layout()
    #fig.savefig('test_avgW.png')
    fig.savefig(plotName(cwd,Re,config,field,idx))
    fig.clf()
    plt.close()
    return

#Plot Streamlines for interpolated velocity fields Ux and Uy
def PlotStreamsU(cwd,time,mx,my,Ux,Uy,pos,pars):
    Re = pars[0]
    config = pars[1]
    field = pars[2]
    idx = pars[3]
    #Here, we will visualize the velocity field on the new coordinate system
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8),dpi=200)
    ax.set_title(r'time = %.2fs: Streamlines'%(time),fontsize=16)
    #Contour
    magU = np.hypot(Ux,Uy)
    #levels = MaxNLocator(nbins=21).tick_values(magU.min(), magU.max())
    #ax.contourf(mx,my,magU,cmap='viridis',levels=levels)
    #Add Discs
    AddDiscsToPlot(ax,pos)
    #Streamlines
    x = np.linspace(-0.05,0.05,1024)
    y = np.linspace(-0.05,0.05,1024)
    lw = 5*magU.T/magU.max()
    #ax.streamplot(x,y,Ux.T,Uy.T,color='white',linewidth=lw,arrowsize=1.0,density=10)
    ax.streamplot(x,y,Ux.T,Uy.T,color='black',linewidth=1.5,arrowsize=1.0,density=3)
    axis = GetAxisBounds(pos)
    #axis = [pos.loc[0,'aXU']-0.02,pos.loc[0,'aXU']+0.02,pos.loc[0,'aYU']-0.02,pos.loc[0,'aYU']+0.02]
    ax.axis(axis)
    fig.tight_layout()
    fig.savefig(plotName(cwd,Re,config,field,idx))
    #fig.savefig('test_streams.png')
    fig.clf()
    plt.close()
    return

def PlotPresField(cwd,time,mx,my,p,pos,pars):
    Re = pars[0]
    config = pars[1]
    field = pars[2]
    idx = pars[3]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200)
    ax.set_title(r'time = %.2fs: Pressure field'%(time),fontsize=16)
    #Contour
    #levels = MaxNLocator(nbins=21).tick_values(p.min(), p.max())
    levels = MaxNLocator(nbins=21).tick_values(-1, 1)
    cf = ax.contourf(mx,my,p,cmap='PRGn',levels=levels,extend='both')
    fig.colorbar(cf,ax=ax,shrink=0.75)
    #Add Discs
    AddDiscsToPlot(ax,pos)
    axis = GetAxisBounds(pos)
    #axis = [pos.loc[0,'aXU']-0.02,pos.loc[0,'aXU']+0.02,pos.loc[0,'aYU']-0.02,pos.loc[0,'aYU']+0.02]
    ax.axis(axis)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(plotName(cwd,Re,config,field,idx))
    #fig.savefig('test_avgP.png')
    fig.clf()
    plt.close()
    return

def PlotVelField(cwd,time,mx,my,Ux,Uy,pos,pars):
    Re = pars[0]
    config = pars[1]
    field = pars[2]
    idx = pars[3]
    magU = np.hypot(Ux,Uy)
    normUx, normUy = Ux/magU, Uy/magU
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200)
    ax.set_title(r'time = %.2fs: Velocity field'%(time),fontsize=16)
    levels = MaxNLocator(nbins=21).tick_values(0.0, 0.0075)
    cf = ax.contourf(mx,my,magU,cmap='viridis',levels=levels,extend='both')
    #Add Discs
    AddDiscsToPlot(ax,pos)
    fig.colorbar(cf,ax=ax,shrink=0.75)
    space = 8 #spacing for quiver plot
    ax.quiver(mx[::space,::space],my[::space,::space],
                        normUx[::space,::space],normUy[::space,::space],
                        color='white',pivot='mid',scale=50,zorder=5)
    axis = GetAxisBounds(pos)
    #axis = [pos.loc[0,'aXU']-0.02,pos.loc[0,'aXU']+0.02,pos.loc[0,'aYU']-0.02,pos.loc[0,'aYU']+0.02]
    ax.axis(axis)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(plotName(cwd,Re,config,field,idx))
    #fig.savefig('test_avgU.png')
    fig.clf()
    plt.close()
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

def PlotCombinedField(cwd,time,mx,my,w,Ux,Uy,pos,pars):
    Re = pars[0]
    field = pars[1]
    idx = pars[2]
    magU = np.hypot(Ux,Uy)
    normUx, normUy = Ux/magU, Uy/magU
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200,num=10)
    ax.set_title(r'time = %.2fs'%(time),fontsize=16)
    #Vorticity Field (Contour)
    levels = MaxNLocator(nbins=21).tick_values(-5, 5)
    ax.contourf(mx,my,w,cmap='bwr',levels=levels,extend='both')
    #Velocity Field
    #Color changes with magnitude
    #opacity = magU/magU.max()
    cm = mpl.cm.binary
    norm.autoscale([0.0,0.003])
    sm = mpl.cm.ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    space = 8 #spacing for quiver plot
    '''
        ax.quiver(mx[::space,::space].flatten(),my[::space,::space].flatten(),
        normUx[::space,::space].flatten(),normUy[::space,::space].flatten(),
        pivot='mid',scale=50,zorder=5,color=cm(norm(magU[::space,::space].flatten())))
        '''
    ax.quiver(mx[::space,::space].flatten(),my[::space,::space].flatten(),
              normUx[::space,::space].flatten(),normUy[::space,::space].flatten(),
              pivot='mid',scale=50,zorder=5,color='black')
    #Add Discs
    AddDiscsToPlot(ax,pos)
    axis = [pos.loc[0,'xU']-0.025,pos.loc[0,'xU']+0.025,pos.loc[0,'yU']-0.025,pos.loc[0,'yU']+0.025]
    ax.axis(axis)
    ax.set_aspect('equal')
    # Turn off tick labels
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    fig.tight_layout()
    ax = set_size(6,6,ax)
    fig.savefig(plotName(cwd,Re,field,idx))
    #fig.savefig('test_avgComb.png')
    fig.clf()
    plt.close()
    return

def PlotCombinedStream(cwd,time,mx,my,w,Ux,Uy,pos,pars):
    Re = pars[0]
    field = pars[1]
    idx = pars[2]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200,num=10)
    #ax.set_title(r'time = %.2fs'%(time),fontsize=16)
    #Vorticity Field (Contour)
    levels = MaxNLocator(nbins=21).tick_values(-5, 5)
    ax.contourf(mx,my,w,cmap='bwr',levels=levels,extend='both')
    #Streamlines
    x = np.linspace(-0.05,0.05,1024)
    y = np.linspace(-0.05,0.05,1024)
    #lw = 5*magU.T/magU.max()
    ax.streamplot(x,y,Ux.T,Uy.T,color='k',linewidth=1.5,arrowsize=1.0,density=1)
    
    #Add Discs
    AddDiscsToPlot(ax,pos)
    xCM = 0.25*(pos.loc[0,'aXU'] + pos.loc[0,'bXU'] + pos.loc[0,'aXL'] + pos.loc[0,'bXL'])
    yCM = 0.25*(pos.loc[0,'aYU'] + pos.loc[0,'bYU'] + pos.loc[0,'aYL'] + pos.loc[0,'bYL'])
    #axis = [pos.loc[0,'xU']-0.025,pos.loc[0,'xU']+0.025,pos.loc[0,'yU']-0.025,pos.loc[0,'yU']+0.025]
    axis = [xCM - 0.025,xCM+0.025,yCM-0.025,yCM+0.025]
    ax.axis(axis)
    ax.set_aspect('equal')
    # Turn off tick labels
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    fig.tight_layout()
    ax = set_size(6,6,ax)
    fig.savefig(plotName(cwd,Re,field,idx))
    #fig.savefig('test_avgComb.png')
    fig.clf()
    plt.close()
    return

def PlotForceField(cwd,time,mx,my,Fx,Fy,pos,pars):
    Re = pars[0]
    field = pars[1]
    idx = pars[2]
    magF = np.hypot(Fx,Fy)
    normFx, normFy = Fx/magF, Fy/magF
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200)
    ax.set_title(r'time = %.2fs: Velocity field'%(time),fontsize=16)
    levels = MaxNLocator(nbins=21).tick_values(0.0, np.amax(magF)/2.0)
    ax.contourf(mx,my,magF,cmap='viridis',levels=levels,extend='both')
    #Add Discs
    AddDiscsToPlot(ax,pos)
    space = 8 #spacing for quiver plot
    ax.quiver(mx[::space,::space],my[::space,::space],
              normFx[::space,::space],normFy[::space,::space],
              color='white',pivot='mid',scale=50,zorder=5)
              axis = [pos.loc[0,'xU']-0.025,pos.loc[0,'xU']+0.025,pos.loc[0,'yU']-0.025,pos.loc[0,'yU']+0.025]
              ax.axis(axis)
              ax.set_aspect('equal')
              # Turn off tick labels
              ax.set_yticklabels([])
              ax.set_xticklabels([])
              fig.tight_layout()
              ax = set_size(6,6,ax)
              fig.savefig(plotName(cwd,Re,field,idx))
              #fig.savefig('test_avgU.png')
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
    pos = data[data['time'] == time]
    pos = pos.reset_index(drop=True)
    return pos

def GetPosDataLength(cwd):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    return len(data['time'])

def GetAvgFieldData(cwd,idx):
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    #cwd = cwd_PYTHON
    fieldData = pd.read_csv(cwd+'AVG_Acc_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    WList  = fieldData['avgW'].values.tolist()
    PList  = fieldData['avgP'].values.tolist()
    UxList = fieldData['avgUx'].values.tolist()
    UyList = fieldData['avgUy'].values.tolist()
    axList = fieldData['avgax'].values.tolist()
    ayList = fieldData['avgay'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 1024, 1024
    mxArr = np.array(mxList).reshape((Nx,Ny))
    myArr = np.array(myList).reshape((Nx,Ny))
    WArr  = np.array(WList).reshape((Nx,Ny))
    PArr  = np.array(PList).reshape((Nx,Ny))
    UxArr = np.array(UxList).reshape((Nx,Ny))
    UyArr = np.array(UyList).reshape((Nx,Ny))
    axArr = np.array(axList).reshape((Nx,Ny))
    ayArr = np.array(ayList).reshape((Nx,Ny))
    return (mxArr, myArr, WArr, PArr, UxArr, UyArr, axArr, ayArr)

if __name__ == '__main__':
    
    #READ ALL AVG FILES IN A SIMULATION DIRECTORY
    #EXTRACT AVERAGE FIELD DATA INTO NUMPY ARRAYS
    #PLOT AVERAGED FIELD DATA
    #Simulation Parameters
    #simList = ['HB','SF','L','V']
    cwd_Re = cwd_PYTHON+'../'+config+'/SweepRe/Re'+Re+'/'
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
        strDir = cwd_PYTHON+"../Figures/AVG/"+config+"/"+Re+"/avgU/"
        #AVGPlot = pathlib.Path(strDir+'avgU_%i.png'%countPer)
        AVGPlot = pathlib.Path(cwd_DATA+'AVG_Acc_%04d.csv'%countPer)
        #if not AVGPlot.exists ():
        if AVGPlot.exists ():    
            start = t.clock()
            
            #Get Avg Field Data
            mx,my,avgW,avgP,avgUx,avgUy,avgax,avgay = GetAvgFieldData(cwd_DATA,countPer)
            #Extract Position and Time Data
            time = np.round(0.05 + countPer*PERIOD,2)
            #print('time = ',time)
            posData = GetPosData(cwd_POS,time)
            #print(posData)
            #Plot Averaged Field Data
            #Vorticity
            pars = [Re,config,'avgW',countPer]
            #PlotVortField(cwd_PYTHON,time,mx,my,avgW,posData,pars)
            #Pressure
            #pars[2] = 'avgP'
            #PlotPresField(cwd_PYTHON,time,mx,my,avgP,posData,pars)
            #Velocity Field
            #pars[2] = 'avgU'
            #PlotVelField(cwd_PYTHON,time,mx,my,avgUx,avgUy,posData,pars)
            #Combination of Vorticity and Velocity field (Quiver)
            #pars[2] = 'combine'
            #PlotCombinedField(cwd_PYTHON,time,mx,my,avgW,avgUx,avgUy,posData,pars)
            #Streamlines                                                                                       
            #pars[2] = 'stream'
            #PlotStreamsU(cwd_PYTHON,time,mx,my,avgUx,avgUy,posData,pars)
            #Force field with quiver
            MASS_FLUID = DX*DX*RHO #mass per length
            avgFx, avgFy = MASS_FLUID*avgax, MASS_FLUID*avgay
            pars[2] = 'forcefield'
            #PlotForceField(cwd,time,mx,my,avgFx,avgFy,posData,pars)
            pars[2] = 'forcestream'
            #PlotCombinedStream(cwd_PYTHON,time,mx,my,avgW,avgFx,avgFy,posData,pars)
    
            stend = t.clock()
            diff = stend - start
            print('Time to run for 1 period = %.5fs'%diff)
            sys.stdout.flush()

