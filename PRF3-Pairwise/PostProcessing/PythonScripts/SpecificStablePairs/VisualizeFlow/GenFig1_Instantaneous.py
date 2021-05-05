#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 17:22:55 2020

@author: thomas
"""

#Import modules
import numpy as np
import pandas as pd
import os, sys
import time as t
import matplotlib as mpl
import pathlib
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import Normalize
norm = Normalize()

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
config = sys.argv[1]
Re = sys.argv[2]
PERIOD = 0.1
maxWin = 0.02
minWin = -1.0*maxWin

#AUXILIARY FUNCTIONS
#Position Data Functions
# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    return cwd+"/pd.txt"

def GetPosData(cwd,time):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    pos = data[data['time'] == time]
    pos = pos.reset_index(drop=True)
    return pos

def GetPosDataLength(cwd):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    return len(data['time'])

### VTK Functions
##### We need to extract a few different things from the VTK file
##### Lines 0-10: Junk
##### Line 11: time value
##### Line 12: Gives the dimensions for coordinates Nx, Ny, Ndim
##### Line 13: Name of variable, #of points, data type
##### Line 14 -> trunc(Nx/9)+1+14: X_coord values
##### For this demonstration, we know Nx = 1024
##### Line 128: var name, #of pts, data type
##### Lines 129 -> 129+trunc(Ny/9)+1: y_coord values
##### For this demonstration, we know Ny = 1024
##### Lines 244-245: Z_coord info = 0
##### Line 246: POINT_DATA, #of pts = Npts
##### Line 247: SCALARS, var Name, data type
##### Line 248: Junk
##### Lines 249 -> 249 + trunc(Npts/9)+1: var values
##### For this demonstration, Npts = 1048576
##### Line 116757: FIELD, FieldData, 2
##### Line 116758: var Name, Ndim, Npts, data type
##### Lines 116759 -> trunc(Npts*Ndim/9)+1+116759: field values
##### For this demonstration: Npts = 1048576
##### Line 233268: var Name, Ndim, data type
##### Line 233269 -> 233269 + trunc(Ndim*Npts/9)+1: field values

#VTK Functions
def ReadVTK(data):
    #Obtain specific values and arrays from read file
    #Time
    timeValue = float(data[11][0])
    print('time = ',timeValue)
    #Obtain Refinement sizes
    Nx = int(data[12][1])
    Ny = int(data[12][2])
    #print('Nx = ',Nx)
    #print('Ny = ',Ny)
    #Obtain xcoord
    Nxline = int(np.trunc(Nx/9.0)+1.0)
    startIdx = 14
    endIdx = startIdx+Nxline
    xList = data[startIdx:endIdx]
    #print(xList[0])
    xFlat = [item for sublist in xList for item in sublist]
    xVals = [float(i) for i in xFlat]
    xArr = np.array(xVals)
    #Obtain ycoord
    Nyline = int(np.trunc(Ny/9.0)+1.0)
    startIdx = endIdx+1
    endIdx = startIdx+Nyline
    yList = data[startIdx:endIdx]
    #print(yList[0])
    yFlat = [item for sublist in yList for item in sublist]
    yVals = [float(i) for i in yFlat]
    yArr = np.array(yVals)
    #Scalar Data
    Npts = int(data[endIdx+2][1])
    Nlines = int(np.trunc(Npts/9.0)+1.0)
    #Obtain Omega
    startIdx = endIdx+5
    endIdx = startIdx+Nlines
    wList = data[startIdx:endIdx]
    wFlat = [item for sublist in wList for item in sublist]
    wVals = [float(i) for i in wFlat]
    wArr = np.array(wVals)
    wArr = wArr.reshape((Nx,Ny))
    #Field Data
    Ndim = int(data[endIdx+1][1])
    Nlines = int(np.trunc(Npts*Ndim/9.0)+1.0)
    #Obtain Pressure
    startIdx = endIdx + 2
    endIdx = startIdx + Nlines
    pList = data[startIdx:endIdx]
    pFlat = [item for sublist in pList for item in sublist]
    pVals = [float(i) for i in pFlat]
    pArr = np.array(pVals)
    pArr = pArr.reshape((Nx,Ny))
    #Obtain Velocity
    Ndim = int(data[endIdx][1])
    Nlines = int(np.trunc(Npts*Ndim/9.0)+1.0)
    startIdx = endIdx + 1
    endIdx = startIdx + Nlines
    uList = data[startIdx:endIdx]
    uFlat = [item for sublist in uList for item in sublist]
    uVals = [float(i) for i in uFlat]
    uArr = np.array(uVals)
    uArr = uArr.reshape((Nx,Ny,Ndim))
    uxArr = uArr[:,:,0]
    uyArr = uArr[:,:,1]
    #print(uArr[0]) 
    
    #Make mesh out of xArr and yArr
    mx,my = np.meshgrid(xArr,yArr,indexing='ij')
    #print(mx[0])

    return (timeValue,mx,my,wArr.T,pArr.T,uxArr.T,uyArr.T)

def ReadAllVTK(cwd,idx):
    count = 0
    with open(cwd+'DATA%05d.vtk'%idx) as fp:
        for line in iter(fp.readline, ''):
            count += 1
    allData = [[] for i in range(count)]
    print(len(allData))
    count = 0
    with open(cwd+'DATA%05d.vtk'%idx) as fp:
        for line in fp:
            allData[count] = line.split()
            count += 1
        print(count)
        return allData

### PLOTTING
#### CONTAINS ALL FUNCTIONS USED FOR PLOTTING

# constructs a filepath for the plot omages of Re=$Re, config=$config, and field=$field
def plotName(cwd,Re,config,field,idx):
    #field = Vort, Pres, Vel, AvgW, AvgP, AvgU
    strDir = cwd+"../Figures/Inst/{0}/{1}/{2}/".format(config,Re,field)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return strDir+"Re{0}_{1}_{2}_{3}".format(Re,config,field,idx)

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
    circle3 = Circle((pos.loc[0,'bXU'], pos.loc[0,'bYU']), RADIUS_LARGE, facecolor=(0.5,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle3)
    circle4 = Circle((pos.loc[0,'bXL'], pos.loc[0,'bYL']), RADIUS_SMALL, facecolor=(0.5,)*3,
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
    config = pars[1]
    field = pars[2]
    idx = pars[3]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200,num=10)
    #ax.set_title(r'time = %.2fs'%(time),fontsize=16)
    #Streamlines
    x = np.linspace(-0.05,0.05,1024)
    y = np.linspace(-0.05,0.05,1024)
    #Vorticity Field (Contour)
    if(field == 'W_Stream'):
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
    fig.savefig(plotName(cwd,Re,config,field,idx)+'.png')
    #fig.savefig(plotName(cwd,Re,config,field,idx)+'.svg')
    #fig.savefig('test_avgComb.png')
    fig.clf()
    plt.close()
    return

def PlotCombinedField(cwd,time,mx,my,w,Ux,Uy,pos,pars):
    Re = pars[0]
    config = pars[1]
    field = pars[2]
    idx = pars[3]
    magU = np.hypot(Ux,Uy)
    normUx, normUy = Ux/magU, Uy/magU
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200)
    #ax.set_title(r'time = %.2fs'%(time),fontsize=16)
    #Vorticity Field (Contour)
    levels = MaxNLocator(nbins=21).tick_values(-5, 5)
    cf = ax.contourf(mx,my,w,cmap='bwr',levels=levels,extend='both')
    #fig.colorbar(cf,ax=ax,shrink=0.75)
    #Velocity Field
    #Color changes with magnitude
    #opacity = magU/magU.max()
    cm = mpl.cm.binary
    norm.autoscale([0.0,0.005])
    sm = mpl.cm.ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    space = 8 #spacing for quiver plot
    ax.quiver(mx[::space,::space].flatten(),my[::space,::space].flatten(),
              normUx[::space,::space].flatten(),normUy[::space,::space].flatten(),
              pivot='mid',scale=50,zorder=5,color=cm(norm(magU[::space,::space].flatten())))
              
    #Add Discs
    AddDiscsToPlot(ax,pos)
    xCM = 0.25*(pos.loc[0,'aXU'] + pos.loc[0,'bXU'] + pos.loc[0,'aXL'] + pos.loc[0,'bXL'])
    yCM = 0.25*(pos.loc[0,'aYU'] + pos.loc[0,'bYU'] + pos.loc[0,'aYL'] + pos.loc[0,'bYL'])
    #axis = [pos.loc[0,'xU']-0.025,pos.loc[0,'xU']+0.025,pos.loc[0,'yU']-0.025,pos.loc[0,'yU']+0.025]
    axis = [xCM - 0.025,xCM+0.025,yCM-0.025,yCM+0.025]
    ax.axis(axis)
    ax.set_aspect('equal')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    fig.tight_layout()
    ax = set_size(6,6,ax)
    fig.savefig(plotName(cwd,Re,config,field,idx)+'.png')
    #fig.savefig('test_avgComb.png')
    fig.clf()
    plt.close()
    return
    
### MAIN SCRIPT
#### READ ALL VTK FILES IN A SIMULATION DIRECTORY
#### CALCULATE AVERAGE FIELD DATA
#### EXPORT AVERAGED FIELD DATA TO CSV

if __name__ == '__main__':
    #Simulation Parameters
    cwd_Re = cwd_PYTHON+'../'+config+'/SweepRe/Re'+Re+'/'
    #Extract Position Data
    cwd_POS = cwd_Re
    #Calculate # Periods
    DUMP_INT = 20.0
    nTime = GetPosDataLength(cwd_POS)
    nPer = int(np.trunc(1.0*nTime/DUMP_INT))
    #nPer = 1
    #Paths to data and plots
    cwd_DATA = cwd_Re+'/VTK/'
    #Test
    #nPer = 1
    for countPer in range(nPer):
        start = t.clock()
        for idx in range(0,int(DUMP_INT)):
            dumpIdx = int(DUMP_INT)*countPer+idx
            #strDir = cwd_PYTHON+"../Figures/Inst/"+config+"/"+Re+"/instU/"
            #InstPlot = pathlib.Path(strDir+'instU_%i.png'%dumpIdx)
            strDir = cwd_DATA
            InstPlot = pathlib.Path(strDir+'DATA%5d.png'%dumpIdx)
            if not InstPlot.exists ():
                #Read All VTK Data into List
                allData = ReadAllVTK(cwd_DATA,dumpIdx)
                #Extract important values from allData
                time, mx, my, wArr, pArr, uxArr, uyArr = ReadVTK(allData)
                #Plot Instantaneous Flow
                posData = GetPosData(cwd_POS,time)
                #print(posData)
                #Plot Instantaneous Fields
                #Combined Stream
                pars = [Re,config,'W_Stream',dumpIdx]
                PlotCombinedStream(cwd_PYTHON,time,mx,my,wArr,uxArr,uyArr,posData,pars)
                #Combined Quiver
                pars[2] = 'W_Quiver'
                PlotCombinedField(cwd_PYTHON,time,mx,my,wArr,uxArr,uyArr,posData,pars)

            stend = t.clock()
            diff = stend - start
            print('Time to run for 1 period = %.5fs'%diff)
            sys.stdout.flush()
        else:
            print('%s: Re = %s: Per = %i: Inst file exists already'%(config,Re,countPer))

