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
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
Par = sys.argv[1]
ParVal = sys.argv[2]

#AUXILIARY FUNCTIONS
#Position Data Functions
# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    return cwd+"/pd.txt"

def GetPosData(cwd,time):
    '''
    data = pd.read_csv(pname(cwd),delimiter=' ')
    pos = data[data['time'] == time]
    pos = pos.reset_index(drop=True)
    '''
    data = np.loadtxt(pname(cwd),delimiter=' ')
    dataDict = {'xU':data[:,0],'yU':data[:,1],
                'xL':data[:,2],'yL':data[:,3],
                'curr':data[:,4],'des':data[:,5],'time':data[:,6]}
    posData = pd.DataFrame(data=dataDict)
    pos = posData[posData['time'] == time]
    pos = pos.reset_index(drop=True)
    return pos

def GetPosDataLength(cwd):
    #data = pd.read_csv(pname(cwd),delimiter=' ')
    data = np.loadtxt(pname(cwd),delimiter=' ')
    #return len(data['time'])
    return len(data[:,6])

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
    Nz = int(data[12][3])
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
    #Obtain zcoord
    Nzline = int(np.trunc(Nz/9.0)+1.0)
    startIdx = endIdx+1
    endIdx = startIdx+Nzline
    zList = data[startIdx:endIdx]
    #print(yList[0])
    zFlat = [item for sublist in zList for item in sublist]
    zVals = [float(i) for i in zFlat]
    zArr = np.array(zVals)
    #Point Data
    Npts = int(data[endIdx][1])
    Nlines = int(np.trunc(Npts*3.0/9.0)+1.0)
    #Obtain Omega (Vector)
    startIdx = endIdx+2
    endIdx = startIdx+Nlines
    wList = data[startIdx:endIdx]
    wFlat = [item for sublist in wList for item in sublist]
    wVals = [float(i) for i in wFlat]
    wArr = np.array(wVals)
    wArr = wArr.reshape((Nx,Ny,Nz,3))
    wxArr = wArr[:,:,:,0]
    wyArr = wArr[:,:,:,1]
    wzArr = wArr[:,:,:,2]
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
    pArr = pArr.reshape((Nx,Ny,Nz))
    #Obtain Velocity
    Ndim = int(data[endIdx][1])
    Nlines = int(np.trunc(Npts*Ndim/9.0)+1.0)
    startIdx = endIdx + 1
    endIdx = startIdx + Nlines
    uList = data[startIdx:endIdx]
    uFlat = [item for sublist in uList for item in sublist]
    uVals = [float(i) for i in uFlat]
    uArr = np.array(uVals)
    uArr = uArr.reshape((Nx,Ny,Nz,Ndim))
    uxArr = uArr[:,:,:,0]
    uyArr = uArr[:,:,:,1]
    uzArr = uArr[:,:,:,2]
    #print(uArr[0]) 
    
    #Make mesh out of xArr and yArr
    mx,my,mz = np.meshgrid(xArr,yArr,zArr,indexing='ij')
    #print(mx[0])

    return (timeValue,mx,my,mz,wxArr.T,wyArr.T,wzArr.T,pArr.T,uxArr.T,uyArr.T,uzArr.T)

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

### AUXILIARY PLOTTING FUNCTIONS
### USED FOR TESTING 1 VTK ONLY
def AddDiscsToPlot(ax,pos):
    RADIUS_LARGE = 0.3
    RADIUS_SMALL = 0.15
    #Add Discs
    circle1 = Circle((pos['xU'], pos['yU']), RADIUS_LARGE, facecolor='k',
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle1)
    circle2 = Circle((pos['xL'], pos['yL']), RADIUS_SMALL, facecolor='k',
                                      linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle2)
    return

def PlotVortField(time,mx,my,w,pos):
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200)
    ax.set_title(r'time = %.2fs: $\omega$ field'%(time),fontsize=16)
    #Contour
    #levels = MaxNLocator(nbins=21).tick_values(w.min(), w.max())
    levels = MaxNLocator(nbins=21).tick_values(-10, 10)
    cf = ax.contourf(mx,my,w,cmap='bwr',levels=levels,extend='both')
    fig.colorbar(cf,ax=ax,shrink=0.75)
    #Add Discs
    AddDiscsToPlot(ax,pos)
    ax.axis([-1.0,1.0,-2.0,1.0])
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'../Figures/testVTK/testVort.png')
    fig.clf()
    plt.close()
    return

def PlotPresField(time,mx,my,p,pos):
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200)
    ax.set_title(r'time = %.2fs: Pressure field'%(time),fontsize=16)
    #Contour
    #levels = MaxNLocator(nbins=21).tick_values(p.min(), p.max())
    levels = MaxNLocator(nbins=21).tick_values(-10, 10)
    cf = ax.contourf(mx,my,p,cmap='PRGn',levels=levels,extend='both')
    fig.colorbar(cf,ax=ax,shrink=0.75)
    #Add Discs
    AddDiscsToPlot(ax,pos)
    ax.axis([-1.0,1.0,-2.0,1.0])
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'../Figures/TestVTK/testPres.png')
    fig.clf()
    plt.close()
    return

def PlotVelField(time,mx,my,Ux,Uy,pos):
    magU = np.hypot(Ux,Uy)
    normUx, normUy = Ux/magU, Uy/magU
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200)
    ax.set_title(r'time = %.2fs: Pressure field'%(time),fontsize=16)
    levels = MaxNLocator(nbins=21).tick_values(magU.min(), magU.max())
    cf = ax.contourf(mx,my,magU,cmap='viridis',levels=levels,extend='both')
    #Add Discs
    AddDiscsToPlot(ax,pos)
    fig.colorbar(cf,ax=ax,shrink=0.75)
    space = 4 #spacing for quiver plot
    ax.quiver(mx[::space,::space],my[::space,::space],
              normUx[::space,::space],normUy[::space,::space],
              color='white',pivot='mid',scale=50,zorder=5)
    ax.axis([-1.0,1.0,-2.0,1.0])
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'../Figures/TestVTK/testVel.png')
    fig.clf()
    plt.close()
    return
    
### MAIN SCRIPT
#### READ ALL VTK FILES IN A SIMULATION DIRECTORY
#### CALCULATE AVERAGE FIELD DATA
#### EXPORT AVERAGED FIELD DATA TO CSV

if __name__ == '__main__':
    
    '''
    ####READ ONE VTK FILE
    
    start = t.clock()
    #Read All VTK Data into List
    cwd_DATA = cwd_PYTHON
    dumpIdx = 850
    allData = ReadAllVTK(cwd_DATA,dumpIdx)
    #Extract important values from allData
    time, mx, my, mz, wxArr, wyArr, wzArr, pArr, uxArr, uyArr, uzArr = ReadVTK(allData)
    stend = t.clock()
    diff = stend - start
    print('Time to Read VTK = %.5fs'%diff)
    #Extract Position Data
    posData = GetPosData(cwd_PYTHON,time)
    #Plot Necessary values
    start = t.clock()
    #Vorticity
    test_mx, test_my = mx[:,:,64], my[:,:,64]
    test_wzArr = wzArr[:,:,64].copy()
    PlotVortField(time,test_mx,test_my,test_wzArr,posData)
    #Streamlines
    #PlotStreamsU(time,mx,my,uxArr,uyArr,posData)
    #Pressure
    test_pArr = pArr[:,:,64].copy()
    PlotPresField(time,test_mx,test_my,test_pArr,posData)
    #Velocity Field
    test_uxArr = uxArr[:,:,64].copy()
    test_uyArr = uyArr[:,:,64].copy()
    PlotVelField(time,test_mx,test_my,test_uxArr,test_uyArr,posData)
    stend = t.clock()
    diff = stend - start
    print('Time to Plot = %.5fs'%diff)
    '''
    
    ####READ ALL VTK FILES
            
    cwd_SIM = cwd_PYTHON+'../'+Par+'/'+ParVal+'/'
    #Extract Position Data
    cwd_POS = cwd_SIM
    #Calculate # Periods
    nTime = GetPosDataLength(cwd_POS)
    nPer = int(np.trunc(1.0*nTime/1000.0))
    nSims = 20
    nPer -= 1
    #Paths to data and plots
    cwd_DATA = cwd_SIM+'VTK/'
    countPer = nPer
    countPlot = 0
    start = t.clock()
    #Create sum arrays
    Nx, Ny, Nz = 128, 512, 128
    mxsum = np.zeros((Nx,Ny,Nz))
    mysum = np.zeros((Nx,Ny,Nz))
    mzsum = np.zeros((Nx,Ny,Nz))
    wxsum = np.zeros((Nx,Ny,Nz))
    wysum = np.zeros((Nx,Ny,Nz))
    wzsum = np.zeros((Nx,Ny,Nz))
    psum = np.zeros((Nx,Ny,Nz))
    uxsum = np.zeros((Nx,Ny,Nz))
    uysum = np.zeros((Nx,Ny,Nz))
    uzsum = np.zeros((Nx,Ny,Nz))
    #Create posData sum database
    for idx in range(0,nSims):
        dumpIdx = nSims*countPer+idx
        #Read All VTK Data into List
        allData = ReadAllVTK(cwd_DATA,dumpIdx)
        #Extract important values from allData
        time,mx,my,mz,wxArr,wyArr,wzArr,pArr,uxArr,uyArr,uzArr, = ReadVTK(allData)
        #Add fields to sum
        mxsum += mx
        mysum += my
        mzsum += mz
        wxsum += wxArr
        wysum += wyArr
        wzsum += wzArr
        psum += pArr
        uxsum += uxArr
        uysum += uyArr
        uzsum += uzArr
        countPlot += 1
    #Calculate Averaged Fields
    mxavg = mxsum/(1.0*nSims)
    myavg = mysum/(1.0*nSims)
    mzavg = mzsum/(1.0*nSims)
    wxavg = wxsum/(1.0*nSims)
    wyavg = wysum/(1.0*nSims)
    wzavg = wzsum/(1.0*nSims)
    pavg = psum/(1.0*nSims)
    uxavg = uxsum/(1.0*nSims)
    uyavg = uysum/(1.0*nSims)
    uzavg = uzsum/(1.0*nSims)
    #Export Averaged Field Data
    #Create a dataframe containing flattened arrays for (mx, my, avgW, avgP, avgU)
    avgDict = {'mx':mxavg.flatten(),'my':myavg.flatten(),'mz':mzavg.flatten(),
               'avgWx':wxavg.flatten(),'avgWy':wyavg.flatten(),'avgWz':wzavg.flatten(),
               'avgP':pavg.flatten(),'avgUx':uxavg.flatten(),'avgUy':uyavg.flatten(),
               'avgUz':uzavg.flatten()}
    avgFieldData = pd.DataFrame(data=avgDict)
    #Export Avg Field Data to .csv file
    pathlib.Path(cwd_DATA+'/AVG/').mkdir(parents=True, exist_ok=True)
    avgFieldData.to_csv(cwd_DATA+'/AVG/AVG_%04d.csv'%countPer,index=False,sep=' ',float_format='%.5e')
        
    stend = t.clock()
    diff = stend - start
    print('Time to run for 1 period = %.5fs'%diff)
    sys.stdout.flush()
    #countPer += 1

