#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 15:33:01 2020

@author: thomas
"""

#Import modules
import numpy as np
import pandas as pd
import os, sys
import time as t
import matplotlib as mpl
mpl.use('Agg')
import pathlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
Par = sys.argv[1]
ParVal = sys.argv[2]

#SIM PARAMETERS
PERIOD = 0.1

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

### AVERAGE FIELD DATA FUNCTIONS
def GetAvgFieldData(cwd,idx):
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    #cwd = cwd_PYTHON
    fieldData = pd.read_csv(cwd+'AVG_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    mzList = fieldData['mz'].values.tolist()
    UxList = fieldData['avgUx'].values.tolist()
    UyList = fieldData['avgUy'].values.tolist()
    UzList = fieldData['avgUz'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny, Nz = 128, 512, 128
    mxArr = np.array(mxList).reshape((Nx,Ny,Nz))
    myArr = np.array(myList).reshape((Nx,Ny,Nz))
    mzArr = np.array(mzList).reshape((Nx,Ny,Nz))
    UxArr = np.array(UxList).reshape((Nx,Ny,Nz))
    UyArr = np.array(UyList).reshape((Nx,Ny,Nz))
    UzArr = np.array(UzList).reshape((Nx,Ny,Nz))
    return (mxArr, myArr, mzArr, UxArr, UyArr, UzArr)

def SubsetZYPlane(my,mz,uy,uz):
    Nx = 128
    mym = my[int(0.5*Nx-1),:,:].copy()
    myp = my[int(0.5*Nx),:,:].copy()
    my0 = 0.5*(mym+myp)
    mzm = mz[int(0.5*Nx-1),:,:].copy()
    mzp = mz[int(0.5*Nx),:,:].copy()
    mz0 = 0.5*(mzm+mzp)
    uym = uy[int(0.5*Nx-1),:,:].copy()
    uyp = uy[int(0.5*Nx),:,:].copy()
    uy0 = 0.5*(uym+uyp)
    uzm = uz[int(0.5*Nx-1),:,:].copy()
    uzp = uz[int(0.5*Nx),:,:].copy()
    uz0 = 0.5*(uzm+uzp)
    
    return (my0.T,mz0.T,uy0.T,uz0.T)

def SubsetXYPlane(mx,my,ux,uy):
    Nz = 128
    mxm = mx[:,:,int(0.5*Nz-1)].copy()
    mxp = mx[:,:,int(0.5*Nz)].copy()
    mx0 = 0.5*(mxm+mxp)
    mym = my[:,:,int(0.5*Nz-1)].copy()
    myp = my[:,:,int(0.5*Nz)].copy()
    my0 = 0.5*(mym+myp)
    uxm = ux[:,:,int(0.5*Nz-1)].copy()
    uxp = ux[:,:,int(0.5*Nz)].copy()
    ux0 = 0.5*(uxm+uxp)
    uym = uy[:,:,int(0.5*Nz-1)].copy()
    uyp = uy[:,:,int(0.5*Nz)].copy()
    uy0 = 0.5*(uym+uyp)
    
    return (mx0,my0,ux0,uy0)

###PLOT AVG FIELD DATA FUNCTIONS
# constructs a filepath for the plot omages of Re=$Re, config=$config, and field=$field
def plotavgName(cwd,par,field,plane,idx):
    #field = Vort, Pres, Vel, AvgW, AvgP, AvgU
    strDir = cwd+"../Figures/AVG/Last/{0}/{1}/{2}/".format(par,field,plane)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return strDir+"{0}_{1}.png".format(field,idx)

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

def PlotAvgVelField(cwd,time,mx,my,Ux,Uy,pos,pars):
    par = pars[0]
    field = pars[1]
    plane = pars[2]
    idx = pars[3]
    magU = np.hypot(Ux,Uy)
    normUx, normUy = Ux/magU, Uy/magU
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200)
    ax.set_title(r'time = %.2fs: Velocity field'%(time),fontsize=16)
    levels = MaxNLocator(nbins=21).tick_values(magU.min(), magU.max())
    cf = ax.contourf(mx,my,magU,cmap='viridis',levels=levels,extend='both')
    #Add Discs
    AddDiscsToPlot(ax,pos)
    fig.colorbar(cf,ax=ax,shrink=0.75)
    space = 4 #spacing for quiver plot
    ax.quiver(mx[::space,::space],my[::space,::space],
              normUx[::space,::space],normUy[::space,::space],
              color='white',pivot='mid',scale=50,zorder=5)
    ax.axis([-1.0,1.0,-1.5,1.5])
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(plotavgName(cwd,par,field,plane,idx))
    fig.clf()
    plt.close()
    return

### HAT_U AUXILIARY FUNCTIONS
def Trapezoidal(data, minBound, maxBound, nTime):
    h = (maxBound-minBound)/float(nTime)
    print(h)
    s = 0.5*(data[0] + data[nTime-1])
    for idx in range(1,nTime-1):
        s = s + data[idx]
    return h*s

def PlothatUMagField(cwd,time,mx,my,U,pos,pars):
    par = pars[0]
    field = pars[1]
    plane = pars[2]
    idx = pars[3]
    Uxc, Uxs, Uyc, Uys = U[0].copy(), U[1].copy(), U[2].copy(), U[3].copy()
    magUx = np.hypot(Uxc,Uxs)
    magUy = np.hypot(Uyc,Uys)
    fig, ax = plt.subplots(nrows=1, ncols=2,figsize=(12,6),dpi=200)
    ax[0].set_title(r'time = %.2fs: $\hat{U}_x$ field'%(time),fontsize=16)
    ax[1].set_title(r'time = %.2fs: $\hat{U}_y$ field'%(time),fontsize=16)
    #hatUx
    levels = MaxNLocator(nbins=21).tick_values(magUx.min(), magUx.max())
    cf = ax[0].contourf(mx,my,magUx,cmap='viridis',levels=levels,extend='both')
    #Add Discs
    AddDiscsToPlot(ax[0],pos)
    fig.colorbar(cf,ax=ax[0],shrink=0.75)
    #hatUy
    levels = MaxNLocator(nbins=21).tick_values(magUy.min(), magUy.max())
    cf = ax[1].contourf(mx,my,magUy,cmap='viridis',levels=levels,extend='both')
    #Add Discs
    AddDiscsToPlot(ax[1],pos)
    fig.colorbar(cf,ax=ax[1],shrink=0.75)
    ax[0].axis([-1.0,1.0,-1.5,1.5])
    ax[0].set_aspect('equal')
    ax[1].axis([-1.0,1.0,-1.5,1.5])
    ax[1].set_aspect('equal')
    fig.tight_layout()
    fig.savefig(plotavgName(cwd,par,field,plane,idx))
    fig.clf()
    plt.close()
    return

def PlothatUCSField(cwd,time,mx,my,U,pos,pars):
    par = pars[0]
    field = pars[1]
    plane = pars[2]
    idx = pars[3]
    Uxc, Uxs, Uyc, Uys = U[0].copy(), U[1].copy(), U[2].copy(), U[3].copy()
    magUc = np.hypot(Uxc,Uyc)
    magUs = np.hypot(Uxs,Uys)
    fig, ax = plt.subplots(nrows=1, ncols=2,figsize=(12,6),dpi=200)
    ax[0].set_title(r'time = %.2fs: $\hat{U}cos$ field'%(time),fontsize=16)
    ax[1].set_title(r'time = %.2fs: $\hat{U}sin$ field'%(time),fontsize=16)
    #hatUx
    levels = MaxNLocator(nbins=21).tick_values(magUc.min(), magUc.max())
    cf = ax[0].contourf(mx,my,magUc,cmap='viridis',levels=levels,extend='both')
    #Vector field for hat{U}cos
    space = 4 #spacing for quiver plot
    normUxc, normUyc = Uxc/magUc, Uyc/magUc
    ax[0].quiver(mx[::space,::space],my[::space,::space],
              normUxc[::space,::space],normUyc[::space,::space],
              color='white',pivot='mid',scale=50,zorder=5)
    #Add Discs
    AddDiscsToPlot(ax[0],pos)
    fig.colorbar(cf,ax=ax[0],shrink=0.75)
    #hatUy
    levels = MaxNLocator(nbins=21).tick_values(magUs.min(), magUs.max())
    cf = ax[1].contourf(mx,my,magUs,cmap='viridis',levels=levels,extend='both')
    #Vector field for hat{U}cos
    space = 4 #spacing for quiver plot
    normUxs, normUys = Uxs/magUs, Uys/magUs
    ax[1].quiver(mx[::space,::space],my[::space,::space],
              normUxs[::space,::space],normUys[::space,::space],
              color='white',pivot='mid',scale=50,zorder=5)
    #Add Discs
    AddDiscsToPlot(ax[1],pos)
    fig.colorbar(cf,ax=ax[1],shrink=0.75)
    ax[0].axis([-1.0,1.0,-1.5,1.5])
    ax[0].set_aspect('equal')
    ax[1].axis([-1.0,1.0,-1.5,1.5])
    ax[1].set_aspect('equal')
    fig.tight_layout()
    fig.savefig(plotavgName(cwd,par,field,plane,idx))
    fig.clf()
    plt.close()
    return


### MAIN SCRIPT
#### READ ALL VTK FILES IN A SIMULATION DIRECTORY
#### CALCULATE AVERAGE FIELD DATA
#### EXPORT AVERAGED FIELD DATA TO CSV

if __name__ == '__main__':
    
    ####READ ALL VTK FILES
            
    cwd_SIM = cwd_PYTHON+'../'+Par+'/'+ParVal+'/'
    #Extract Position Data
    cwd_POS = cwd_SIM
    #Calculate # Periods
    nTime = GetPosDataLength(cwd_POS)
    nPer = int(np.trunc(1.0*nTime/1000.0))
    print('nPer = ',nPer)
    sys.stdout.flush()
    nSims = 20
    nPer -= 1
    #Paths to data and plots
    cwd_DATA = cwd_SIM+'VTK/'
    countPer = nPer
    countPlot = 0
    ###Obtain AVG Field Data for Per countPer
    cwd_AVG = cwd_DATA+'AVG/'
    avgmx,avgmy,avgmz,avgUx,avgUy,avgUz = GetAvgFieldData(cwd_AVG,countPer)
    #Extract Position and Time Data
    time = np.round(0.05 + countPer*PERIOD,2)
    #print('time = ',time)
    posData = GetPosData(cwd_SIM,time)#Extract Position and Time Data
    '''###Plot AVG Field Data
    #XY Plane
    #Velocity Field
    #pars = [ParVal,'avgU','XY',countPer]
    #PlotAvgVelField(cwd_PYTHON,time,avgmx_z0,avgmy_z0,avgUx_z0,avgUy_z0,posData,pars)
    #ZY Plane
    #Velocity Field
    #pars[2] = 'ZY'
    #PlotAvgVelField(cwd_PYTHON,time,avgmz_x0,avgmy_x0,avgUz_x0,avgUy_x0,posData,pars)
    #Difference
    #Velocity Field
    #diffUpar = abs(avgUy_x0 - avgUy_z0)
    #diffUperp = abs(avgUz_x0 - avgUx_z0)
    #pars[2] = 'Diff'
    #PlotAvgVelField(cwd_PYTHON,time,avgmx_z0,avgmy_z0,diffUperp,diffUpar,posData,pars)
    '''
    
    start = t.clock()
    #Create sum arrays
    Nx, Ny, Nz = 128, 512, 128
    #Create Arrays for XY and ZY Planes
    #hatU terms
    #0 = xcos; 1 = x(isin); 2 = ycos; 3 = y(isin); 4 = zcos; 5 = z(isin)
    hatU_terms = np.zeros((nSims,6,Nx,Ny,Nz))

    #Create posData sum database
    for idx in range(0,nSims):
        vtkStart = t.clock()
        dumpIdx = nSims*countPer+idx
        #Read All VTK Data into List
        allData = ReadAllVTK(cwd_DATA,dumpIdx)
        #Extract important values from allData
        time,mx,my,mz,wxArr,wyArr,wzArr,pArr,uxArr,uyArr,uzArr, = ReadVTK(allData)
        #Calculate hatU terms
        hatU_terms[idx,0] = (uxArr - avgUx)*np.cos(2.0*time*np.pi/PERIOD)
        hatU_terms[idx,2] = (uyArr - avgUy)*np.cos(2.0*time*np.pi/PERIOD)
        hatU_terms[idx,4] = (uzArr - avgUz)*np.cos(2.0*time*np.pi/PERIOD)
        hatU_terms[idx,1] = (uxArr - avgUx)*np.sin(2.0*time*np.pi/PERIOD)
        hatU_terms[idx,3] = (uyArr - avgUy)*np.sin(2.0*time*np.pi/PERIOD)
        hatU_terms[idx,5] = (uzArr - avgUz)*np.sin(2.0*time*np.pi/PERIOD)
        
        #Add fields to sum
        vtkEnd = t.clock()
        diff = vtkEnd - vtkStart
        print('Time to read and plot 1 vtk = %.5fs'%diff)
        sys.stdout.flush()
        countPlot += 1
    #Calculate hatU Fields
    hatU = (1.0/PERIOD)*Trapezoidal(hatU_terms,0,PERIOD,nSims)
    '''###Plot hatU Field Data
    #time = np.round(0.05 + countPer*PERIOD,2)
    #XY
    #pars = [ParVal,'hatU','XY',countPer]
    #PlothatUMagField(cwd_PYTHON,time,avgmx_z0,avgmy_z0,hatU_z0,posData,pars)
    #pars[1] = 'cossin'
    #PlothatUCSField(cwd_PYTHON,time,avgmx_z0,avgmy_z0,hatU_z0,posData,pars)
    #ZY
    #pars = [ParVal,'hatU','ZY',countPer]
    #PlothatUMagField(cwd_PYTHON,time,avgmz_x0,avgmy_z0,hatU_x0,posData,pars)
    #pars[1] = 'cossin'
    #PlothatUCSField(cwd_PYTHON,time,avgmz_x0,avgmy_x0,hatU_x0,posData,pars)
    '''

    #Export hatU Field Data
    #Create a dataframe containing flattened arrays for (mx, my, avgW, avgP, avgU)
    hatUDict = {'mx':avgmx.flatten(),'my':avgmy.flatten(),'mz':avgmz.flatten(),
                'hatUxc':hatU[0].flatten(),'hatUxs':hatU[1].flatten(),
                'hatUyc':hatU[2].flatten(),'hatUys':hatU[3].flatten(),
                'hatUzc':hatU[4].flatten(),'hatUzs':hatU[5].flatten()}
    hatUFieldData = pd.DataFrame(data=hatUDict)
    #Export Avg Field Data to .csv file
    pathlib.Path(cwd_DATA+'/AVG/').mkdir(parents=True, exist_ok=True)
    hatUFieldData.to_csv(cwd_DATA+'/AVG/allhatU_%04d.csv'%countPer,index=False,sep=' ',float_format='%.5e')
        
    stend = t.clock()
    diff = stend - start
    print('Time to run for 1 period = %.5fs'%diff)
    sys.stdout.flush()
    #countPer += 1


