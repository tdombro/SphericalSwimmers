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
from scipy import interpolate

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
    return cwd+"pd.txt"

def GetPosData(cwd,time):
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

### AVERAGE FIELD DATA FUNCTIONS
def GetAvgFieldData(cwd,idx):
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    #cwd = cwd_PYTHON
    fieldData = pd.read_csv(cwd+'AVG_XY_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    UxList = fieldData['avgUx'].values.tolist()
    UyList = fieldData['avgUy'].values.tolist()
    WList = fieldData['avgWz'].values.tolist()
    PList = fieldData['avgP'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 128, 512
    mxArr = np.array(mxList).reshape((Nx,Ny))
    myArr = np.array(myList).reshape((Nx,Ny))
    UxArr = np.array(UxList).reshape((Nx,Ny))
    UyArr = np.array(UyList).reshape((Nx,Ny))
    WArr = np.array(WList).reshape((Nx,Ny))
    PArr = np.array(PList).reshape((Nx,Ny))
    return (mxArr, myArr, UxArr, UyArr, WArr, PArr)

def GetAvgBiPFieldData(cwd,idx):
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    #cwd = cwd_PYTHON
    fieldData = pd.read_csv(cwd+'AVGB_XY_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    UxList = fieldData['avgUx'].values.tolist()
    UyList = fieldData['avgUy'].values.tolist()
    PList = fieldData['avgP'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 256, 256
    mxArr = np.array(mxList).reshape((Nx,Ny))
    myArr = np.array(myList).reshape((Nx,Ny))
    UxArr = np.array(UxList).reshape((Nx,Ny))
    UyArr = np.array(UyList).reshape((Nx,Ny))
    PArr = np.array(PList).reshape((Nx,Ny))
    return (mxArr, myArr, UxArr, UyArr, PArr)

###INSTANTANEOUS FIELD FUNCTIONS
def GetInstBiPFieldData(cwd,idx):
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    #cwd = cwd_PYTHON
    fieldData = pd.read_csv(cwd+'Transformed_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mxBip'].values.tolist()
    myList = fieldData['myBip'].values.tolist()
    UxList = fieldData['UxBip'].values.tolist()
    UyList = fieldData['UyBip'].values.tolist()
    PList = fieldData['pBip'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 256, 256
    mxArr = np.array(mxList).reshape((Nx,Ny))
    myArr = np.array(myList).reshape((Nx,Ny))
    UxArr = np.array(UxList).reshape((Nx,Ny))
    UyArr = np.array(UyList).reshape((Nx,Ny))
    PArr = np.array(PList).reshape((Nx,Ny))
    return (mxArr, myArr, UxArr, UyArr, PArr)

def GetInstFieldData(cwd,idx):
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    #cwd = cwd_PYTHON
    fieldData = pd.read_csv(cwd+'RAW_XY_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    UxList = fieldData['Ux'].values.tolist()
    UyList = fieldData['Uy'].values.tolist()
    PList = fieldData['P'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 128, 512
    mxArr = np.array(mxList).reshape((Nx,Ny))
    myArr = np.array(myList).reshape((Nx,Ny))
    UxArr = np.array(UxList).reshape((Nx,Ny))
    UyArr = np.array(UyList).reshape((Nx,Ny))
    PArr = np.array(PList).reshape((Nx,Ny))
    return (mxArr, myArr, UxArr, UyArr, PArr)

###INTERPOLATE BIPOLAR DATA TO CARTESIAN
###ONLY AVG: INSTANTANEOUS GIVES NO NEW INFO

def InterpolateToNewCoordinateSystem(mx,my,UxCart,UyCart,pCart,Nx,Ny):
    #Interpolate Ux and Uy from original cartesian coordainates to new ones
    x = np.linspace(-1.0,1.0,Nx)
    y = np.linspace(-1.5,1.5,Ny)
    mx_new, my_new = np.meshgrid(x,y)
    #Griddata
    Ux_new=interpolate.griddata((mx.flatten(),my.flatten()),UxCart.flatten() , (mx_new,my_new),method='linear')
    Uy_new=interpolate.griddata((mx.flatten(),my.flatten()),UyCart.flatten() , (mx_new,my_new),method='linear')
    p_new=interpolate.griddata((mx.flatten(),my.flatten()),pCart.flatten() , (mx_new,my_new),method='linear')
    
    print('Coordinate Transformation Complete!')
    return (mx_new,my_new,Ux_new,Uy_new,p_new)

###PLOT AVG FIELD DATA FUNCTIONS
# constructs a filepath for the plot images of Re=$Re, config=$config, and field=$field
def plotavgName(cwd,eps,s,field):
    #field = Vort, Pres, Vel, AvgW, AvgP, AvgU
    strDir = cwd+"../Figures/Transformed/AVG/{0}/{1}/".format(eps,s)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return strDir+"AVG_{0}_.png".format(field)

def plotinstName(cwd,eps,s,field,idx):
    #field = Vort, Pres, Vel, AvgW, AvgP, AvgU
    strDir = cwd+"../Figures/Transformed/Inst/{0}/{1}/".format(eps,s)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return strDir+"Inst_{0}_{1}_.png".format(field,idx)

def AddDiscsToPlot(ax,pos,CartoBip):
    RADIUS_LARGE = 0.3
    RADIUS_SMALL = 0.15
    #Add Discs
    if(CartoBip == 'Cart'):
        circle1 = Circle((pos['xU'], pos['yU']), RADIUS_LARGE, facecolor='k',
                         linewidth=1,alpha=1.0,zorder=6)
        ax.add_patch(circle1)
        circle2 = Circle((pos['xL'], pos['yL']), RADIUS_SMALL, facecolor='k',
                         linewidth=1,alpha=1.0,zorder=6)
        ax.add_patch(circle2)
    else:
        circle1 = Circle((0.0, 0.45), RADIUS_LARGE, facecolor='gray',
                         linewidth=1,alpha=1.0,zorder=6)
        ax.add_patch(circle1)
        circle2 = Circle((0.0, -0.45), RADIUS_SMALL, facecolor='gray',
                         linewidth=1,alpha=1.0,zorder=6)
        ax.add_patch(circle2)
    return

def PlotVelocityField(cwd,time,mx,my,Ux,Uy,pos,pars):
    eps = pars[0]
    sval = pars[1]
    AvgoInst = pars[2]
    field = pars[3]
    CartoBip = pars[4]
    idx = pars[5]
    magU = np.hypot(Ux,Uy)
    print('Ux = ',Ux[0:10,0])
    print('Uy = ',Uy[0:10,0])
    print('magU = ',magU[0:10,0])
    print('maxMagU = ',np.nanmax(magU))
    print('minMagU = ',np.nanmin(magU))
    normUx, normUy = Ux/magU, Uy/magU
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,9),dpi=200)
    ax.set_title(r'time = %.2fs: Velocity field'%(time),fontsize=16)
    levels = MaxNLocator(nbins=21).tick_values(np.nanmin(magU), 0.5*np.nanmax(magU))
    cf = ax.contourf(mx,my,magU,cmap='viridis',levels=levels,extend='both')
    #Add Discs
    AddDiscsToPlot(ax,pos,CartoBip)
    fig.colorbar(cf,ax=ax,shrink=0.75)
    space = 6 #spacing for quiver plot
    ax.quiver(mx[::space,::space],my[::space,::space],
              normUx[::space,::space],normUy[::space,::space],
              color='white',pivot='mid',scale=25,zorder=5,width=0.005*1.0,headwidth=2.5,headlength=3.0)
    '''
    if(field == 'U'):
        space = 6 #spacing for quiver plot
        ax.quiver(mx[::space,::space],my[::space,::space],
                  normUx[::space,::space],normUy[::space,::space],
                  color='white',pivot='mid',scale=25,zorder=5,width=0.005*2.0,headwidth=2.0,headlength=3.0)
    else:
        space = 6 #spacing for quiver plot
        ax.quiver(mx[::space,::space],my[::space,::space],
                  normUx[::space,::space],normUy[::space,::space],
                  color='red',pivot='mid',scale=40,zorder=5,width=0.005*1.0,headwidth=2.0,headlength=3.0)
    '''
    ax.axis([-1.0,1.0,-1.5,1.5])
    #ax.axis([-1.0,1.0,-4.0,4.0])
    ax.set_aspect('equal')
    fig.tight_layout()
    if(AvgoInst == 'AVG'):
        fig.savefig(plotavgName(cwd,eps,sval,field))
    else:
        fig.savefig(plotinstName(cwd,eps,sval,field,idx))
    fig.clf()
    plt.close()
    return

#Plot Streamlines for interpolated velocity fields Ux and Uy
def PlotVelocityStream(cwd,time,mx,my,Ux,Uy,pos,pars):
    eps = pars[0]
    sval = pars[1]
    AvgoInst = pars[2]
    field = pars[3]
    CartoBip = pars[4]
    idx = pars[5]
    #Here, we will visualize the velocity field on the new coordinate system
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,9),dpi=200)
    ax.set_title(r'time = %.2fs: Streamlines'%(time),fontsize=16)
    #Contour
    magU = np.hypot(Ux,Uy)
    levels = MaxNLocator(nbins=21).tick_values(magU.min(), 0.5*magU.max())
    ax.contourf(mx,my,magU,cmap='viridis',levels=levels,extend='both')
    #Add Discs
    AddDiscsToPlot(ax,pos,CartoBip)
    #Streamlines
    x = np.linspace(-1.0,1.0,128)
    y = np.linspace(-4.0,4.0,512)
    ax.streamplot(x,y,Ux.T,Uy.T,color='white',linewidth=1.5,arrowsize=1.0,density=3)
    ax.axis([-1.0,1.0,-1.5,1.5])
    ax.set_aspect('equal')
    fig.tight_layout()
    if(AvgoInst == 'AVG'):
        fig.savefig(plotavgName(cwd,eps,sval,field))
    else:
        fig.savefig(plotinstName(cwd,eps,sval,field,idx))
    fig.clf()
    plt.close()
    return

def PlotPressureField(cwd,time,mx,my,p,pos,pars):
    eps = pars[0]
    sval = pars[1]
    AvgoInst = pars[2]
    field = pars[3]
    CartoBip = pars[4]
    idx = pars[5]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,9),dpi=200)
    ax.set_title(r'time = %.2fs: Pressure field'%(time),fontsize=16)
    #Contour
    levels = MaxNLocator(nbins=21).tick_values(p.min(), p.max())
    #levels = MaxNLocator(nbins=21).tick_values(-1, 1)
    cf = ax.contourf(mx,my,p,cmap='PRGn',levels=levels,extend='both')
    fig.colorbar(cf,ax=ax,shrink=0.75)
    #Add Discs
    AddDiscsToPlot(ax,pos,CartoBip)
    ax.axis([-1.0,1.0,-1.5,1.5])
    ax.set_aspect('equal')
    fig.tight_layout()
    if(AvgoInst == 'AVG'):
        fig.savefig(plotavgName(cwd,eps,sval,field))
    else:
        fig.savefig(plotinstName(cwd,eps,sval,field,idx))
    fig.clf()
    plt.close()
    return

def PlotVorticityStream(cwd,time,mx,my,w,Ux,Uy,pos,pars):
    eps = pars[0]
    sval = pars[1]
    AvgoInst = pars[2]
    field = pars[3]
    CartoBip = pars[4]
    idx = pars[5]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,9),dpi=200,num=10)
    ax.set_title(r'time = %.2fs'%(time),fontsize=16)
    #Vorticity Field (Contour)
    levels = MaxNLocator(nbins=21).tick_values(-5, 5)
    ax.contourf(mx,my,w,cmap='bwr',levels=levels,extend='both')
    #Streamlines
    x = np.linspace(-1.0,1.0,128)
    y = np.linspace(-4.0,4.0,512)
    ax.streamplot(x,y,Ux.T,Uy.T,color='k',linewidth=1.5,arrowsize=1.0,density=3)
    #Add Discs
    AddDiscsToPlot(ax,pos,CartoBip)
    ax.axis([-1.0,1.0,-1.5,1.5])
    ax.set_aspect('equal')
    # Turn off tick labels
    #ax.set_yticklabels([])
    #ax.set_xticklabels([])
    fig.tight_layout()
    #ax = set_size(6,6,ax)
    fig.savefig(plotavgName(cwd,eps,sval,field))
    fig.clf()
    plt.close()
    return

### MAIN SCRIPT
#### READ ALL CSV FILES
#### PLOT AVERAGE FIELD DATA (RAW AND TRANSFORMED)
#### PLOT INSTANTANEOUS FIELD DATA (RAW AND TRANSFORMED)

if __name__ == '__main__':
    
    ####READ ALL CSV FILES
    start = t.clock()
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
    ###Obtain Instantaneous Field Data for idx
    for idx in range(nSims):
        mx,my,Ux,Uy,P = GetInstFieldData(cwd_DATA,idx)
        mxB,myB,UxB,UyB,PB = GetInstBiPFieldData(cwd_DATA,idx)
        time = np.round(countPer*PERIOD + 0.005*idx,3)
        posData = GetPosData(cwd_SIM,time)#Extract Position and Time Data
        #Plot Velocity Field and Pressure
        pars = [Par,ParVal,'Inst','U','Cart',idx]
        PlotVelocityField(cwd_PYTHON,time,mx,my,Ux,Uy,posData,pars)
        pars[3] = 'UStream'
        PlotVelocityStream(cwd_PYTHON,time,mx,my,Ux,Uy,posData,pars)
        pars[3] = 'P'
        PlotPressureField(cwd_PYTHON,time,mx,my,P,posData,pars)
        pars[4] = 'Bip'
        pars[3] = 'UBiP'
        PlotVelocityField(cwd_PYTHON,time,mxB,myB,UxB,UyB,posData,pars)
        pars[3] = 'PBiP'
        PlotPressureField(cwd_PYTHON,time,mxB,myB,PB,posData,pars)
        
    ###Obtain AVG Field Data for Per countPer
    cwd_AVG = cwd_DATA+'AVG/'
    avgmx,avgmy,avgUx,avgUy,avgW, avgP = GetAvgFieldData(cwd_AVG,countPer)
    avgmxB,avgmyB,avgUxB,avgUyB, avgPB = GetAvgBiPFieldData(cwd_AVG,countPer)
    #Interpolate Bipolar values onto a cartesian grid for plotting and visualization purposes
    avgmxB_cart,avgmyB_cart,avgUxB_cart,avgUyB_cart,avgPB_cart = InterpolateToNewCoordinateSystem(avgmxB,avgmyB,avgUxB,avgUyB,avgPB,128,192)
    #Extract Position and Time Data
    time = np.round(0.05 + countPer*PERIOD,2)
    #print('time = ',time)
    posData = GetPosData(cwd_SIM,time)#Extract Position and Time Data
    #Plot Velocity Field and Pressure
    pars = [Par,ParVal,'AVG','U','Cart',countPer]
    PlotVelocityField(cwd_PYTHON,time,avgmx,avgmy,avgUx,avgUy,posData,pars)
    pars[3] = 'UStream'
    PlotVelocityStream(cwd_PYTHON,time,avgmx,avgmy,avgUx,avgUy,posData,pars)
    pars[3] = 'P'
    PlotPressureField(cwd_PYTHON,time,avgmx,avgmy,avgP,posData,pars)
    pars[3] = 'W'
    PlotVorticityStream(cwd_PYTHON,time,avgmx,avgmy,avgW,avgUx,avgUy,posData,pars)
    pars[4] = 'Bip'
    pars[3] = 'UBiP'
    #PlotVelocityField(cwd_PYTHON,time,avgmxB,avgmyB,avgUxB,avgUyB,posData,pars)
    PlotVelocityField(cwd_PYTHON,time,avgmxB_cart,avgmyB_cart,avgUxB_cart,avgUyB_cart,posData,pars)
    pars[3] = 'UBiPStream'
    #PlotVelocityStream(cwd_PYTHON,time,avgmxB_cart,avgmyB_cart,avgUxB_cart,avgUyB_cart,posData,pars)
    pars[3] = 'PBiP'
    PlotPressureField(cwd_PYTHON,time,avgmxB,avgmyB,avgPB,posData,pars)
    #PlotPressureField(cwd_PYTHON,time,avgmxB_cart,avgmyB_cart,avgPB_cart,posData,pars)

    stend = t.clock()
    diff = stend - start
    print('Time to run for 1 period = %.5fs'%diff)
    sys.stdout.flush()
    #countPer += 1


