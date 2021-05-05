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
from scipy import interpolate

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
Par = sys.argv[1]
ParVal = sys.argv[2]

#AUXILIARY FUNCTIONS
#Position Data Functions
# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    return cwd+"pd.txt"

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

### DATA PROCESSING ###
### BIPOLAR TRANSFORMATION, ONLY NEED X AND Y ###
def SubsetXYPlane_Vector(mx,my,xvar,yvar,zvar,NZ):
    mxm = mx[:,:,int(0.5*NZ-1)].copy()
    mxp = mx[:,:,int(0.5*NZ)].copy()
    mx0 = 0.5*(mxm+mxp)
    mym = my[:,:,int(0.5*NZ-1)].copy()
    myp = my[:,:,int(0.5*NZ)].copy()
    my0 = 0.5*(mym+myp)
    xvarm = xvar[:,:,int(0.5*NZ-1)].copy()
    xvarp = xvar[:,:,int(0.5*NZ)].copy()
    xvar0 = 0.5*(xvarm+xvarp)
    yvarm = yvar[:,:,int(0.5*NZ-1)].copy()
    yvarp = yvar[:,:,int(0.5*NZ)].copy()
    yvar0 = 0.5*(yvarm+yvarp)
    zvarm = zvar[:,:,int(0.5*NZ-1)].copy()
    zvarp = zvar[:,:,int(0.5*NZ)].copy()
    zvar0 = 0.5*(zvarm+zvarp)
    
    return (mx0,my0,xvar0,yvar0,zvar0)

def SubsetXYPlane_Scalar(var,NZ):
    varm = var[:,:,int(0.5*NZ-1)].copy()
    varp = var[:,:,int(0.5*NZ)].copy()
    var0 = 0.5*(varm+varp)
    
    return (var0)

def InterpolateToNewCoordinateSystem(cwd,time,mx,my,UxCart,UyCart,pCart,Nx,Ny):
    #Create a uniform bipolar mesh for the interpolated fields!
    ###BIPOLAR###
    #Given r_1, r_2, d
    #eta = np.linspace(-np.pi,np.pi,Nx)
    #xi_bar = np.linspace(0.0,1.0,Ny)
    #print('xi_bar_Int = ',xi_bar)
    #meta, mxi_bar = np.meshgrid(eta,xi_bar)
    #Calculate a, xi1, xi2
    r_2 = 0.15
    r_1 = 2.0*r_2
    pos = GetPosData(cwd,time) #xU,xL,yU,yL,curr,des,time
    dx = pos.loc[0,'xU'] - pos.loc[0,'xL']
    dy = pos.loc[0,'yU'] - pos.loc[0,'yL']
    d = np.hypot(dx,dy)
    a,meta,mxi,xi_1,xi_2 = CalculateBipolarMesh(r_1,r_2,d,Nx,Ny)
    h = a/(np.cosh(mxi) - np.cos(meta))
    f_1 = 1.0/abs(np.tanh(xi_1))
    f_2 = 1.0/abs(np.tanh(xi_2))
    shift_x = pos.loc[0,'xU'] + (f_1/(f_1+f_2))*(pos.loc[0,'xL'] - pos.loc[0,'xU'])#TEMP
    shift_y = pos.loc[0,'yU'] + (f_1/(f_1+f_2))*(pos.loc[0,'yL'] - pos.loc[0,'yU'])#TEMP
    my_new = h*np.sinh(mxi) + shift_y
    mx_new = h*np.sin(meta) + shift_x
    ###END BIPOLAR
    #print('mx_new = ',mx_new)
    #print('my_new = ',my_new)
    #Interpolate Ux and Uy from original cartesian coordainates to new ones
    #Griddata
    UxBip=interpolate.griddata((mx.flatten(),my.flatten()),UxCart.flatten() , (mx_new,my_new),method='cubic')
    UyBip=interpolate.griddata((mx.flatten(),my.flatten()),UyCart.flatten() , (mx_new,my_new),method='cubic')
    pBip=interpolate.griddata((mx.flatten(),my.flatten()),pCart.flatten() , (mx_new,my_new),method='cubic')
    
    print('Coordinate Transformation Complete!')
    return (mx_new,my_new,UxBip,UyBip,pBip)

def CalculateBipolarMesh(r_1,r_2,d,Nx,Ny):
    eta = np.linspace(-np.pi,np.pi,Nx)
    xi_bar = np.linspace(0.0,1.0,Ny)
    v_2 = 2*r_1*r_2/d**2
    v_1 = (r_1+r_2)/d
    v = 1 - v_1**2 + v_2
    a = 0.5*d*np.sqrt(v**2 - v_2**2)
    xi_1 = np.arcsinh(a/r_1)
    xi_2 = -np.arcsinh(a/r_2)
    #Calc non-normalized xi
    #print('xi_bar = ',xi_bar)
    xi = xi_bar*(xi_2 - xi_1) + xi_1 #renormalization (squishes/expands mesh) (okay with this)
    meta, mxi = np.meshgrid(eta,xi)

    return (a,meta,mxi,xi_1,xi_2)
    
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
    nSims = 20
    nPer -= 1
    #Paths to data and plots
    cwd_DATA = cwd_SIM+'VTK/'
    countPer = nPer
    countPlot = 0
    start = t.clock()
    #Create sum arrays
    Nx, Ny, Nz = 128, 512, 128
    mxsum = np.zeros((Nx,Ny))
    mysum = np.zeros((Nx,Ny))
    wxsum = np.zeros((Nx,Ny))
    wysum = np.zeros((Nx,Ny))
    wzsum = np.zeros((Nx,Ny))
    psum = np.zeros((Nx,Ny))
    uxsum = np.zeros((Nx,Ny))
    uysum = np.zeros((Nx,Ny))
    uzsum = np.zeros((Nx,Ny))
    #Bipolar
    NBP = 256
    pBsum = np.zeros((NBP,NBP))
    uxBsum = np.zeros((NBP,NBP))
    uyBsum = np.zeros((NBP,NBP))
    #Create posData sum database
    for idx in range(0,nSims):
        dumpIdx = nSims*countPer+idx
        #Read All VTK Data into List
        allData = ReadAllVTK(cwd_DATA,dumpIdx)
        #Extract important values from allData
        time,mx,my,mz,wxArr,wyArr,wzArr,pArr,uxArr,uyArr,uzArr, = ReadVTK(allData)
        
        #Given Cartesian Coordinates and values for U
        #Subset to XY plane
        mx_z0,my_z0,Ux_z0,Uy_z0,Uz_z0 = SubsetXYPlane_Vector(mx,my,uxArr,uyArr,uzArr,Nz)
        mx_z0,my_z0,wx_z0,wy_z0,wz_z0 = SubsetXYPlane_Vector(mx,my,wxArr,wyArr,wzArr,Nz)
        p_z0 = SubsetXYPlane_Scalar(pArr,Nz)
        #Convert to Bipolar coordinates (128x512)?
        mxBip,myBip,UxBip,UyBip,pBip = InterpolateToNewCoordinateSystem(cwd_SIM,time,mx_z0,my_z0,Ux_z0,Uy_z0,p_z0,NBP,NBP)
        #Save cartesian/bipolar XY mid-plane data in a dataframe
        coordDict = {'mxBip':mxBip.flatten(),'myBip':myBip.flatten(),
                     'UxBip':UxBip.flatten(),'UyBip':UyBip.flatten(),'pBip':pBip.flatten()}
        coordData = pd.DataFrame(data=coordDict)
        coordData.to_csv(cwd_DATA+'/Transformed_%04d.csv'%idx,index=False,sep=' ',float_format='%.6e')
        rawDict = {'mx':mx_z0.flatten(),'my':my_z0.flatten(),
                   'Ux':Ux_z0.flatten(),'Uy':Uy_z0.flatten(),'P':p_z0.flatten()}
        rawData = pd.DataFrame(data=rawDict)
        rawData.to_csv(cwd_DATA+'/RAW_XY_%04d.csv'%idx,index=False,sep=' ',float_format='%.6e')

        #Add fields to sum
        mxsum += mx_z0
        mysum += my_z0
        wxsum += wx_z0
        wysum += wy_z0
        wzsum += wz_z0
        psum += p_z0
        uxsum += Ux_z0
        uysum += Uy_z0
        uzsum += Uz_z0
        #Add Bipolar
        pBsum += pBip
        uxBsum += UxBip
        uyBsum += UyBip
        countPlot += 1
        print('Transformed and Stored Data: idx = %i'%idx)
        sys.stdout.flush()
            
    #Calculate Averaged Fields
    mxavg = mxsum/(1.0*nSims)
    myavg = mysum/(1.0*nSims)
    wxavg = wxsum/(1.0*nSims)
    wyavg = wysum/(1.0*nSims)
    wzavg = wzsum/(1.0*nSims)
    pavg = psum/(1.0*nSims)
    uxavg = uxsum/(1.0*nSims)
    uyavg = uysum/(1.0*nSims)
    uzavg = uzsum/(1.0*nSims)
    #Bipolar
    #eta = np.linspace(-np.pi,np.pi,128)
    #xi_bar = np.linspace(0.0,1.0,128)
    r_1, r_2, d = 0.3, 0.15, 0.9
    a,meta,mxi,xi_1,xi_2 = CalculateBipolarMesh(r_1,r_2,d,NBP,NBP)
    h = a/(np.cosh(mxi) - np.cos(meta))
    myBavg = h*np.sinh(mxi)
    mxBavg = h*np.sin(meta)
    pBavg  = pBsum/(1.0*nSims)
    uxBavg = uxBsum/(1.0*nSims)
    uyBavg = uyBsum/(1.0*nSims)

    #Export Averaged Field Data
    #Create a dataframe containing flattened arrays for (mx, my, avgW, avgP, avgU)
    avgDict = {'mx':mxavg.flatten(),'my':myavg.flatten(),
               'avgWx':wxavg.flatten(),'avgWy':wyavg.flatten(),'avgWz':wzavg.flatten(),
               'avgP':pavg.flatten(),'avgUx':uxavg.flatten(),'avgUy':uyavg.flatten(),
               'avgUz':uzavg.flatten()}
    avgFieldData = pd.DataFrame(data=avgDict)
    #Export Avg Field Data to .csv file
    pathlib.Path(cwd_DATA+'/AVG/').mkdir(parents=True, exist_ok=True)
    avgFieldData.to_csv(cwd_DATA+'/AVG/AVG_XY_%04d.csv'%countPer,index=False,sep=' ',float_format='%.6e')
    #Bipolar
    #Create a dataframe containing flattened arrays for (mx, my, avgW, avgP, avgU)
    avgBDict = {'mx':mxBavg.flatten(),'my':myBavg.flatten(),
                'avgP':pBavg.flatten(),'avgUx':uxBavg.flatten(),'avgUy':uyBavg.flatten()}
    avgBFieldData = pd.DataFrame(data=avgBDict)
    #Export Avg Field Data to .csv file
    pathlib.Path(cwd_DATA+'/AVG/').mkdir(parents=True, exist_ok=True)
    avgBFieldData.to_csv(cwd_DATA+'/AVG/AVGB_XY_%04d.csv'%countPer,index=False,sep=' ',float_format='%.6e')
        
    stend = t.clock()
    diff = stend - start
    print('Time to run for 1 period = %.5fs'%diff)
    sys.stdout.flush()
    #countPer += 1

