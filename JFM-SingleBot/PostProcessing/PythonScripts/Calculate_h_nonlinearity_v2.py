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
DX = 2.0/128.0
RHO = 2.0
FREQ = 10.0
OMEGA = 2.0*np.pi*FREQ
RADIUSLARGE = 0.3


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
    PList  = fieldData['avgP'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny, Nz = 128, 512, 128
    mxArr = np.array(mxList).reshape((Nx,Ny,Nz))
    myArr = np.array(myList).reshape((Nx,Ny,Nz))
    mzArr = np.array(mzList).reshape((Nx,Ny,Nz))
    UxArr = np.array(UxList).reshape((Nx,Ny,Nz))
    UyArr = np.array(UyList).reshape((Nx,Ny,Nz))
    UzArr = np.array(UzList).reshape((Nx,Ny,Nz))
    PArr  = np.array(PList).reshape((Nx,Ny,Nz))
    return (mxArr, myArr, mzArr, UxArr, UyArr, UzArr,PArr)

def GetHatUFieldData(cwd,idx):
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    #cwd = cwd_PYTHON
    fieldData = pd.read_csv(cwd+'allhatU_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    mzList = fieldData['mz'].values.tolist()
    UxcList = fieldData['hatUxc'].values.tolist()
    UxsList = fieldData['hatUxs'].values.tolist()
    UycList = fieldData['hatUyc'].values.tolist()
    UysList = fieldData['hatUys'].values.tolist()
    UzcList = fieldData['hatUzc'].values.tolist()
    UzsList = fieldData['hatUzs'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny, Nz = 128, 512, 128
    mxArr = np.array(mxList).reshape((Nx,Ny,Nz))
    myArr = np.array(myList).reshape((Nx,Ny,Nz))
    mzArr = np.array(mzList).reshape((Nx,Ny,Nz))
    UxcArr = np.array(UxcList).reshape((Nx,Ny,Nz))
    UxsArr = np.array(UxsList).reshape((Nx,Ny,Nz))
    UycArr = np.array(UycList).reshape((Nx,Ny,Nz))
    UysArr = np.array(UysList).reshape((Nx,Ny,Nz))
    UzcArr = np.array(UzcList).reshape((Nx,Ny,Nz))
    UzsArr = np.array(UzsList).reshape((Nx,Ny,Nz))
    return (mxArr, myArr, mzArr, UxcArr, UxsArr, UycArr, UysArr, UzcArr, UzsArr)

# caculates the gradient of a grid of values, provided an
# isotropic grid spacing h. Throws away outer skin of points
def grad(mat,h=1):
    gy = 0.5*(np.roll(mat,-1,axis=1) - np.roll(mat,1,axis=1))/h
    gx = 0.5*(np.roll(mat,-1,axis=0) - np.roll(mat,1,axis=0))/h
    gz = 0.5*(np.roll(mat,-1,axis=2) - np.roll(mat,1,axis=2))/h
    return (gx[1:-1,1:-1,1:-1],gy[1:-1,1:-1,1:-1],gz[1:-1,1:-1,1:-1])

# calculates the laplacian of a grid of values, provided
# an isotropic grid spacing h. Throws away outer skin of points
def lapl(mat,h=1):
    l = -6*mat + np.roll(mat,-1,axis=0) + np.roll(mat,1,axis=0) +\
        np.roll(mat,-1,axis=1) + np.roll(mat,1,axis=1) +\
        np.roll(mat,-1,axis=2) + np.roll(mat,1,axis=2)
    return l[1:-1,1:-1,1:-1]/(h*h)


### MAIN SCRIPT
#### READ ALL VTK FILES IN A SIMULATION DIRECTORY
#### CALCULATE AVERAGE FIELD DATA
#### EXPORT AVERAGED FIELD DATA TO CSV

if __name__ == '__main__':
    
    ####READ ALL VTK FILES
    
    start = t.clock()
    cwd_SIM = cwd_PYTHON+'../'+Par+'/s'+ParVal+'/'
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
    avgmx,avgmy,avgmz,avgUx,avgUy,avgUz,avgP = GetAvgFieldData(cwd_AVG,countPer)
    ###Obtain hatU Field Data for Per countPer
    hatmx,hatmy,hatmz,hatUxc,hatUxs,hatUyc,hatUys,hatUzc,hatUzs = GetHatUFieldData(cwd_AVG,countPer)
    
    udotgradusum_x = np.zeros((126,510,126))
    udotgradusum_y = np.zeros((126,510,126))
    udotgradusum_z = np.zeros((126,510,126))
    ###Obtain U field data for Per countPer
    for idx in range(0,nSims):
        vtkStart = t.clock()
        dumpIdx = nSims*countPer+idx
        #Read All VTK Data into List
        allData = ReadAllVTK(cwd_DATA,dumpIdx)
        #Extract important values from allData
        time,mx,my,mz,wxArr,wyArr,wzArr,pArr,uxArr,uyArr,uzArr, = ReadVTK(allData)
        #Calculate gradient of u, then find vdotgradv
        gradUx_x, gradUx_y, gradUx_z = grad(uxArr,h=DX)
        gradUy_x, gradUy_y, gradUy_z = grad(uyArr,h=DX)
        gradUz_x, gradUz_y, gradUz_z = grad(uzArr,h=DX)
        uxdotgradu = uxArr[1:-1,1:-1,1:-1]*gradUx_x + uyArr[1:-1,1:-1,1:-1]*gradUx_y + uzArr[1:-1,1:-1,1:-1]*gradUx_z
        uydotgradu = uxArr[1:-1,1:-1,1:-1]*gradUy_x + uyArr[1:-1,1:-1,1:-1]*gradUy_y + uzArr[1:-1,1:-1,1:-1]*gradUy_z
        uzdotgradu = uxArr[1:-1,1:-1,1:-1]*gradUz_x + uyArr[1:-1,1:-1,1:-1]*gradUz_y + uzArr[1:-1,1:-1,1:-1]*gradUz_z
        #Add to the sum
        udotgradusum_x += uxdotgradu
        udotgradusum_y += uydotgradu
        udotgradusum_z += uzdotgradu

    #Calculate Avg vdotgradv
    udotgraduavg_x = udotgradusum_x/(1.0*nSims)
    udotgraduavg_y = udotgradusum_y/(1.0*nSims)
    udotgraduavg_z = udotgradusum_z/(1.0*nSims)
    #Calculate nonlinear g term with instantaneous flow
    g2x = -1.0*RHO*udotgraduavg_x
    g2y = -1.0*RHO*udotgraduavg_y
    g2z = -1.0*RHO*udotgraduavg_z
    
    ###Calculate Gradients and Laplacians for avgU, avgP, and hatU
    ##lap(U) - del(P)
    lapUx = lapl(avgUx,DX)
    lapUy = lapl(avgUy,DX)
    lapUz = lapl(avgUz,DX)
    gradP = grad(avgP,h=DX)
    ##hatU dot (del(hatU_star))
    gradhatUxc_x,gradhatUxc_y,gradhatUxc_z = grad(hatUxc,DX)
    gradhatUxs_x,gradhatUxs_y,gradhatUxs_z = grad(hatUxs,DX)
    gradhatUyc_x,gradhatUyc_y,gradhatUyc_z = grad(hatUyc,DX)
    gradhatUys_x,gradhatUys_y,gradhatUys_z = grad(hatUys,DX)
    gradhatUzc_x,gradhatUzc_y,gradhatUzc_z = grad(hatUzc,DX)
    gradhatUzs_x,gradhatUzs_y,gradhatUzs_z = grad(hatUzs,DX)
    ###Calculate by dimension
    #Msquared = 2.0*float(ParVal)
    MU = RHO*OMEGA*RADIUSLARGE**(2)/(2.0*float(ParVal))
    gx = gradP[0] - MU*lapUx
    fx = -0.5*RHO*(hatUxc[1:-1,1:-1,1:-1]*gradhatUxc_x - hatUxs[1:-1,1:-1,1:-1]*gradhatUxs_x +
                   hatUyc[1:-1,1:-1,1:-1]*gradhatUxc_y - hatUys[1:-1,1:-1,1:-1]*gradhatUxs_y +
                   hatUzc[1:-1,1:-1,1:-1]*gradhatUxc_z - hatUzs[1:-1,1:-1,1:-1]*gradhatUxs_z)
    #fx = -0.5*Msquared*(np.dot(hatUxc[1:-1,1:-1,1:-1],gradhatUxc) + np.dot(hatUxs[1:-1,1:-1,1:-1],gradhatUxs))
    gy = gradP[1] - MU*lapUy
    fy = -0.5*RHO*(hatUxc[1:-1,1:-1,1:-1]*gradhatUyc_x - hatUxs[1:-1,1:-1,1:-1]*gradhatUys_x +
                   hatUyc[1:-1,1:-1,1:-1]*gradhatUyc_y - hatUys[1:-1,1:-1,1:-1]*gradhatUys_y +
                   hatUzc[1:-1,1:-1,1:-1]*gradhatUyc_z - hatUzs[1:-1,1:-1,1:-1]*gradhatUys_z)
    #fy = -0.5*Msquared*(np.dot(hatUyc[1:-1,1:-1,1:-1],gradhatUyc) + np.dot(hatUys[1:-1,1:-1,1:-1],gradhatUys))
    gz = gradP[2] - MU*lapUz
    fz = -0.5*RHO*(hatUxc[1:-1,1:-1,1:-1]*gradhatUzc_x - hatUxs[1:-1,1:-1,1:-1]*gradhatUzs_x +
                   hatUyc[1:-1,1:-1,1:-1]*gradhatUzc_y - hatUys[1:-1,1:-1,1:-1]*gradhatUzs_y +
                   hatUzc[1:-1,1:-1,1:-1]*gradhatUzc_z - hatUzs[1:-1,1:-1,1:-1]*gradhatUzs_z)
    #fz = -0.5*Msquared*(np.dot(hatUzc[1:-1,1:-1,1:-1],gradhatUzc) + np.dot(hatUzs[1:-1,1:-1,1:-1],gradhatUzs))
    hx = fx - gx
    hy = fy - gy
    hz = fz - gz
    h2x = fx - g2x
    h2y = fy - g2y
    h2z = fz - g2z
    
    nonlinearDict = {'fx':fx.flatten(),'fy':fy.flatten(),'fz':fz.flatten(),
                     'gx':gx.flatten(),'gy':gy.flatten(),'gz':gz.flatten(),
                     'hx':hx.flatten(),'hy':hy.flatten(),'hz':hz.flatten(),
                     'mx':hatmx[1:-1,1:-1,1:-1].flatten(),'my':hatmy[1:-1,1:-1,1:-1].flatten(),'mz':hatmz[1:-1,1:-1,1:-1].flatten(),
                     'lapux':lapUx.flatten(),'lapuy':lapUy.flatten(),'lapuz':lapUz.flatten(),
                     'gradpx':gradP[0].flatten(),'gradpy':gradP[1].flatten(),'gradpz':gradP[2].flatten(),
                     'g2x':g2x.flatten(),'g2y':g2y.flatten(),'g2z':g2z.flatten(),
                     'h2x':h2x.flatten(),'h2y':h2y.flatten(),'h2z':h2z.flatten()
                    }
    nonlinearData = pd.DataFrame(data=nonlinearDict)
    #Export Avg Field Data to .csv file
    pathlib.Path(cwd_DATA+'/AVG/').mkdir(parents=True, exist_ok=True)
    #nonlinearData.to_csv(cwd_DATA+'/AVG/h_nonlinear_%04d.csv'%countPer,index=False,sep=' ',float_format='%.5e')
    nonlinearData.to_csv(cwd_DATA+'/AVG/h_nonlinear_v2.csv',index=False,sep=' ',float_format='%.5e')

    stend = t.clock()
    diff = stend - start
    print('Time to run for 1 period = %.5fs'%diff)
    sys.stdout.flush()
 
