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

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
parTheta = sys.argv[1]
parHx = sys.argv[2]
parHy = sys.argv[3]

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
    #Obtain Pressure
    startIdx = endIdx+5
    endIdx = startIdx+Nlines
    pList = data[startIdx:endIdx]
    pFlat = [item for sublist in pList for item in sublist]
    pVals = [float(i) for i in pFlat]
    pArr = np.array(pVals)
    pArr = pArr.reshape((Nx,Ny))
    #Field Data
    Ndim = int(data[endIdx+1][1])
    Nlines = int(np.trunc(Npts*Ndim/9.0)+1.0)
    #Obtain Omega
    startIdx = endIdx + 2
    endIdx = startIdx + Nlines
    wList = data[startIdx:endIdx]
    wFlat = [item for sublist in wList for item in sublist]
    wVals = [float(i) for i in wFlat]
    wArr = np.array(wVals)
    wArr = wArr.reshape((Nx,Ny))
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

# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    return cwd+"pd.txt"

def GetPosDataLength(cwd):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    UAdata = data[data['idx'] == 6].copy()
    return len(UAdata['time'])
    
### MAIN SCRIPT
#### READ ALL VTK FILES IN A SIMULATION DIRECTORY
#### CALCULATE AVERAGE FIELD DATA
#### EXPORT AVERAGED FIELD DATA TO CSV

if __name__ == '__main__':
    #Simulation Parameters
    cwd_Re = cwd_PYTHON+'../Fluid/Theta{0}/Hx{1}/Hy{2}/'.format(parTheta,parHx,parHy)
    #Calculate # Periods
    DUMP_INT = 20.0
    DT, PERIOD = 0.005, 0.1
    cwd_POS = cwd_Re
    nTime = GetPosDataLength(cwd_POS)
    nPer = int(np.trunc(1.0*nTime/DUMP_INT))
    #nPer = 1
    #Paths to data and plots
    cwd_DATA = cwd_Re+'VTK/'

    for countPer in range(nPer):
        AVGfile = pathlib.Path(cwd_DATA+'AVG/AVG_%04d.csv'%countPer)
        if not AVGfile.exists ():
            start = t.clock()
            #Create sum arrays
            Nx, Ny = 1024, 1024
            mxsum = np.zeros((Nx,Ny))
            mysum = np.zeros((Nx,Ny))
            wsum = np.zeros((Nx,Ny))
            psum = np.zeros((Nx,Ny))
            uxsum = np.zeros((Nx,Ny))
            uysum = np.zeros((Nx,Ny))
            #Create posData sum database
            for idx in range(0,int(DUMP_INT)):
                dumpIdx = int(DUMP_INT)*countPer+idx
                #Read All VTK Data into List
                allData = ReadAllVTK(cwd_DATA,dumpIdx)
                #Extract important values from allData
                time, mx, my, wArr, pArr, uxArr, uyArr = ReadVTK(allData)
                print('time = ',time)
                sys.stdout.flush()
                #Add fields to sum
                mxsum += mx
                mysum += my
                wsum += wArr
                psum += pArr
                uxsum += uxArr
                uysum += uyArr
            #Calculate Averaged Fields
            mxavg = mxsum/DUMP_INT
            myavg = mysum/DUMP_INT
            wavg = wsum/DUMP_INT
            pavg = psum/DUMP_INT
            uxavg = uxsum/DUMP_INT
            uyavg = uysum/DUMP_INT
            #Export Averaged Field Data
            #Create a dataframe containing flattened arrays for (mx, my, avgW, avgP, avgU)
            avgDict = {'mx':mxavg.flatten(),'my':myavg.flatten(),'avgW':wavg.flatten(),
                        'avgP':pavg.flatten(),'avgUx':uxavg.flatten(),'avgUy':uyavg.flatten()}
            avgFieldData = pd.DataFrame(data=avgDict)
            #Export Avg Field Data to .csv file
            pathlib.Path(cwd_DATA+'AVG/').mkdir(parents=True, exist_ok=True)
            avgFieldData.to_csv(cwd_DATA+'AVG/AVG_%04d.csv'%countPer,index=False,sep=' ',float_format='%.5e')
            
            stend = t.clock()
            diff = stend - start
            print('Time to run for 1 period = %.5fs'%diff)
            sys.stdout.flush()
        else:
            print('Re = %s: Per = %i: AVG file exists already'%(Re,countPer))

