#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 13:54:26 2019

@author: thomas
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#This script will be used to compare 2 time steps of
#flow tank data. We will have inputs (U_x,U_y,U_z) from 
#time step t and t + 0.25 s
#The input files will be in the directory path 
#'../VelocityFields/time'+str(timeValue)+'s/' 
#
#Input data comes from a uniform grid of the flow tank
#The script will output the relative error for each cell
#The final output will be an average of all the cells (1 number)
#If the average relative error < 5%, we have reached steady state!
#
#The following code will be specific to using data from time steps at 4.00s and 4.25s
#The code will be generalized after it is shown to be working properly
#

#CONSTANTS
    
#Obtain python directory    
cwd_PYTHON = os.getcwd()
TINY_NUMBER = 1.0e-15
FIGNUM = 1

def getInputData(timeValue):
    #Obtain directory of data for timeValue
    cwd_DATA = cwd_PYTHON + '/'+'../VelocityFields/time'+timeValue+'s/ZNormal/'
    #Create an array for files U_x, U_y, U_z
    #arrays will be of size (1024,256) where it is the data of coords (y,z) along plane x=0
    #U[0,0] -> (-40,160)
    arrayU_x = np.loadtxt(cwd_DATA+'U_x_clean.txt')
    arrayU_y = np.loadtxt(cwd_DATA+'U_y_clean.txt')
    arrayU_z = np.loadtxt(cwd_DATA+'U_z_clean.txt')
    print('time = ',timeValue)
    print(arrayU_y.shape)
    print(arrayU_y[1023,1])
    #Reverse the order of y coord data
    #(y,z)
    #U[0,0] -> (-40,-160)
    flipU_x = np.flip(arrayU_x,0)
    flipU_y = np.flip(arrayU_y,0)
    flipU_z = np.flip(arrayU_z,0)
    print(flipU_y[0,1])
    #Transpose the array to get y on x-axis and z on y-axis for plotting
    #(z,y)
    #U[0,0] -> (-160,-40)
    TflipU_x = flipU_x.T
    TflipU_y = flipU_y.T
    TflipU_z = flipU_z.T
    print(TflipU_y.shape)
    print(TflipU_y[1,0])
    return (TflipU_x, TflipU_y, TflipU_z)

def calcRelativeError():
    
    return

def PlotVelData(data,tVL,plotName):
    global FIGNUM
    #View Data with imshow
    #View U_? data
    nTime=len(tVL)
    fig = plt.figure(num=FIGNUM,figsize=(8,3),dpi=200)
    for idxTime in range(nTime):
        U = data[idxTime]
        ax = fig.add_subplot(1,1,idxTime+1)
        ax.set_title('Z= 10.0 mm Slice: time = '+tVL[idxTime] + ' s')
        if(plotName == 'U_y' or plotName == 'U_mag'):
            vminY = 0.0
            vmaxY = 70.0
            img = ax.imshow(U,interpolation='nearest',cmap='viridis',
                       extent=[-160,160,-40,40],
                       vmin = vminY, vmax = vmaxY,
                       origin='lower',aspect='equal')
        else:
            img = ax.imshow(U,interpolation='nearest',cmap='viridis',
                       extent=[-160,160,-40,40],
                       origin='lower',aspect='equal')
        ax.set_xlabel('Y-Axis (mm)')
        ax.set_ylabel('X-Axis (mm)')
        fig.colorbar(img,ax=ax)
    fig.tight_layout()
    fig.savefig('../Figures/ZNorm_'+plotName+'.png')
    FIGNUM += 1
    fig.clf()
    
def PlotRelError(data,nDim):
    global FIGNUM
    fig = plt.figure(num=FIGNUM,figsize=(8,6),dpi=200)
    for idxDim in range(nDim):
        ax = fig.add_subplot(3,1,idxDim+1)
        Err = data[idxDim]
        ax.set_title('Relative Diff for U[%i]'%(idxDim))
        img = ax.imshow(Err,interpolation='nearest',cmap='viridis',
                        vmin=0.01,vmax=10.0,
                       extent=[-160,160,-40,40],
                       origin='lower',aspect='equal')
        fig.colorbar(img,ax=ax)
    fig.tight_layout()
    fig.savefig('../Figures/err.png')
    FIGNUM += 1
    fig.clf()
    
    return

def main():
    global FIGNUM
    #List of time values to be compared (strings)
    timeValueList = ['8.5']
    nTime = len(timeValueList)
    nDim = 3
    nY = 1024
    nZ = 256  
    #Allocate Array for Velocity Data
    #there will be 1 extra dim for time. It will be (t=4,t=4.25,error)
    #there will be 1 extra dim for space (x, y, z, mag)
    dataU = np.zeros((nTime+1,nDim+1,nZ,nY))
    #Loop over time and store U_x,U_y,U_z data
    for idxTime in range(nTime):
        timeValue = timeValueList[idxTime]
        dataU[idxTime,0], dataU[idxTime,1], dataU[idxTime,2] = getInputData(timeValue)
        dataU[idxTime,3] = np.sqrt(dataU[idxTime,0]**2 + dataU[idxTime,1]**2 + dataU[idxTime,2]**2)
        #print(dataU[idxTime,3,0,1])
        
    #View Data with imshow
    plotNameList = ['U_x','U_y','U_z','U_mag']
    for idxDim in range(len(plotNameList)):
        pN = plotNameList[idxDim]
        PlotVelData(dataU[:,idxDim],timeValueList,pN)
    
    #Data for both time steps is now stored in dataU
    return
    '''#Calculate the error for each cell
    #calculations will be stored in dataU[2]
    dataU[2] = abs(dataU[1] - dataU[0])/(dataU[1]+TINY_NUMBER)
    #1) Calculate out the i) maximum error and ii) the average error
    maxError = np.amax(dataU[2])
    avgError = np.mean(dataU[2])
    print('max error: %.5e \t avg error: %.5e\n'%(maxError,avgError))
    #Calculate each cell 1 by 1
    for idxDim in range(nDim):
        for idxY in range(nY):
            for idxZ in range(nZ):
                #dataU[2,idxDim,idxZ,idxY] = abs(dataU[1,idxDim,idxZ,idxY] - dataU[0,idxDim,idxZ,idxY])/(dataU[1,idxDim,idxZ,idxY]+TINY_NUMBER)
                dataU[2,idxDim,idxZ,idxY] = abs(dataU[1,idxDim,idxZ,idxY] - dataU[0,idxDim,idxZ,idxY])
        maxError = np.amax(dataU[2,idxDim])
        avgError = np.mean(dataU[2,idxDim])
        print('max error: %.5e \t avg error: %.5e\n'%(maxError,avgError))
        
    #View Error with imshow
    PlotRelError(dataU[2],nDim)'''
    
    
    
    
    
    
    return

#------------------------__END_MAIN__------------------------------------
main()