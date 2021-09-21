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
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)
import pathlib
from matplotlib.colors import Normalize
from scipy import interpolate
norm = Normalize()
from resource import getrusage, RUSAGE_SELF

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
RHO = 1000.0
DX = 0.025/256.0
PERIOD = 0.1
FREQUENCY = 1.0/PERIOD
OMEGA = 2.0*np.pi*FREQUENCY
RADIUS_LARGE = 0.002
AMPLITUDE = 0.8*RADIUS_LARGE
EPSILON = AMPLITUDE/AMPLITUDE
maxWin = 0.03
minWin = -1.0*maxWin
config = 'single'
Re = sys.argv[1]

csfont = {'fontname':'Times New Roman'}

# constructs a filepath for the plot omages of Re=$Re, config=$config, and field=$field
def plotName(cwd,Re,config,field,idx):
    #field = Vort, Pres, Vel, AvgW, AvgP, AvgU
    strDir = cwd+"../Figures/VelDecay/{0}/".format(config)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return strDir+"{0}_{1}_{2}_{3}".format(config,Re,field,idx)

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

def Rotate(xy, theta):
    # https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions
    #First Rotate based on Theta
    #Allocate Arrays
    rotationMatrix = np.zeros((2,2))
    #Calculate rotation matrix
    rotationMatrix[0,0] = np.cos(theta)
    rotationMatrix[0,1] = -1.0*np.sin(theta)
    rotationMatrix[1,0] = np.sin(theta)
    rotationMatrix[1,1] = np.cos(theta)
    return rotationMatrix.dot(xy)

def CalcLabAngle(pos):
    #Find swimming axis (normal y-axis)
    xU, xL = pos.loc[0,'xU'], pos.loc[0,'xL']
    yU, yL = pos.loc[0,'yU'], pos.loc[0,'yL']
    labX   = xU - xL
    labY   = yU - yL
    length = np.hypot(labX,labY)
    normX = labX/length
    normY = labY/length
    #2) Calculate Theta
    if(normX <= 0.0):
        theta = np.arccos(normY)
    else:
        theta = -1.0*np.arccos(normY)+2.0*np.pi
    print('theta = ',theta*180.0/np.pi)
    return 2.0*np.pi - theta

def RotateVectorField(pos,mx,my,Ux,Uy,NX,NY):
    #Shift field by CM
    #Calculate angle of swimmer 1 from y-axis
    #Rotate field by 2pi - theta
    #Shift x and y by the CM location
    xCM = 0.8*pos.loc[0,'xU'] + 0.2*pos.loc[0,'xL']
    yCM = 0.8*pos.loc[0,'yU'] + 0.2*pos.loc[0,'yL']
    #Do the same for mx and my
    mx -= xCM
    my -= yCM
    #Shift pos data by xCM and yCM
    pos['xU'] -= xCM
    pos['xL'] -= xCM
    pos['yU'] -= yCM
    pos['yL'] -= yCM
    #Rotate Reference frame by swimmer 1's axis
    #Calculate Theta (Rotate by -Theta)
    theta_rotate = CalcLabAngle(pos)
    print('theta_rotate = ',theta_rotate*180.0/np.pi)
    mxy = np.array([mx.flatten(),my.flatten()])
    mxy_rot = np.zeros((2,NX*NY))
    #Do the same for the U field
    Uxy = np.array([Ux.flatten(),Uy.flatten()])
    Uxy_rot = np.zeros((2,NX*NY))
    for jdx in range(NX*NY):
        mxy_rot[:,jdx] = Rotate(mxy[:,jdx],theta_rotate)
        Uxy_rot[:,jdx] = Rotate(Uxy[:,jdx],theta_rotate)
    mx_rot = mxy_rot[0,:].reshape((NX,NY))
    my_rot = mxy_rot[1,:].reshape((NX,NY))
    Ux_rot = Uxy_rot[0,:].reshape((NX,NY))
    Uy_rot = Uxy_rot[1,:].reshape((NX,NY))
    aU_pos = np.array([pos.loc[0,'xU'],pos.loc[0,'yU']])
    aL_pos = np.array([pos.loc[0,'xL'],pos.loc[0,'yL']])
    aU_rot = Rotate(aU_pos,theta_rotate)
    print('aU = ',aU_pos)
    print('aU_rot = ',aU_rot)
    aL_rot = Rotate(aL_pos,theta_rotate)
    pos['aXU_rot'], pos['aYU_rot'] = aU_rot[0], aU_rot[1]
    pos['aXL_rot'], pos['aYL_rot'] = aL_rot[0], aL_rot[1]
    
    return (pos,mx_rot,my_rot,Ux_rot,Uy_rot)

def InterpolateToNewCoordinateSystem(mx,my,arrayUx,arrayUy,mx_new,my_new):
    
    #Interpolate Ux and Uy from original cartesian coordainates to new ones
    #Griddata
    print('About to inteprolate field data')
    print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)
    sys.stdout.flush()
    arrayUx_new=interpolate.griddata((mx.flatten(),my.flatten()),arrayUx.flatten() , (mx_new,my_new),method='linear')
    print('X transformation complete')
    print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)
    sys.stdout.flush()
    arrayUy_new=interpolate.griddata((mx.flatten(),my.flatten()),arrayUy.flatten() , (mx_new,my_new),method='linear')
    print('Coordinate Transformation Complete!')
    print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)
    sys.stdout.flush()
    return (arrayUx_new,arrayUy_new)

def PlotVelocityDecay(cwd,time,mx,my,w,Ux,Uy,pos,pars):
    global OMEGA, EPSILON, RADIUS_LARGE
    Re = pars[0]
    config = pars[1]
    field = pars[2]
    idx = pars[3]
    fig, ax = plt.subplots(nrows=2, ncols=2,figsize=(12,12),dpi=200,num=10)
    ax[0,0].set_ylabel(r'$U_x$ (U/fR)',fontsize=25,**csfont)
    ax[1,0].set_ylabel(r'$U_y$ (U/fR)',fontsize=25,**csfont)
    ax[1,0].set_xlabel(r'$r$ (R)',fontsize=25,**csfont)
    ax[1,1].set_xlabel(r'$r$ (R)',fontsize=25,**csfont)
    #ax.set_title(r'time = %.2fs'%(time),fontsize=16)

    #Rotate by Swimmer 1 axis around CM of pair
    pos,mx_rot,my_rot,Ux_rot,Uy_rot = RotateVectorField(pos,mx,my,Ux,Uy,1022,1022)

    #Interpolate onto a new coordinate system
    x = np.linspace(-0.05,0.05,512)
    y = np.linspace(-0.05,0.05,512)
    mx_stream, my_stream = np.meshgrid(x,y)
    interpUx, interpUy = InterpolateToNewCoordinateSystem(mx_rot,my_rot,Ux_rot,Uy_rot,mx_stream,my_stream)
    
    #Now that we have the interpolated (Rotated) velocity field
    #We can calculate the lineouts for each dir (+-x, +-y)
    #First, get lineouts from mesh
    #For neg, loop over 0:256, get avg of 255,256 values for Ux or Uy
    #For pos, loop over 256:512, get avg of 255,256 values for Ux or Uy
    line_x_neg = mx_stream[0,0:256]/RADIUS_LARGE
    line_x_pos = mx_stream[0,256:512]/RADIUS_LARGE
    print(line_x_pos)
    sys.stdout.flush()
    line_y_neg = my_stream[0:256,0]/RADIUS_LARGE
    line_y_pos = my_stream[256:512,0]/RADIUS_LARGE
    U_x_neg = 0.5*(interpUx[255,0:256] + interpUx[256,0:256])/(FREQUENCY*RADIUS_LARGE)
    U_y_neg = 0.5*(interpUy[0:256,255] + interpUy[0:256,256])/(FREQUENCY*RADIUS_LARGE)
    U_x_pos = 0.5*(interpUx[255,256:512] + interpUx[256,256:512])/(FREQUENCY*RADIUS_LARGE)
    U_y_pos = 0.5*(interpUy[256:512,255] + interpUy[256:512,256])/(FREQUENCY*RADIUS_LARGE)
    
    #Plot lineout for each axis direction
    #Plot y = 0
    ax[0,0].plot([np.amin(line_x_neg),np.amax(line_x_neg)],[0.0,0.0],c='k')
    ax[0,1].plot([np.amin(line_x_pos),np.amax(line_x_pos)],[0.0,0.0],c='k')
    ax[1,0].plot([np.amin(line_y_neg),np.amax(line_y_neg)],[0.0,0.0],c='k')
    ax[1,1].plot([np.amin(line_y_pos),np.amax(line_y_pos)],[0.0,0.0],c='k')
    #Plot Velocities
    ax[0,0].plot(line_x_neg,U_x_neg,c='k')
    ax[0,1].plot(line_x_pos,U_x_pos,c='k')
    ax[1,0].plot(line_y_neg,U_y_neg,c='k')
    ax[1,1].plot(line_y_pos,U_y_pos,c='k')
    
    #axis = [- 0.025,0.025,-0.025,0.025]
    #ax.axis(axis)
    #ax.set_aspect('equal')
    
    ax[0,0] = SetTickParams(ax[0,0])
    ax[0,1] = SetTickParams(ax[0,1])
    ax[1,0] = SetTickParams(ax[1,0])
    ax[1,1] = SetTickParams(ax[1,1])
    #ax[0,0].yaxis.set_minor_locator(MultipleLocator(0.0025))
    #ax[0,0].yaxis.set_major_locator(MultipleLocator(0.005))
    #ax[0,1].yaxis.set_minor_locator(MultipleLocator(0.0025))
    #ax[0,1].yaxis.set_major_locator(MultipleLocator(0.005))
    #ax[1,0].yaxis.set_minor_locator(MultipleLocator(0.04))
    #ax[1,0].yaxis.set_major_locator(MultipleLocator(0.02))
    #ax[1,1].yaxis.set_minor_locator(MultipleLocator(0.005))
    #ax[1,1].yaxis.set_major_locator(MultipleLocator(0.0025))
    for jdx in range(2):
        for kdx in range(2):
            for label in (ax[jdx,kdx].get_xticklabels() + ax[jdx,kdx].get_yticklabels()):
                label.set_fontsize(20)

    #ax[0,0] = set_size(6,6,ax[0,0])
    #ax[0,1] = set_size(6,6,ax[0,1])
    #ax[1,0] = set_size(6,6,ax[1,0])
    #ax[1,1] = set_size(6,6,ax[1,1])
    fig.tight_layout()

    fig.savefig(plotName(cwd,Re,config,field,idx)+'.png')
    fig.clf()
    plt.close()
    return

def SetTickParams(ax):
    ax.tick_params(which='major',axis='both',direction='in',length=14,width=1,zorder=10)
    ax.tick_params(which='minor',axis='both',direction='in',length=8,width=0.75)
    ax.xaxis.set_major_locator(MultipleLocator(5.0))
    ax.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=2))
    return ax

# constructs a filepath for the pos data of Re = $Re
def pname(cwd,Re):
    #return cwd+"/pd.txt"
    #cwd = cwd_PYTHON
    return cwd+"/pd_Re{0}.txt".format(Re)

def GetPosData(cwd,time,Re):
    data = pd.read_csv(pname(cwd,Re),delimiter=' ')
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
    cwd_Re = cwd_PYTHON+'../Fig4_VisitFiles_Single/Re'+Re+'/'
    #Extract Position Data
    cwd_POS = cwd_POS = cwd_PYTHON+'../PosData/'
    #Calculate # Periods
    DUMP_INT = 20.0
    #nTime = GetPosDataLength(cwd_POS)
    nPer = 200
    #nPer = 2
    #Paths to data and plots
    cwd_DATA = cwd_Re+'/VTK/AVG/'
    countPer = 0
    for countPer in range(nPer):
        if(countPer == 100):
            AVGPlot = pathlib.Path(cwd_DATA+'AVG_%04d.csv'%countPer)
            #if not AVGPlot.exists ():
            if AVGPlot.exists ():
                start = t.clock()
                #Get Avg Field Data
                mx,my,avgW,avgP,avgUx,avgUy = GetAvgFieldData(cwd_DATA,countPer)
                #Extract Position and Time Data
                time = np.round(0.05 + countPer*PERIOD,2)
                #print('time = ',time)
                posData = GetPosData(cwd_POS,time,Re)
                #print(posData)
                #Plot Averaged Field Data
                #Vorticity And Streamlines
                pars = [Re,'single','VelDecay',countPer]
                PlotVelocityDecay(cwd_PYTHON,time,mx,my,avgW,avgUx,avgUy,posData,pars)

                stend = t.clock()
                diff = stend - start
                print('Time to run for 1 period = %.5fs'%diff)
                sys.stdout.flush()

