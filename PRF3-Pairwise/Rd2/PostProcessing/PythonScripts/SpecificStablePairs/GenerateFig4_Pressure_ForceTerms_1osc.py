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
import copy
norm = Normalize()

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
RHO = 1000.0
DX = 0.025/256.0
PERIOD = 0.1
OMEGA = 2.0*np.pi/PERIOD
RADIUS_LARGE = 0.002
RADIUS_SMALL = 0.001
AMPLITUDE = 0.8*RADIUS_LARGE
AMPLITUDE_SMALL = 0.8*AMPLITUDE
maxWin = 0.03
minWin = -1.0*maxWin
config = sys.argv[1]
Re = sys.argv[2]
MU = RHO*OMEGA*AMPLITUDE_SMALL*RADIUS_SMALL/float(Re)

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

# constructs a filepath for the plot omages of Re=$Re, config=$config, and field=$field
def plotName(cwd,Re,config,field,idx):
    #field = Vort, Pres, Vel, AvgW, AvgP, AvgU
    strDir = cwd+"../Figures/AVG/{0}/".format(config)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    #return strDir+"{0}_{1}_{2}_{3}.png".format(config,Re,field,idx)
    #return strDir+"{0}_{1}_{2}_{3}.svg".format(config,Re,field,idx)
    return strDir+"{0}_{1}_{2}_{3}".format(config,Re,field,idx)

def AddDiscsToPlot(ax,pos):
    global RADIUS_LARGE, RADIUS_SMALL
    #Add Discs
    circle1 = Circle((pos.loc[0,'aXU_rot'], pos.loc[0,'aYU_rot']), RADIUS_LARGE, facecolor=(0.0,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle1)
    circle2 = Circle((pos.loc[0,'aXL_rot'], pos.loc[0,'aYL_rot']), RADIUS_SMALL, facecolor=(0.0,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle2)
    circle3 = Circle((pos.loc[0,'bXU_rot'], pos.loc[0,'bYU_rot']), RADIUS_LARGE, facecolor=(0.5,)*3,
                             linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle3)
    circle4 = Circle((pos.loc[0,'bXL_rot'], pos.loc[0,'bYL_rot']), RADIUS_SMALL, facecolor=(0.5,)*3,
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

def PlotPressure(cwd,time,mx,my,p,pos,pars):
    Re = pars[0]
    config = pars[1]
    field = pars[2]
    idx = pars[3]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200,num=10)
    #Pressure Contour
    #levels = MaxNLocator(nbins=21).tick_values(0.5*np.amin(p), -0.5*np.amin(p))
    levels = MaxNLocator(nbins=21).tick_values(-1.0,1.0)
    print('min P = ',0.5*np.amin(p))
    sys.stdout.flush()
    ax.contourf(mx,my,p,cmap='PuOr',levels=levels,extend='both')
    
    #Add Discs
    AddDiscsToPlot(ax,pos)
    xCM = 0.25*(pos.loc[0,'aXU'] + pos.loc[0,'bXU'] + pos.loc[0,'aXL'] + pos.loc[0,'bXL'])
    yCM = 0.25*(pos.loc[0,'aYU'] + pos.loc[0,'bYU'] + pos.loc[0,'aYL'] + pos.loc[0,'bYL'])
    axis = [xCM - 0.015,xCM+0.015,yCM-0.015,yCM+0.015]
    ax.axis(axis)
    ax.set_aspect('equal')
    # Turn off tick labels
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    fig.tight_layout()
    ax = set_size(6,6,ax)
    fig.savefig(plotName(cwd,Re,config,field,idx)+'.png')
    fig.clf()
    plt.close()
    return

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
    xU, xL = pos.loc[0,'aXU'], pos.loc[0,'aXL']
    yU, yL = pos.loc[0,'aYU'], pos.loc[0,'aYL']
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
    xCM = 0.25*(pos.loc[0,'aXU'] + pos.loc[0,'bXU'] + pos.loc[0,'aXL'] + pos.loc[0,'bXL'])
    yCM = 0.25*(pos.loc[0,'aYU'] + pos.loc[0,'bYU'] + pos.loc[0,'aYL'] + pos.loc[0,'bYL'])
    #Do the same for mx and my
    mx -= xCM
    my -= yCM
    #Shift pos data by xCM and yCM
    pos['aXU'] -= xCM
    pos['aXL'] -= xCM
    pos['bXU'] -= xCM
    pos['bXL'] -= xCM
    pos['aYU'] -= yCM
    pos['aYL'] -= yCM
    pos['bYU'] -= yCM
    pos['bYL'] -= yCM
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
    aU_pos = np.array([pos.loc[0,'aXU'],pos.loc[0,'aYU']])
    aL_pos = np.array([pos.loc[0,'aXL'],pos.loc[0,'aYL']])
    bU_pos = np.array([pos.loc[0,'bXU'],pos.loc[0,'bYU']])
    bL_pos = np.array([pos.loc[0,'bXL'],pos.loc[0,'bYL']])
    aU_rot = Rotate(aU_pos,theta_rotate)
    print('aU = ',aU_pos)
    print('aU_rot = ',aU_rot)
    aL_rot = Rotate(aL_pos,theta_rotate)
    bU_rot = Rotate(bU_pos,theta_rotate)
    bL_rot = Rotate(bL_pos,theta_rotate)
    pos['aXU_rot'], pos['aYU_rot'] = aU_rot[0], aU_rot[1]
    pos['aXL_rot'], pos['aYL_rot'] = aL_rot[0], aL_rot[1]
    pos['bXU_rot'], pos['bYU_rot'] = bU_rot[0], bU_rot[1]
    pos['bXL_rot'], pos['bYL_rot'] = bL_rot[0], bL_rot[1]
    
    return (pos,mx_rot,my_rot,Ux_rot,Uy_rot)

def PlotForceDensity(cwd,time,mx,my,forcex,forcey,pos,pars):
    FORCETOL = 1.0e-5
    SMALL_NUM = 1.0e-25
    Re = pars[0]
    config = pars[1]
    field = pars[2]
    print('field = ',field)
    sys.stdout.flush()
    idx = pars[3]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,8),dpi=200,num=10)
    #Force Contours
    forcemag = np.hypot(forcex,forcey)
    #Rotate by Swimmer 1 axis around CM of pair
    pos,mx_rot,my_rot,forcex_rot,forcey_rot = RotateVectorField(pos,mx,my,forcex,forcey,1022,1022)
    
    forcemag = np.where(forcemag == 0.0, SMALL_NUM,forcemag) #Avoid undefined log numbers
    #levels = MaxNLocator(nbins=21).tick_values(0.0, 0.5*np.amax(forcemag))
    levels = MaxNLocator(nbins=21).tick_values(-2.0, 3.0)
    ax.contourf(mx_rot,my_rot,np.log10(forcemag),cmap='viridis',levels=levels,extend='both')
    
    #Add quiver to magnitude plot
    normFx, normFy = np.zeros((1022,1022)), np.zeros((1022,1022))
    print('forcemagmax = ',np.amax(forcemag))
    forcemag = np.where(forcemag/np.amax(forcemag) <= FORCETOL, SMALL_NUM,forcemag)
    normFx = np.where(forcemag == SMALL_NUM, 0.0, forcex_rot/forcemag)
    normFy = np.where(forcemag == SMALL_NUM, 0.0, forcey_rot/forcemag)
    #normFx = forcex_rot/forcemag
    #normFy = forcey_rot/forcemag
    space=8
    ax.quiver(mx_rot[::space,::space],my_rot[::space,::space],
              normFx[::space,::space],normFy[::space,::space],
              color='crimson',pivot='mid',scale=40,zorder=5,minlength=0)
    
    #Add Discs
    AddDiscsToPlot(ax,pos)
    axis = [-0.02,0.02,-0.02,0.02]
    ax.axis(axis)
    ax.set_aspect('equal')
    # Turn off tick labels
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    fig.tight_layout()
    ax = set_size(6,6,ax)
    
    fig.savefig(plotName(cwd,Re,config,field,idx)+'.png')
    fig.clf()
    plt.close()
    return

# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    #return cwd+"/pd.txt"
    #cwd = cwd_PYTHON
    return cwd+"/pd.txt"

def GetPosData(cwd,time,config):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    if(config == 'V' or config == 'O'):
        pos = data[data['time'] == time*2.0]
    else:
        pos = data[data['time'] == time]
    #pos = data[data['time'] == time]
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
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 1024, 1024
    mxArr = np.array(mxList).reshape((Nx,Ny))
    myArr = np.array(myList).reshape((Nx,Ny))
    WArr  = np.array(WList).reshape((Nx,Ny))
    PArr  = np.array(PList).reshape((Nx,Ny))
    UxArr = np.array(UxList).reshape((Nx,Ny))
    UyArr = np.array(UyList).reshape((Nx,Ny))
    #return (mxArr, myArr, WArr, PArr, UxArr, UyArr, axArr, ayArr)
    return (mxArr, myArr, WArr, PArr, UxArr, UyArr)

# caculates the gradient of a grid of values, provided an
# isotropic grid spacing h. Throws away outer skin of points
def grad(mat,h=1):
    gy = 0.5*(np.roll(mat,-1,axis=1) - np.roll(mat,1,axis=1))/h
    gx = 0.5*(np.roll(mat,-1,axis=0) - np.roll(mat,1,axis=0))/h
    return (gx[1:-1,1:-1],gy[1:-1,1:-1])

# calculates the laplacian of a grid of values, provided
# an isotropic grid spacing h. Throws away outer skin of points
def lapl(mat,h=1):
    l = -4*mat + np.roll(mat,-1,axis=0) + np.roll(mat,1,axis=0) +\
        np.roll(mat,-1,axis=1) + np.roll(mat,1,axis=1)
    return l[1:-1,1:-1]/(h*h)

def distance_between(vec1,vec2):
    return np.hypot(vec1[0] - vec2[0], vec1[1] - vec2[1])

def RemoveCellData(cwd,mx,my,force_conv_x2,force_conv_y2,force_vari_x,force_vari_y,force_pres_x,force_pres_y,force_diff_x,force_diff_y):
    global PERIOD, DX, RADIUS_LARGE, RADIUS_SMALL
    #In this function, we will be setting the values of forces to zero if they are within the path of the spherobot
    #For each set of pos data, change value of field cell to zero if not already done so
    forces = [force_conv_x2,force_conv_y2,force_vari_x,force_vari_y,force_pres_x,force_pres_y,force_diff_x,force_diff_y]
    #Get Pos Data
    time_start, time_end = 10.0, 10.0+PERIOD
    data = pd.read_csv(pname(cwd),delimiter=' ')
    if(config == 'V' or config == 'O'):
        pos = data[data['time'] >= time_start*2.0].copy()
        pos = pos[pos['time'] < time_end*2.0]
    else:
        pos = data[data['time'] >= time_start].copy()
        pos = pos[pos['time'] < time_end]
    pos = pos.reset_index(drop=True)
    assert len(pos['time']) == 20

    NU = int(np.ceil(RADIUS_LARGE/DX))+1
    NL = int(np.ceil(RADIUS_SMALL/DX))+1
    for idx in range(len(pos['time'])):
        #Large Sphere 1
        aXU, aYU = pos.loc[idx,'aXU'], pos.loc[idx,'aYU']
        NaXU, NaYU = int((aXU + 0.05)/DX) - 1, int((aYU + 0.05)/DX) - 1
        centeraU = np.array([mx[NaXU,NaYU],my[NaXU,NaYU]])
        for jdx in range(NaXU - (NU+1),NaXU + (NU+1)):
            for kdx in range(NaYU - (NU+1), NaYU + (NU+1)):
                cell = np.array([mx[jdx,kdx],my[jdx,kdx]])
                dist = distance_between(centeraU,cell)
                for force in forces:
                    if dist <= NU*DX and force[jdx,kdx] != 0.0:
                        force[jdx,kdx] = 0.0
        #Large Sphere 2
        bXU, bYU = pos.loc[idx,'bXU'], pos.loc[idx,'bYU']
        NbXU, NbYU = int((bXU + 0.05)/DX) - 1, int((bYU + 0.05)/DX) - 1
        centerbU = np.array([mx[NbXU,NbYU],my[NbXU,NbYU]])
        for jdx in range(NbXU - (NU+1),NbXU + (NU+1)):
            for kdx in range(NbYU - (NU+1), NbYU + (NU+1)):
                cell = np.array([mx[jdx,kdx],my[jdx,kdx]])
                dist = distance_between(centerbU,cell)
                for force in forces:
                    if dist <= NU*DX and force[jdx,kdx] != 0.0:
                        force[jdx,kdx] = 0.0
        #Small Sphere 1
        aXL, aYL = pos.loc[idx,'aXL'], pos.loc[idx,'aYL']
        NaXL, NaYL = int((aXL + 0.05)/DX) - 1, int((aYL + 0.05)/DX) - 1
        centeraL = np.array([mx[NaXL,NaYL],my[NaXL,NaYL]])
        for jdx in range(NaXL - (NL+1),NaXL + (NL+1)):
            for kdx in range(NaYL - (NL+1), NaYL + (NL+1)):
                cell = np.array([mx[jdx,kdx],my[jdx,kdx]])
                dist = distance_between(centeraL,cell)
                for force in forces:
                    if dist <= NL*DX and force[jdx,kdx] != 0.0:
                        force[jdx,kdx] = 0.0
        #Small Sphere 2
        bXL, bYL = pos.loc[idx,'bXL'], pos.loc[idx,'bYL']
        NbXL, NbYL = int((bXL + 0.05)/DX) - 1, int((bYL + 0.05)/DX) - 1
        centerbL = np.array([mx[NbXL,NbYL],my[NbXL,NbYL]])
        for jdx in range(NbXL - (NL+1),NbXL + (NL+1)):
            for kdx in range(NbYL - (NL+1), NbYL + (NL+1)):
                cell = np.array([mx[jdx,kdx],my[jdx,kdx]])
                dist = distance_between(centerbL,cell)
                for force in forces:
                    if dist <= NL*DX and force[jdx,kdx] != 0.0:
                        force[jdx,kdx] = 0.0

    force_conv_x2 = forces[0]
    force_conv_y2 = forces[1]
    force_vari_x = forces[2]
    force_vari_y = forces[3]
    force_pres_x = forces[4]
    force_pres_y = forces[5]
    force_diff_x = forces[6]
    force_diff_y = forces[7]
    
    return (force_conv_x2,force_conv_y2,force_vari_x,force_vari_y,force_pres_x,force_pres_y,force_diff_x,force_diff_y)

if __name__ == '__main__':
    
    #READ ALL AVG FILES IN A SIMULATION DIRECTORY
    #EXTRACT AVERAGE FIELD DATA INTO NUMPY ARRAYS
    #PLOT AVERAGED FIELD DATA
    #Simulation Parameters
    #simList = ['HB','SF','L','V']
    cwd_Re = cwd_PYTHON+'../'+config+'/Re'+Re+'/'
    #Extract Position Data
    cwd_POS = cwd_PYTHON+'../PosData/'+config+'/Re'+Re+'/'
    #Calculate # Periods
    DUMP_INT = 20.0
    nTime = GetPosDataLength(cwd_POS)
    nPer = int(np.trunc(1.0*nTime/DUMP_INT))
    #nPer = 2
    #Paths to data and plots
    cwd_DATA = cwd_Re+'/VTK/AVG/'
    cwd_FORCE = cwd_PYTHON + '../ForceData/'+config+'/'
    countPer = 0
    #Create Text file with maximum magnitude values for net force
    
    for countPer in range(nPer):
        if(countPer == 100):
            #strDir = cwd_PYTHON+"../Figures/AVG/"+config+"/"+Re+"/avgU/"
            #AVGPlot = pathlib.Path(strDir+'avgU_%i.png'%countPer)
            #AVGPlot = pathlib.Path(cwd_DATA+'AVG_Acc_%04d.csv'%countPer)
            AVGPlot = pathlib.Path(cwd_DATA+'AVG_%04d.csv'%countPer)
            #if not AVGPlot.exists ():
            if AVGPlot.exists ():
                start = t.clock()
                #Get Avg Field Data
                mx,my,avgW,avgP,avgUx,avgUy = GetAvgFieldData(cwd_DATA,countPer)
                #Calculate gradient and laplacian terms
                gradavgUx_x, gradavgUx_y = grad(avgUx,h=DX)
                gradavgUy_x, gradavgUy_y = grad(avgUy,h=DX)
                lapUx = lapl(avgUx,DX)
                lapUy = lapl(avgUy,DX)
                gradPx, gradPy = grad(avgP,h=DX)
                #Calculate force density terms
                force_conv_x1 = RHO*(avgUx[1:-1,1:-1]*gradavgUx_x + avgUy[1:-1,1:-1]*gradavgUx_y)
                force_conv_y1 = RHO*(avgUx[1:-1,1:-1]*gradavgUy_x + avgUy[1:-1,1:-1]*gradavgUy_y)
                force_diff_x = -MU*lapUx
                force_diff_y = -MU*lapUy
                force_pres_x = gradPx
                force_pres_y = gradPy
                nSims = 20
                cwd_VTK = cwd_Re+'/VTK/'
                udotgradusum_x = np.zeros((1022,1022))
                udotgradusum_y = np.zeros((1022,1022))
                avgu_deriv_t_x = np.zeros((1024,1024))
                avgu_deriv_t_y = np.zeros((1024,1024))
                for idxSim in range(0,nSims):
                    dumpIdx = countPer*nSims + idxSim
                    allData = ReadAllVTK(cwd_VTK,dumpIdx)
                    #Extract important values from allData
                    time,mx,my,wArr,pArr,uxArr,uyArr = ReadVTK(allData)
                    #Calculate gradient of u, then find vdotgradv
                    gradUx_x, gradUx_y = grad(uxArr,h=DX)
                    gradUy_x, gradUy_y = grad(uyArr,h=DX)
                    uxdotgradu = uxArr[1:-1,1:-1]*gradUx_x + uyArr[1:-1,1:-1]*gradUx_y
                    uydotgradu = uxArr[1:-1,1:-1]*gradUy_x + uyArr[1:-1,1:-1]*gradUy_y
                    #Add to the sum
                    udotgradusum_x += uxdotgradu
                    udotgradusum_y += uydotgradu
                #Calculate avg(du/dt)
                dumpIdx = countPer*nSims
                allData = ReadAllVTK(cwd_VTK,dumpIdx)
                #Extract important values from allData
                time_init,scrap1,scrap2,scrap3,scrap4,ux_init,uy_init = ReadVTK(allData)
                dumpIdx = countPer*(nSims+1)
                allData = ReadAllVTK(cwd_VTK,dumpIdx)
                #Extract important values from allData
                time_fin,scrap1,scrap2,scrap3,scrap4,ux_fin,uy_fin = ReadVTK(allData)
                avgu_deriv_t_x = (ux_fin - ux_init)/(time_fin - time_init)
                avgu_deriv_t_y = (uy_fin - uy_init)/(time_fin - time_init)
                #Calculate Avg vdotgradv
                force_conv_x2 = RHO*udotgradusum_x/(1.0*nSims)
                force_conv_y2 = RHO*udotgradusum_y/(1.0*nSims)
                force_vari_x = RHO*avgu_deriv_t_x
                force_vari_y = RHO*avgu_deriv_t_y
                
                #End of Force Density Calculations
                #Extract Position and Time Data
                time = np.round(0.05 + countPer*PERIOD,2)
                #print('time = ',time)
                posData = GetPosData(cwd_POS,time,config)
                #print(posData)
                #REMOVE CELL DATA WHERE SPHERE HAS BEEN LOCATED
                force_conv_x2,force_conv_y2,force_vari_x,force_vari_y,force_pres_x,force_pres_y,force_diff_x,force_diff_y = RemoveCellData(cwd_POS,mx[1:-1,1:-1],my[1:-1,1:-1],force_conv_x2,force_conv_y2,force_vari_x,force_vari_y,force_pres_x,force_pres_y,force_diff_x,force_diff_y)

                force_stre_x = force_pres_x + force_diff_x
                force_stre_y = force_pres_y + force_diff_y
                force_iner_x = force_conv_x2 + force_vari_x[1:-1,1:-1]
                force_iner_y = force_conv_y2 + force_vari_y[1:-1,1:-1]
                force_net_x = force_iner_x + force_stre_x
                force_net_y = force_iner_y + force_stre_y
                #Store Data in csv
                avgDict = {'mx':mx[1:-1,1:-1].flatten(),'my':my[1:-1,1:-1].flatten(),
                'fvx':force_vari_x[1:-1,1:-1].flatten(),'fvy':force_vari_y[1:-1,1:-1].flatten(),
                    'fcx':force_conv_x2.flatten(),'fcy':force_conv_y2.flatten(),
                    'fpx':force_pres_x.flatten(),'fpy':force_pres_y.flatten(),
                    'fdx':force_diff_x.flatten(),'fdy':force_diff_y.flatten(),
                    'fix':force_iner_x.flatten(),'fiy':force_iner_y.flatten(),
                    'fsx':force_stre_x.flatten(),'fsy':force_stre_y.flatten(),
                    'fnx':force_net_x.flatten(),'fny':force_net_y.flatten()}
                avgFieldData = pd.DataFrame(data=avgDict)
                #Export Avg Field Data to .csv file
                pathlib.Path(cwd_FORCE).mkdir(parents=True, exist_ok=True)
                avgFieldData.to_csv(cwd_FORCE+'/Force_%s_Re%s_%04d.csv'%(config,Re,countPer),index=False,sep=' ',float_format='%.5e')
                '''
                #Plot Averaged Field Data
                #Plot Pressure Field Data
                pars = [Re,config,'P',countPer]
                #PlotPressure(cwd_PYTHON,time,mx,my,avgP,posData,pars)
                pars[2] = 'force_vari'
                PlotForceDensity(cwd_PYTHON,time,mx[1:-1,1:-1],my[1:-1,1:-1],force_vari_x[1:-1,1:-1],force_vari_y[1:-1,1:-1],posData,pars)
                pars[2] = 'force_conv_1'
                PlotForceDensity(cwd_PYTHON,time,mx[1:-1,1:-1],my[1:-1,1:-1],force_conv_x1,force_conv_y1,posData,pars)
                pars[2] = 'force_conv_2'
                PlotForceDensity(cwd_PYTHON,time,mx[1:-1,1:-1],my[1:-1,1:-1],force_conv_x2,force_conv_y2,posData,pars)
                pars[2] = 'force_diff'
                PlotForceDensity(cwd_PYTHON,time,mx[1:-1,1:-1],my[1:-1,1:-1],force_diff_x,force_diff_y,posData,pars)
                pars[2] = 'force_pres'
                PlotForceDensity(cwd_PYTHON,time,mx[1:-1,1:-1],my[1:-1,1:-1],force_pres_x,force_pres_y,posData,pars)
                pars[2] = 'force_stre'
                PlotForceDensity(cwd_PYTHON,time,mx[1:-1,1:-1],my[1:-1,1:-1],force_stre_x,force_stre_y,posData,pars)
                pars[2] = 'force_iner'
                PlotForceDensity(cwd_PYTHON,time,mx[1:-1,1:-1],my[1:-1,1:-1],force_iner_x,force_iner_y,posData,pars)
                pars[2] = 'force_net'
                PlotForceDensity(cwd_PYTHON,time,mx[1:-1,1:-1],my[1:-1,1:-1],force_net_x,force_net_y,posData,pars)
                '''

                stend = t.clock()
                diff = stend - start
                print('Time to run for 1 period = %.5fs'%diff)
                sys.stdout.flush()

