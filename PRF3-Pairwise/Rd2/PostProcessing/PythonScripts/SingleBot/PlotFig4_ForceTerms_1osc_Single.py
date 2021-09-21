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
Re = sys.argv[1]
config = 'single'
MU = RHO*OMEGA*AMPLITUDE_SMALL*RADIUS_SMALL/float(Re)

# constructs a filepath for the plot omages of Re=$Re, config=$config, and field=$field
def plotName(cwd,Re,config,field,idx):
    strDir = cwd+"../Figures/AVGForce/{0}/".format(config)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
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
    return

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
    ax.contourf(mx_rot,my_rot,np.log10(forcemag),cmap='YlOrRd',levels=levels,extend='both')
    
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
              color='black',pivot='mid',scale=40,zorder=5,minlength=0)
    
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

def GetAvgFieldData(cwd,config,Re,idx):
    #Columns
    fieldData = pd.read_csv(cwd+'Force_%s_Re%s_%04d.csv'%(config,Re,idx),delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    fvxList  = fieldData['fvx'].values.tolist()
    fvyList  = fieldData['fvy'].values.tolist()
    fcxList = fieldData['fcx'].values.tolist()
    fcyList = fieldData['fcy'].values.tolist()
    fpxList = fieldData['fpx'].values.tolist()
    fpyList = fieldData['fpy'].values.tolist()
    fdxList = fieldData['fdx'].values.tolist()
    fdyList = fieldData['fdy'].values.tolist()
    fixList = fieldData['fix'].values.tolist()
    fiyList = fieldData['fiy'].values.tolist()
    fsxList = fieldData['fsx'].values.tolist()
    fsyList = fieldData['fsy'].values.tolist()
    fnxList = fieldData['fnx'].values.tolist()
    fnyList = fieldData['fny'].values.tolist()
    
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 1022, 1022
    mxArr = np.array(mxList).reshape((Nx,Ny))
    myArr = np.array(myList).reshape((Nx,Ny))
    fvxArr  = np.array(fvxList).reshape((Nx,Ny))
    fvyArr  = np.array(fvyList).reshape((Nx,Ny))
    fcxArr = np.array(fcxList).reshape((Nx,Ny))
    fcyArr = np.array(fcyList).reshape((Nx,Ny))
    fpxArr = np.array(fpxList).reshape((Nx,Ny))
    fpyArr = np.array(fpyList).reshape((Nx,Ny))
    fdxArr = np.array(fdxList).reshape((Nx,Ny))
    fdyArr = np.array(fdyList).reshape((Nx,Ny))
    fixArr = np.array(fixList).reshape((Nx,Ny))
    fiyArr = np.array(fiyList).reshape((Nx,Ny))
    fsxArr = np.array(fsxList).reshape((Nx,Ny))
    fsyArr = np.array(fsyList).reshape((Nx,Ny))
    fnxArr = np.array(fnxList).reshape((Nx,Ny))
    fnyArr = np.array(fnyList).reshape((Nx,Ny))
    return (mxArr, myArr, fvxArr, fvyArr, fcxArr, fcyArr,
            fpxArr, fpyArr, fdxArr, fdyArr, fixArr, fiyArr,
            fsxArr, fsyArr, fnxArr, fnyArr)

if __name__ == '__main__':
    
    #READ ALL AVG FILES IN A SIMULATION DIRECTORY
    #EXTRACT AVERAGE FIELD DATA INTO NUMPY ARRAYS
    #PLOT AVERAGED FIELD DATA
    #Paths to data and plots
    cwd_POS = cwd_PYTHON+'../PosData/'
    cwd_FORCE = cwd_PYTHON + '../ForceData/'
    countPer = 100
    mx,my,f_vari_x,f_vari_y,f_conv_x,f_conv_y,f_pres_x,f_pres_y,f_diff_x,f_diff_y,f_iner_x,f_iner_y,f_stre_x,f_stre_y,f_net_x,f_net_y = GetAvgFieldData(cwd_FORCE,config,Re,countPer)
    
    #Extract Position and Time Data
    time = np.round(0.05 + countPer*PERIOD,2)
    #print('time = ',time)
    posData = GetPosData(cwd_POS,time,Re)
    
    #Plot Averaged Field Data
    pars = [Re,config,'force_vari',countPer]
    PlotForceDensity(cwd_PYTHON,time,mx,my,f_vari_x,f_vari_y,posData,pars)
    pars[2] = 'force_conv_2'
    PlotForceDensity(cwd_PYTHON,time,mx,my,f_conv_x,f_conv_y,posData,pars)
    pars[2] = 'force_diff'
    PlotForceDensity(cwd_PYTHON,time,mx,my,f_diff_x,f_diff_y,posData,pars)
    pars[2] = 'force_pres'
    PlotForceDensity(cwd_PYTHON,time,mx,my,f_pres_x,f_pres_y,posData,pars)
    pars[2] = 'force_stre'
    PlotForceDensity(cwd_PYTHON,time,mx,my,f_stre_x,f_stre_y,posData,pars)
    pars[2] = 'force_iner'
    PlotForceDensity(cwd_PYTHON,time,mx,my,f_iner_x,f_iner_y,posData,pars)
    pars[2] = 'force_net'
    PlotForceDensity(cwd_PYTHON,time,mx,my,f_net_x,f_net_y,posData,pars)
    

