#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:54:13 2019

@author: thomas
"""

#MODULES
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)
from scipy.signal import savgol_filter
import pathlib

mpl.rcParams['axes.linewidth'] = 1.5 #set the value globally

#Gather position data located in ../PosData/$RLength/PI$Theta/$ANTIoPARA/$SSLoLSL/pd.txt
#Create a pandas database from each of them
#a,b = bot #; X,Y = dir; U,L = up/low spheres
#Create new variables: CM1, CM2, and distbw2CM
#Plot Trajectories of CM, color changes by time
#Plot distbw2CM vs time
#Look at data and decide how to handle calculating velocity

#CONSTANTS
cwd_PYTHON = os.getcwd()
PERIOD = 0.1
DT = 1.0e-2
RADIUSLARGE = 0.002
RADIUSSMALL = 0.001

#Lists
#RLength
R = ["3","5","6","7"]
#Periodic
#Theta = ["0", "1", "2", "3", "4"]
#Configurations
#ANTIoPARA = ["Anti", "Parallel","PerpL","PerpS"]
#Re
SSLoLSL = ["SSL", "LSL", "Stat"]


def StoreData(strR,strTheta, strConfig, strRe):
    #Reset position data every Re
    pdData = []
    #Load position data
    #Periodic
    pdData = pd.read_csv(cwd_PYTHON+'/../Periodic/PhiPI/'+strR+'/'+strTheta+'/'+strConfig+'/'+strRe+'/pd.txt',delimiter = ' ')
    #Create CM variables
    pdData["aXCM"] = 0.8*pdData.aXU + 0.2*pdData.aXL
    pdData["aYCM"] = 0.8*pdData.aYU + 0.2*pdData.aYL
    pdData["bXCM"] = 0.8*pdData.bXU + 0.2*pdData.bXL
    pdData["bYCM"] = 0.8*pdData.bYU + 0.2*pdData.bYL
    PlotTrajectory(pdData,strR,strTheta,strConfig,strRe)
    #Find distance between 2 swimmers LS
    pdData['distU'], pdData['angleU'] = FindDistBW(pdData['aXU'],pdData['aYU'],
                                                     pdData['bXU'],pdData['bYU'])
    #Find distance between 2 swimmers SS
    pdData['distL'], pdData['angleL'] = FindDistBW(pdData['aXL'],pdData['aYL'],
                                                     pdData['bXL'],pdData['bYL'])
    #Find distance between 2 swimmers CMs
    pdData['distCM'], pdData['angleCM'] = FindDistBW(pdData['aXCM'],pdData['aYCM'],
                                                     pdData['bXCM'],pdData['bYCM'])
    #Plot distBW vs time and angleBW vs time
    PlotDistAngle(pdData['distU'],pdData['angleU'],pdData['distL'],pdData['angleL'],
                  pdData['distCM'],pdData['angleCM'],pdData['time'],
                  strR,strTheta,strConfig,strRe)
    
    return

def PlotTrajectory(data,strR,strTheta,strConfig,strRe):
    #Create Folder for Images
    #Periodic
    pathlib.Path('../TrajectoryImages/Periodic/PhiPI/'+strR+'/'+strConfig+'/').mkdir(parents=True, exist_ok=True)
    cwd_TRAJ = cwd_PYTHON + '/../TrajectoryImages/Periodic/PhiPI/'+strR+'/'+strConfig+'/'
    #GENERATE FIGURE
    csfont = {'fontname':'Times New Roman'}
    fig = plt.figure(num=0,figsize=(4,4),dpi=120)
    ax = fig.add_subplot(111)
    ax.set_title('Trajectory: '+strR+'R: '+strTheta+': '+strConfig+': '+strRe,fontsize=16,**csfont)
    ax.set_ylabel('y (m)',fontsize=14,**csfont)
    ax.set_xlabel('x (m)',fontsize=14,**csfont)
    ax.scatter(data.aXCM,data.aYCM,c=data.time,cmap='viridis',s=1)
    ax.scatter(data.bXCM,data.bYCM,c=data.time,cmap='rainbow',s=1)
    #cbar = plt.colorbar(orientation="vertical")
    #ax.scatter(data.loc[len(data.xCM) - 1,"xCM"],data.loc[len(data.yCM) - 1,"yCM"],color='k',s=3,zorder=5)
    #ax.plot(data.xLow,data.yLow,color='b')
    #ax.axis([-0.05,0.05,-0.05,0.05])
    #SetAxesParameters(ax,-0.05,0.05,0.02,0.01,0.01)
    ax.axis('equal')
    
    fig.tight_layout()
    fig.savefig(cwd_TRAJ+'Traj_PI_'+strR+'R_'+strTheta+'_'+strConfig+'_'+strRe+'.png')
    fig.clf()
    plt.close()
    return
    
def SetAxesParameters(ax,xMin,xMax,xStep,xMinor,yMinor):
    #Axes Parameters
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    ax.set_axisbelow(False)
    ax.set_xticks(np.arange(xMin,xMax+xStep,step=xStep))
    ax.xaxis.set_minor_locator(MultipleLocator(xMinor))
    ax.yaxis.set_minor_locator(MultipleLocator(yMinor))
    return
    
def FindDistBW(aX,aY,bX,bY):
    #Find Distance b/w the 2 swimmers
    distX = bX - aX
    distY = bY - aY
    distBW = np.sqrt(distX*distX + distY*distY)
    #Find angle formed by 2 swimmers w/ respect to x-axis
    angleBW = np.arctan(distY/distX)
    for idx in range(len(angleBW)):
        if(distY[idx]/distX[idx] >= 0.0 and distX[idx] < 0.0):
            #Quadrant 3
            angleBW[idx] += np.pi
        elif(distY[idx]/distX[idx] <= 0.0 and distX[idx] < 0.0):
            #Quadrant 2
            angleBW[idx] += np.pi
        elif(distY[idx]/distX[idx] <= 0.0 and distX[idx] > 0.0):
            #Quadrant 4
            angleBW[idx] += 2.0*np.pi
    
    return (distBW, angleBW)

def PlotDistAngle(distU,angleU,distL,angleL,distCM,angleCM,time,strR,strTheta,strConfig,strRe):
    #Create Folder for Plots
    #Periodic
    pathlib.Path('../DistAnglePlots/Periodic/PhiPI/'+strR+'/'+strConfig+'/').mkdir(parents=True, exist_ok=True)
    cwd_DIST = cwd_PYTHON + '/../DistAnglePlots/Periodic/PhiPI/'+strR+'/'+strConfig+'/'
    #GENERATE FIGURE
    csfont = {'fontname':'Times New Roman'}
    fig = plt.figure(num=0,figsize=(8,8),dpi=120)
    ax1 = fig.add_subplot(221)
    ax1.set_title('DistBW: '+strR+'R: '+strTheta+': '+strConfig+': '+strRe,fontsize=16,**csfont)
    ax1.set_ylabel('distBW (R)',fontsize=14,**csfont)
    ax1.set_xlabel('time (s)',fontsize=14,**csfont)
    ax1.plot(time,distU/RADIUSLARGE,color='r',label='LS-LS',alpha=0.7)
    ax1.plot(time,distL/RADIUSLARGE,color='b',label='SS-SS',alpha=0.7)
    ax1.plot(time,distCM/RADIUSLARGE,color='k',label='CM-CM',alpha=0.7)
    ax1.legend(loc='best',fontsize='small')
    ax1.axis([0.0,np.amax(time)+0.1,0.0,np.amax([distU,distL,distCM])/RADIUSLARGE])
    ax2 = fig.add_subplot(222,projection='polar')
    ax2.set_title('Large Sphere',fontsize=16,**csfont)
    #ax3.set_yticklabels([])
    ax2.scatter(angleU,distU/RADIUSLARGE,c=time,cmap='rainbow',s=2)
    ax2.set_ylim(0.0,np.amax(distU)/RADIUSLARGE)
    ax3 = fig.add_subplot(223,projection='polar')
    ax3.set_title('Small Sphere',fontsize=16,**csfont)
    #ax3.set_yticklabels([])
    ax3.scatter(angleL,distL/RADIUSLARGE,c=time,cmap='rainbow',s=2)
    ax3.set_ylim(0.0,np.amax(distL)/RADIUSLARGE)
    ax4 = fig.add_subplot(224,projection='polar')
    ax4.set_title('Center of Mass',fontsize=16,**csfont)
    #ax3.set_yticklabels([])
    ax4.scatter(angleCM,distCM/RADIUSLARGE,c=time,cmap='rainbow',s=2)
    ax4.set_ylim(0.0,np.amax(distCM)/RADIUSLARGE)
    #ax3.axis([0.0,2.0,np.amin(distBW),np.amax(distBW)])
    #SetAxesParameters(ax,-0.05,0.05,0.02,0.01,0.01)
    
    fig.tight_layout()
    fig.savefig(cwd_DIST+'Dist_PI_'+strR+'R_'+strTheta+'_'+strConfig+'_'+strRe+'.png')
    fig.clf()
    plt.close()
    return

if __name__ == '__main__':
    
    #dirsVisc = [d for d in os.listdir(cwd_STRUCT) if os.path.isdir(d)]
    
    for idxR in range(len(R)):
        cwd_R = cwd_PYTHON + '/../Periodic/PhiPI/' + R[idxR]
        dirsTheta = [d for d in os.listdir(cwd_R) if not d.startswith('.')]
        for Theta in dirsTheta:
            cwd_THETA = cwd_R + '/' + Theta
            dirsConfig = [d for d in os.listdir(cwd_THETA) if not d.startswith('.')]
            for Config in dirsConfig:
                for idxRe in range(len(SSLoLSL)):
                    StoreData(R[idxR],Theta,Config,SSLoLSL[idxRe])