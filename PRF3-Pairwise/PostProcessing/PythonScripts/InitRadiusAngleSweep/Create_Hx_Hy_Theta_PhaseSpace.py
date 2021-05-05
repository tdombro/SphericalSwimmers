#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 15:06:23 2020

@author: thomas
"""

#MODULES
import os,sys
import numpy as np
import pandas as pd
from mpl_toolkits import mplot3d
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

#CONSTANTS
cwd_PYTHON = os.getcwd()
PERIOD = 0.1
DT = 1.0e-2
RADIUSLARGE = 0.002
RADIUSSMALL = 0.001

#Lists
#RLength
R = ["3","5","6","7"]
SwimDirList = ["SSL", "LSL", "Stat"]

allData = []

def StoreData(strR,strTheta, strConfig, strRe):
    #global axAll
    #Reset position data every Re
    pdData = []
    #Load position data
    #Periodic
    pdData = pd.read_csv(cwd_PYTHON+'/../Periodic/'+strR+'/'+strTheta+'/'+strConfig+'/'+strRe+'/pd.txt',delimiter = ' ')
    #Save only every 10 rows (Every period)
    pdData = pdData.iloc[::10]
    #Reindex so index corresponds to period number
    pdData = pdData.reset_index(drop=True)
    #Create CM variables
    pdData['aXCM'] = 0.8*pdData.aXU + 0.2*pdData.aXL
    pdData['aYCM'] = 0.8*pdData.aYU + 0.2*pdData.aYL
    pdData['bXCM'] = 0.8*pdData.bXU + 0.2*pdData.bXL
    pdData['bYCM'] = 0.8*pdData.bYU + 0.2*pdData.bYL
    #Find separation distance and relative heading for LS
    pdData['Hx'], pdData['Hy'], pdData['Theta'] = FindDistAngleBW(pdData['aXU'],pdData['aYU'],
                                                                  pdData['bXU'],pdData['bYU'],
                                                                  pdData['aXL'],pdData['aYL'],
                                                                  pdData['bXL'],pdData['bYL'])
    #Calculate deltaH and deltaTheta
    #pdData['dH'] = CalcDelta(pdData['H'])
    #pdData['dTheta'] = CalcDelta(pdData['Theta'])
    
    '''#Plot distBW vs time and angleBW vs time
    PlotHTheta(pdData['H'],pdData['Theta'],pdData['dH'],pdData['dTheta'],
               pdData['time'],strR,strTheta,strConfig,strRe)'''
    
    #PlotAllHTheta(pdData['H'],pdData['Theta'],pdData['dH'],pdData['dTheta'],pdData['time'])
    
    return pdData
    
def FindDistAngleBW(aXU,aYU,bXU,bYU,aXL,aYL,bXL,bYL):
    #Find Distance b/w the 2 swimmers (Large Spheres)
    distXU = bXU - aXU
    distYU = bYU - aYU
    distBW = np.hypot(distXU,distYU)
    #Find relative heading Theta (angle formed by 2 swimmers)
    #1) Find normal orientation for swimmer A
    dist_aX = aXU - aXL
    dist_aY = aYU - aYL
    dist_a = np.hypot(dist_aX,dist_aY)
    norm_aX, norm_aY = dist_aX/dist_a, dist_aY/dist_a
    #Find normal orientation for swimmer B
    dist_bX = bXU - bXL
    dist_bY = bYU - bYL
    dist_b = np.hypot(dist_bX,dist_bY)
    norm_bX, norm_bY = dist_bX/dist_b, dist_bY/dist_b
    #2) Calculate Theta_a
    Theta_a = np.zeros(len(norm_aX))
    for idx in range(len(norm_aX)):
        if(norm_aY[idx] >= 0.0):
            Theta_a[idx] = np.arccos(norm_aX[idx])
        else:
            Theta_a[idx] = -1.0*np.arccos(norm_aX[idx])+2.0*np.pi
    #print('Theta_a = ',Theta_a*180.0/np.pi)
    #3) Rotate A and B ccw by 2pi - Theta_a
    Angle = 2.0*np.pi - Theta_a
    #print('Angle = ',Angle*180.0/np.pi)
    #Create rotation matrix
    rotationMatrix = np.zeros((len(Angle),2,2))
    rotationMatrix[:,0,0] = np.cos(Angle)
    rotationMatrix[:,0,1] = -1.0*np.sin(Angle)
    rotationMatrix[:,1,0] = np.sin(Angle)
    rotationMatrix[:,1,1] = np.cos(Angle)
    #print('rotationMatrix[-1] = ',rotationMatrix[10,:,:])
    #Create swimmer position arrays
    pos_a = np.zeros((len(norm_aX),2))
    pos_b = np.zeros((len(norm_bX),2))
    pos_a[:,0],pos_a[:,1] = norm_aX,norm_aY
    pos_b[:,0],pos_b[:,1] = norm_bX,norm_bY
    #print('pos_a = ',pos_a)
    #print('pos_b = ',pos_b)
    #Apply rotation operator
    rotpos_a, rotpos_b = np.zeros((len(norm_aX),2)), np.zeros((len(norm_bX),2))
    for idx in range(len(norm_aX)):
        #print('pos_a = ',pos_a[idx,:])
        #print('pos_b = ',pos_b[idx,:])
        #print('rotationMatrix = ',rotationMatrix[idx])
        rotpos_a[idx,:] = rotationMatrix[idx,:,:].dot(pos_a[idx,:])
        rotpos_b[idx,:] = rotationMatrix[idx,:,:].dot(pos_b[idx,:])
        #print('rotpos_a = ',rotpos_a[idx,:])
        #print('rotpos_b = ',rotpos_b[idx,:])
    #print('rotpos_a = ',rotpos_a)
    #print('rotpos_b = ',rotpos_b)
    #Calculate angleBW
    angleBW = np.zeros(len(norm_bY))
    for idx in range(len(norm_bY)):
        if(rotpos_b[idx,1] >= 0.0):
            angleBW[idx] = np.arccos(rotpos_a[idx,:].dot(rotpos_b[idx,:]))
        else:
            angleBW[idx] = -1.0*np.arccos(rotpos_a[idx,:].dot(rotpos_b[idx,:]))+2.0*np.pi
    #print('angleBW = ',angleBW*180/np.pi)
    
    return (distXU,distYU,angleBW)

def CalcDelta(data):
    #Calculate the change in either H or Theta for every period
    delta = data.copy()
    delta[0] = 0.0
    for idxPer in range(1,len(delta)):
        delta[idxPer] = data[idxPer] - data[idxPer-1]
        
    return delta

def PlotHTheta(H,Theta,dH,dTheta,time,strR,strTheta,strConfig,strRe):
    #Create Folder for Plots
    #Periodic
    pathlib.Path('../HThetaPlots/Periodic/'+strR+'/'+strConfig+'/').mkdir(parents=True, exist_ok=True)
    cwd_DIST = cwd_PYTHON + '/../HThetaPlots/Periodic/'+strR+'/'+strConfig+'/'
    #GENERATE FIGURE
    csfont = {'fontname':'Times New Roman'}
    fig = plt.figure(num=0,figsize=(8,8),dpi=250)
    ax = fig.add_subplot(111,projection='polar')
    ax.set_title('H-Theta Space',fontsize=16,**csfont)
    #Create Mesh of Delta H and Theta
    mTheta, mH = np.meshgrid(Theta,H/RADIUSLARGE)
    mdTheta, mdH = np.meshgrid(dTheta,dH/RADIUSLARGE)
    #Plot Chnages in H and Theta in quiver plot (Vec Field)
    for idx in range(0,len(Theta),20):
        ax.quiver(mTheta[idx,idx],mH[idx,idx],
                      mdH[idx,idx]*np.cos(mTheta[idx,idx])-mdTheta[idx,idx]*np.sin(mTheta[idx,idx]),
                      mdH[idx,idx]*np.sin(mTheta[idx,idx])+mdTheta[idx,idx]*np.cos(mTheta[idx,idx]),
                      color='k')#,pivot='mid',scale=50)
    ax.scatter(Theta,H/RADIUSLARGE,c=time,cmap='rainbow',s=9)
    ax.set_ylim(0.0,np.amax(H)/RADIUSLARGE) #[Theta,H/RADIUSLARGE],
    #ax.set_xlim(90.0*np.pi/180.0,270.0*np.pi/180.0)
    #SetAxesParameters(ax,-0.05,0.05,0.02,0.01,0.01)
    
    fig.tight_layout()
    #fig.savefig(cwd_DIST+'Dist_'+strR+'R_'+strTheta+'_'+strConfig+'_'+strRe+'.png')
    fig.savefig(cwd_DIST+'HTheta'+strTheta+'_'+strRe+'.png')
    fig.clf()
    #plt.close()
    return

def PlotAllHTheta(data,name):
    #global figAll,axAll

    #Figure with every pt for each Re
    csfont = {'fontname':'Times New Roman'}
    figAll = plt.figure(num=1,figsize=(8,8),dpi=250)
    #figAll, axAll = plt.subplots(nrows=1, ncols=1,num=0,figsize=(16,16),dpi=250)
    axAll = figAll.add_subplot(111,projection='3d')
    axAll.set_title('H-Theta Space: '+name,fontsize=20,**csfont)
    axAll.set_zlabel(r'$\Theta/\pi$',fontsize=18,**csfont)
    axAll.set_ylabel('Hy/R\tR=2.0mm',fontsize=18,**csfont)
    axAll.set_xlabel('Hx/R\tR=2.0mm',fontsize=18,**csfont)
    
    print('data length = ',len(data))
    for idx in range(len(data)):
        #print('idx = ',idx)
        tempData = data[idx].copy()
        Hx = tempData['Hx']/RADIUSLARGE
        Hy = tempData['Hy']/RADIUSLARGE
        Theta = tempData['Theta']/np.pi
        #dH = tempData['dH']/RADIUSLARGE
        #dTheta = tempData['dTheta']
        time = tempData['time']
        #GENERATE FIGURE
        #print('Theta = ',Theta[0])
        #axAll.scatter(Theta[0],H[0],c='k')
        axAll.scatter3D(Hx,Hy,Theta,c=time,cmap='rainbow',s=9)
        '''#Create Mesh of Delta H and Theta
        mTheta, mH = np.meshgrid(Theta,H)
        mdTheta, mdH = np.meshgrid(dTheta,dH)
        normdH, normdTheta = mdH/np.hypot(mdH,mdTheta), mdTheta/np.hypot(mdH,mdTheta)
        #Plot Chnages in H and Theta in quiver plot (Vec Field)
        
        for idx in range(0,len(Theta),10):
            #Polar
            axAll.quiver(mTheta[idx,idx],mH[idx,idx],
                          normdH[idx,idx]*np.cos(mTheta[idx,idx])-normdTheta[idx,idx]*np.sin(mTheta[idx,idx]),
                          normdH[idx,idx]*np.sin(mTheta[idx,idx])+normdTheta[idx,idx]*np.cos(mTheta[idx,idx]),
                          color='k',scale=25,headwidth=3,minshaft=2)
            
            
            axAll.quiver(mTheta[idx,idx],mH[idx,idx],
                          mdH[idx,idx]*np.cos(mTheta[idx,idx])-mdTheta[idx,idx]*np.sin(mTheta[idx,idx]),
                          mdH[idx,idx]*np.sin(mTheta[idx,idx])+mdTheta[idx,idx]*np.cos(mTheta[idx,idx]),
                          color='k',scale=25,headwidth=3,minshaft=2)
            
        
            #Cartesian
            axAll.quiver(mTheta[idx,idx],mH[idx,idx],
                          normdTheta[idx,idx],normdH[idx,idx],
                          color='k',scale=25,headwidth=3,minshaft=2)
        '''
        #print('Theta: ',len(Theta))
        #print(Theta)
        #axAll.set_ylim(0.0,np.amax(H)/RADIUSLARGE) #[Theta,H/RADIUSLARGE],
        #ax.set_xlim(90.0*np.pi/180.0,270.0*np.pi/180.0)
        #SetAxesParameters(ax,-0.05,0.05,0.02,0.01,0.01)
    axAll.set_xlim(0.0,8.0)
    axAll.set_ylim(0.0,8.0)
    axAll.set_zlim(0.0,2.0)
    figAll.tight_layout()
    figAll.savefig(cwd_PYTHON + '/../HThetaPlots/Periodic/HThetaAll_3D_'+name+'.png')
    figAll.clf()
    plt.close()
    print('About to exit PlotAll')

    return

if __name__ == '__main__':
    
    #dirsVisc = [d for d in os.listdir(cwd_STRUCT) if os.path.isdir(d)]
    #The main goal of this script is to create an H-Theta Phase Space of all
    #simulations for every period elapsed.
    #Example: 20s of sim time = 200 periods. 200 H-Theta Plots
    #We may find that we can combine the H-Theta data after steady state has been reached
    #1) For each simulation, store position data for every period (nothing in between)
    #2) Calculate separation distance between large spheres (H)
    #3) Calculate relative heading (Theta)
    #4) Calculate deltaH and deltaTheta (change after one Period)
    #5) Plot H vs Theta (Polar) for each period

    count = 0
    for idxR in range(len(R)):
        cwd_R = cwd_PYTHON + '/../Periodic/' + R[idxR]
        dirsTheta = [d for d in os.listdir(cwd_R) if not d.startswith('.')]
        for Theta in dirsTheta:
            cwd_THETA = cwd_R + '/' + Theta
            dirsConfig = [d for d in os.listdir(cwd_THETA) if not d.startswith('.')]
            for Config in dirsConfig:
                for idxRe in range(len(SwimDirList)):
                    if(idxRe == 2):
                        count += 1
    print('count = ',count)
    
    allData = [None]*(count)
    
    #For each Reynolds number (aka Swim Direction)
    for idxRe in range(len(SwimDirList)):
        count = 0
        #For each Initial Separation Distance
        #H(0) = Rval*R; R = radius of large sphere = 0.002m = 2mm
        for idxR in range(len(R)):
            cwd_R = cwd_PYTHON + '/../Periodic/' + R[idxR]
            #For each orientation at the specified length
            dirsTheta = [d for d in os.listdir(cwd_R) if not d.startswith('.')]
            for Theta in dirsTheta:
                cwd_THETA = cwd_R + '/' + Theta
                #For each configuration
                dirsConfig = [d for d in os.listdir(cwd_THETA) if not d.startswith('.')]
                for Config in dirsConfig:
                    print(R[idxR]+'\t'+Theta+'\t'+Config)
                    allData[count] = StoreData(R[idxR],Theta,Config,SwimDirList[idxRe])
                    count += 1
                    #PlotAllHTheta(axAll,data['H'],data['Theta'],data['dH'],data['dTheta'],data['time'])
                    #figAll.savefig(cwd_PYTHON + '/../HThetaPlots/Periodic/HThetaAll_SSL.png')
                    #sys.exit(0)
        #Plot All H-Theta space
        PlotAllHTheta(allData,SwimDirList[idxRe]) 
        #sys.exit(0)               