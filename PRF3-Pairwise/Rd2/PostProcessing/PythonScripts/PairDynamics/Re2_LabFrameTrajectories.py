#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 21:59:46 2020

@author: thomas
"""

#MODULES
import os,sys
import re
import numpy as np
import pandas as pd
from mpl_toolkits import mplot3d
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)
from scipy.signal import savgol_filter
import pathlib
from matplotlib import animation
from IPython.display import display, Image, HTML
import subprocess
import time as tme

mpl.rcParams['axes.linewidth'] = 1.5 #set the value globally

#CONSTANTS
cwd_PYTHON = os.getcwd()
PERIOD = 0.1
DT = 5.0e-3
RADIUSLARGE = 0.002
RADIUSSMALL = 0.001
csfont = {'fontname':'Times New Roman'}

#Lists
#RLength
ThetaList=['0.0','22.5','45.0','67.5','90.0','112.5','135.0','157.5','180.0',
           '202.5','225.0','247.5','270.0','292.5','315.0','337.5']
HxList=['-0.5','0.5','1.5','2.5','3.5','4.5','5.5','6.5','7.5','8.5','9.5','10.5','11.5','12.5']
HyList=['-13','-11','-9','-7','-5','-3','-1','1','3','5','7','9']

allData = []

def PlotSwimmerTrajectories(data,trajData,ax):
    global PERIOD,RADIUSLARGE,RADIUSSMALL
    #Data Values
    axU      = data['aXU']
    axL      = data['aXL']
    ayU      = data['aYU']
    ayL      = data['aYL']
    bxU      = data['bXU']
    bxL      = data['bXL']
    byU      = data['bYU']
    byL      = data['bYL']
    
    #GENERATE FIGURE
    #Phase Space  2D (Hy vs Hx)
    #Swimmer 1
    ax.plot(trajData['aXCM'],trajData['aYCM'],color='g',zorder=5)
    ax.plot(trajData['bXCM'],trajData['bYCM'],color='b',zorder=5)
    ax.plot(trajData['aXL'],trajData['aYL'],color='gray',zorder=5)
    ax.plot(trajData['bXL'],trajData['bYL'],color='k',zorder=5)
    ax.scatter(data.loc[0,'aXCM'],data.loc[0,'aYCM'],c='tab:green',s=12,marker='s',zorder=6)
    ax.scatter(data.loc[0,'bXCM'],data.loc[0,'bYCM'],c='tab:green',s=12,marker='s',zorder=6)
    ax.scatter(data.loc[0,'parHx'],data.loc[0,'parHy'],c='tab:red',s=12,marker='s',zorder=6)
    
    #Add Swimmer locations
    Circle1 = plt.Circle((axU,ayU),1.0,color='k', clip_on=True)
    ax.add_artist(Circle1)
    Circle2 = plt.Circle((axL,ayL),0.5,color='k', clip_on=True)
    ax.add_artist(Circle2)
    Circle3 = plt.Circle((bxU,byU),1.0,color=(0.5,)*3, clip_on=True)
    ax.add_artist(Circle3)
    Circle4 = plt.Circle((bxL,byL),0.5,color=(0.5,)*3, clip_on=True)
    ax.add_artist(Circle4)
    
    #Plot boundary layer thickness around spheres
    #OMEGA = 2.0*np.pi/PERIOD
    AMP_SMALL = RADIUSLARGE*0.8*0.8
    delta = np.sqrt(AMP_SMALL*RADIUSSMALL/10.0)
    delta_norm = delta/RADIUSLARGE
    Circle5 = plt.Circle((axU,ayU),1.0+delta_norm, clip_on=True,fill=False,edgecolor='orange')
    ax.add_artist(Circle5)
    Circle6 = plt.Circle((axL,ayL),0.5+delta_norm, clip_on=True,fill=False,edgecolor='orange')
    ax.add_artist(Circle6)
    Circle7 = plt.Circle((bxU,byU),1.0+delta_norm, clip_on=True,fill=False,edgecolor='orange')
    ax.add_artist(Circle7)
    Circle8 = plt.Circle((bxL,byL),0.5+delta_norm, clip_on=True,fill=False,edgecolor='orange')
    ax.add_artist(Circle8)
    
    return ax

def FindMinValues(data):
    #Xmin
    aXUmin = np.amin(data['aXU'])-1.0
    aXLmin = np.amin(data['aXL'])-1.0
    bXUmin = np.amin(data['bXU'])-1.0
    bXLmin = np.amin(data['bXL'])-1.0
    xmin = min(-2.0,aXUmin,aXLmin,bXUmin,bXLmin)
    #Ymin
    aYUmin = np.amin(data['aYU'])-1.0
    aYLmin = np.amin(data['aYL'])-1.0
    bYUmin = np.amin(data['bYU'])-1.0
    bYLmin = np.amin(data['bYL'])-1.0
    ymin = min(-8.0,aYUmin,aYLmin,bYUmin,bYLmin)
    
    return xmin, ymin

def FindMaxValues(data):
    #Xmax
    aXUmax = np.amax(data['aXU'])+1.0
    aXLmax = np.amax(data['aXL'])+1.0
    bXUmax = np.amax(data['bXU'])+1.0
    bXLmax = np.amax(data['bXL'])+1.0
    xmax = max(7.0,aXUmax,aXLmax,bXUmax,bXLmax)
    #Ymax
    aYUmax = np.amax(data['aYU'])+1.0
    aYLmax = np.amax(data['aYL'])+1.0
    bYUmax = np.amax(data['bYU'])+1.0
    bYLmax = np.amax(data['bYL'])+1.0
    ymax = max(5.5,aYUmax,aYLmax,bYUmax,bYLmax)
    return xmax, ymax

if __name__ == '__main__':
    #Now that we have the simulation data, we can make trajectory plots
    #1: (In lab frame)
    #2: In lab frame centered on CM_pair at origin

    allData = pd.read_csv(cwd_PYTHON+'/LabTrajData_Re2.csv',delimiter=' ')

    allData['xCM'] = 0.5*(allData['aXCM'] + allData['bXCM'])
    allData['yCM'] = 0.5*(allData['aYCM'] + allData['bYCM'])

    lastData = allData[allData['boolStop'] == 1].copy()
    lastData = lastData.reset_index(drop=True)
    parHxList = lastData['parHx'].values.tolist()
    parHyList = lastData['parHy'].values.tolist()
    parThetaList = lastData['parThetaBW'].values.tolist()

    ThetaValue = float(sys.argv[1])
    for idx in range(len(parHxList)):
        Hx = parHxList[idx]
        Hy = parHyList[idx]
        Theta = parThetaList[idx]
        if Theta == ThetaValue:
            start = tme.clock()
            print('Hx = {0}: Hy = {1}: Theta = {2}'.format(Hx,Hy,Theta))
            sys.stdout.flush()
            thetaData = allData[allData['parThetaBW'] == Theta].copy()
            HxData = thetaData[thetaData['parHx'] == Hx].copy()
            HyData = HxData[HxData['parHy'] == Hy].copy()
            HyData = HyData.sort_values(by=['time'])
            HyData = HyData.reset_index(drop=True)
            print(len(HyData['time']))
            for idx in range(len(HyData['time'])):
                HyData.loc[idx,'time'] = np.round(HyData.loc[idx,'time'],1)
            #print(HyData['time'])
            
            strMovie = cwd_PYTHON+"/../Figures/LabTrajectories/T{0}/TrajMov_T{1}_Hx{2}_Hy{3}_.mp4".format(Theta,Theta,Hx,Hy)
            boolMovie = os.path.isfile(strMovie)
            print('boolMovie = ',boolMovie)
            if boolMovie == False:
                for idxTime in range(len(HyData['time'])):
                    time = np.round(PERIOD*idxTime,1)
                    data = HyData[HyData['time'] == time].copy()
                    data = data.reset_index(drop=True)
                    #print(data)
                    assert len(data['time']) == 1
                    #Get traj data
                    trajData = HyData[HyData['time'] <= time].copy()
                    trajData = trajData.reset_index(drop=True)
                    #2D PHASE SPACE PLOTS (Hy vs Hx) Configurations and Dynamics
                    #Figure and subplot for each Theta
                    fig, ax = plt.subplots(nrows=1,ncols=1,num=1,figsize=(6,6),dpi=250)
                    ax.set_title(r'time =%.1f: $\theta$ = %s: $H_y$ = %.1f vs. $H_x$ = %.2f'%(time,Theta,Hy,Hx),fontsize=15,**csfont)
                    ax.set_xlabel(r'$H_x$',fontsize=12,**csfont)
                    ax.set_ylabel(r'$H_y$',fontsize=12,**csfont)
                    xmin, ymin = FindMinValues(trajData)
                    xmax, ymax = FindMaxValues(trajData)
                    ax.axis([xmin,xmax,ymin,ymax])
                    ax.set_aspect('equal')
                    #Plot Swimmers and Traj Data
                    ax = PlotSwimmerTrajectories(data,trajData,ax)
                    #2D PLOTS
                    fig.tight_layout()
                    #plt.show()
                    #sys.exit(0)
                    strDir = cwd_PYTHON+"/../Figures/LabTrajectories/T{0}/".format(Theta)
                    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
                    fig.savefig(strDir+'Traj_T{0}_Hx{1}_Hy{2}_{3}_.png'.format(Theta,Hx,Hy,idxTime))
                    fig.clf()
                    plt.close()
                    if idxTime % 20 == 0:
                        print('time = ',time)
                os.chdir(cwd_PYTHON+"/../Figures/LabTrajectories/T{0}/".format(Theta))
                movieString = "ffmpeg -r 10 -i Traj_T{0}_Hx{1}_Hy{2}_%d_.png -vcodec libx264 -pix_fmt yuv420p -y TrajMov_T{0}_Hx{1}_Hy{2}_.mp4".format(Theta,Hx,Hy,Theta,Hx,Hy)
                os.system(movieString)
                fileName = "rm Traj_T{0}_Hx{1}_Hy{2}_*".format(Theta,Hx,Hy)
                os.system(fileName)
                os.chdir(cwd_PYTHON)
                stend = tme.clock()
                diff = stend - start
                print('Time to run for 1 period = %.5fs'%diff)
                #sys.exit(0)
