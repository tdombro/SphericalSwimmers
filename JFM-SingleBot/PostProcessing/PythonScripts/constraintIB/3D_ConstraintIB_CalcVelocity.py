#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 13:48:11 2020

@author: thomas
"""

#Calculate Steady State Velocity for all Re
#Store pos data
#Store pos data for each period
#Calc vavg for each period
#Plot vavg vs time for each Re
#Extract where it plateaus
#Plot vavg vs. Re

#Visualize Oscillations and CoM Placement throughout

#MODULES
import os, sys
import numpy as np
import pandas as pd
import matplotlib
#mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy import stats
import pathlib
from matplotlib.colors import Normalize
norm = Normalize()

#CONSTANTS
PERIOD = 0.1
FREQ = 10.0
OMEGA = 2.0*np.pi*FREQ
DENS = 2.0
RADIUS_LARGE = 3.0e-1
RADIUS_SMALL = 1.5e-1
RSL_DEFAULT = 3.0*RADIUS_LARGE
EPSILON = 0.2
AMPLITUDE = EPSILON*RADIUS_LARGE
LENGTH = 8.0/RADIUS_LARGE
dt = 1.0e-4
FIGNUM=0
csfont = {'fontname':'Times New Roman'}

# constructs a filepath for the pos data of Re = $Re
def pname(cwd,eps,sval,Nres):
    return cwd+"/../PosData/eps{0}/s{1}/N{2}/pdCM.txt".format(eps,sval,Nres)

def getSpringLength(data):
    global RSL_DEFAULT, AMPLITUDE, OMEGA
    data['d(t)'] = RSL_DEFAULT + AMPLITUDE*np.sin(OMEGA*data['time'])
    return data

def StoreData(cwd,eps,sval,Nres,axLast,axAvgVel):
    global RADIUS_LARGE,OMEGA
    posData = pd.read_csv(pname(cwd,eps,sval,Nres),delimiter="\t",names=["time","xCM","yCM","zCM"],dtype={'time':np.float64,'xCM': np.float64,'yCM':np.float64,'zCM':np.float64})
    initDict = {'time':[0.0],'xCM':[posData.loc[0,'xCM']],'yCM':[posData.loc[0,'yCM']],'zCM':[posData.loc[0,'zCM']]} #[-0.45*8.0/9.0 + 0.45*1.0/9.0]
    initData = pd.DataFrame(data=initDict)
    posData = posData.append(initData)
    posData = posData.sort_values(by=['time'])
    posData = posData.reset_index(drop=True)
    posData = getSpringLength(posData)
    #RENORMALIZE
    posData['xCM'] /= RADIUS_LARGE
    posData['yCM'] /= RADIUS_LARGE
    posData['zCM'] /= RADIUS_LARGE
    posData['d(t)'] /= RADIUS_LARGE
    #Plot the instantaneous position for the last oscillation
    lastOsc = posData[posData['time'] >= 0.9]
    lastOsc = lastOsc.reset_index(drop=True)
    axLast = PlotLastOscPosition(axLast,lastOsc['xCM'],lastOsc['time'],Nres)
    #Plot the average velocity for each oscillation (10 total)
    #Get data every oscillation (0,0.1,0.2,etc.)
    lenData = len(posData['time'])
    timestep=1.0e-4
    '''
    if(eps == 1.0):
        if(sval == 10.0 or sval == 20.0 or sval == 30.0):
            timestep = 5.0e-5
    '''
    itPer = int(PERIOD/timestep)
    nPer = int(np.trunc(lenData/float(itPer)))
    if(nPer == 0):
        return []
    print('nPer = ',nPer)
    perIdx = [itPer*idx for idx in range(nPer)]
    perData = posData.iloc[perIdx]
    perData = perData.sort_values(by=['time'])
    perData = perData.reset_index(drop=True)
    #Calculate Average Velocity
    #Shift xCM by init pos
    perData['xCM_shifted'] = perData['xCM'] - perData.loc[0,'xCM']
    perData['yCM_shifted'] = perData['yCM'] - perData.loc[0,'yCM']
    perData['zCM_shifted'] = perData['zCM'] - perData.loc[0,'zCM']
    #Calc vavg for each period
    perData['avgvx'] = 0.0
    perData['avgvy'] = 0.0
    perData['avgvz'] = 0.0
    perData['avgvCM'] = 0.0
    for idx in range(1,len(perData['time'])):
        perData.loc[idx,'avgvx'] = (perData.loc[idx,'xCM_shifted'] - perData.loc[idx-1,'xCM_shifted'])/PERIOD
        perData.loc[idx,'avgvy'] = (perData.loc[idx,'yCM_shifted'] - perData.loc[idx-1,'yCM_shifted'])/PERIOD
        perData.loc[idx,'avgvz'] = (perData.loc[idx,'zCM_shifted'] - perData.loc[idx-1,'zCM_shifted'])/PERIOD
    '''
    perData['avgvCM'] = np.sqrt(perData['avgvx']*perData['avgvx']+
                              perData['avgvy']*perData['avgvy']+
                              perData['avgvz']*perData['avgvz'])
    '''
    perData['avgvCM'] = perData['avgvx'].copy()
    perData['avgvCM'] /= OMEGA
    perData['nPer'] = perData['time']*FREQ
    axAvgVel = PlotAverageVelocity(axAvgVel,perData['avgvCM'],perData['nPer'],Nres)
    print('Final Position X = ',posData.loc[len(posData['time'])-1,'xCM'])

    return (axLast,axAvgVel)

def PlotLastOscPosition(ax,xCM,time,Nres):
    if(Nres == 256):
        ax.plot(time,xCM,color='tab:blue',label='N256')
        ax.scatter(time[::100],xCM[::100],c='tab:blue',label='')
    elif(Nres == 512):
        ax.plot(time[::100],xCM[::100],color='tab:red',label='N512')
        ax.scatter(time[::100],xCM[::100],c='tab:red',label='')
    elif(Nres == 768):
        ax.plot(time[::100],xCM[::100],color='tab:green',label='N768')
        ax.scatter(time[::100],xCM[::100],c='tab:green',label='')
    return ax

def PlotAverageVelocity(ax,vCM,nPer,Nres):
    if(Nres == 256):
        ax.plot(nPer,vCM*1.0e3,color='tab:blue',label='N256')
        ax.scatter(nPer,vCM*1.0e3,c='tab:blue',label='')
    elif(Nres == 512):
        ax.plot(nPer,vCM*1.0e3,color='tab:red',label='N512')
        ax.scatter(nPer,vCM*1.0e3,c='tab:red',label='')
    elif(Nres == 768):
        ax.plot(nPer,vCM*1.0e3,color='tab:green',label='N768')
        ax.scatter(nPer,vCM*1.0e3,c='tab:green',label='')
    return ax

if __name__ == '__main__':
    #Current Directory
    cwd_PYTHON = os.getcwd()
    epsList = [0.2,0.01]
    NresList = [256,512,768]
    svalList = [10,200]
    for eps in epsList:
        if(eps == 0.2):
            NresList = [256,512,768]
            svalList = [10,200]
        elif(eps == 0.01):
            NresList = [512,768]
            svalList = [200,4000]
        for sval in svalList:
            figLast, axLast = plt.subplots(nrows=1,ncols=1,figsize=(6,6),dpi=200)
            axLast.set_title(r'Last Oscillation xCM Comparison: $s^2$={0}'.format(sval),**csfont,fontsize=16)
            axLast.set_xlabel('time (s)',**csfont,fontsize=12)
            axLast.set_ylabel(r'$x_{CM}$ (R)',**csfont,fontsize=12)
            figAvgVel, axAvgVel = plt.subplots(nrows=1,ncols=1,figsize=(6,6),dpi=200)
            axAvgVel.set_title(r'Avg Velocity vCM Comparison: $s^2$={0}'.format(sval),**csfont,fontsize=16)
            axAvgVel.set_xlabel('nPer',**csfont,fontsize=12)
            axAvgVel.set_ylabel(r'$v_{CM}$ ($10^3$v/R$\omega$)',**csfont,fontsize=12)
            for Nres in NresList:
                axLast,axAvgVel = StoreData(cwd_PYTHON,eps,sval,Nres,axLast,axAvgVel)
            #Format Last Oscillation Plot
            axLast.plot([0.8,1.1],[0.0,0.0],c='k')
            axLast.set_xlim(0.89,1.01)
            #axLast.set_ylim(-5.0,5.0)
            axLast.legend(loc='best',fontsize='x-small')
            figLast.tight_layout()
            strDir = cwd_PYTHON+"/../Figures/eps{0}/LastOscillation".format(eps)
            pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
            figLast.savefig(strDir+'/lastOsc_s2_{0}.png'.format(sval))
            figLast.clf()
            #Format Average Velocity Plot
            axAvgVel.plot([-0.01,20.01],[0.0,0.0],c='k')
            axAvgVel.set_xlim(0.0,20.0)
            #axLast.set_ylim(-5.0,5.0)
            axAvgVel.legend(loc='best',fontsize='x-small')
            figAvgVel.tight_layout()
            strDir = cwd_PYTHON+"/../Figures/eps{0}/AvgVel".format(eps)
            pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
            figAvgVel.savefig(strDir+'/AvgVel_s2_{0}.png'.format(sval))
            figAvgVel.clf()
        plt.close()


