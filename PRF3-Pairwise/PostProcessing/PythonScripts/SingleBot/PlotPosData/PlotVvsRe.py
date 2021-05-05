#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 17:29:39 2020

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
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy import stats
import pathlib

#CONSTANTS
PERIOD = 0.1
FREQ = 10.0
DENS = 1000.0
RSMALL = 1.0e-3
RLARGE = 2.0*RSMALL
MAXAMP = 1.25*RSMALL
RSL_DEFAULT = RSMALL+RLARGE+1.5e-3
dt = 5.0e-3
FIGNUM=0
csfont = {'fontname':'Times New Roman'}

# constructs a filepath for the pos data of Re = $Re
def pname(Re):
    return "/../PosData/pd_Re{%.1f}.txt".format(Re)

def springName(Re,idx):
    strDir = "../Figures/SpringLength/Re{0}/".format(Re)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return "/../Figures/SpringLength/Re{0}/SpringLengthCheck_{1}.png".format(Re,idx)

# constructs a filepath for the pos data of Re = $Re
def plotName(Re):
    return "/../Figures/VvsTime/Re{:.1f}_VvsTime.png".format(Re)

def VvsReName():
    return "/../Figures/VvsRe.png"

def PlotSpringLength(cwd,data,Re,idx):
    global RSMALL
    RADIUS = RSMALL
    fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(8,3),dpi=200)
    ax[0].set_title(r'Re = %.1f: time = %.1f'%(Re,idx),**csfont,fontsize=16)
    ax[0].set_xlabel(r'time $\tau$',**csfont,fontsize=12)
    ax[0].set_ylabel(r'A/r',**csfont,fontsize=12)
    ax[1].set_xlabel(r'time $\tau$',**csfont,fontsize=12)
    ax[1].set_ylabel(r'% error',**csfont,fontsize=12)
    #Amplitude Comparison
    ax[0].plot(data['tau'],data['A_exp']/RADIUS,lw=1.5,c='k')
    ax[0].plot(data['tau'],data['A_sim']/RADIUS,lw=1.5,c='r')
    #Percent Error
    ax[1].plot(data['tau'],data['d_error'],lw=1.5,c='k')
    ax[0].axis([0.0,1.0,-0.0015/RADIUS,0.0015/RADIUS])
    minE, maxE = np.amin(data['d_error']), np.amax(data['d_error'])
    ax[1].axis([0.0,1.0,minE-0.02,maxE+0.02])
    fig.tight_layout()
    fig.savefig(cwd+springName(str(Re),str(idx)))
    fig.clf()
    plt.close()
    
    return

def StoreData(cwd,Re):
    global RSMALL, MAXAMP
    #Load position data
    #Columns
    #xL yL xS yS time
    data = pd.read_csv(cwd+pname(Re),delimiter=' ')
    #Add initial positioning
    #initDict = {'xS':[0.0],'yS':[-4.5e-3],'xL':[0.0],'yL':[0.0],'time':[0.0]}
    initDict = {'xS':[0.0],'yS':[-5.0e-3],'xL':[0.0],'yL':[0.0],'time':[0.0]}
    initData = pd.DataFrame(data=initDict)
    data = data.append(initData)
    data = data.sort_values(by=['time'])
    data = data.reset_index(drop=True)
    #print(data.head())
    #print(data.tail())
    #Calculate CM
    data['xCM'] = 0.8*data['xL']+0.2*data['xS']
    data['yCM'] = 0.8*data['yL']+0.2*data['yS']
    
    #Calculate Spring Length
    data['d_x']   = data['xL'] - data['xS']
    data['d_y']   = data['yL'] - data['yS']
    data['d_sim'] = np.hypot(data['d_x'],data['d_y'])
    AMPLITUDE = MAXAMP
    data['d_exp'] = RSL_DEFAULT + AMPLITUDE*np.sin(2.0*np.pi*FREQ*data['time'])
    data['A_sim'] = (data['d_sim'] - RSL_DEFAULT)
    data['A_exp'] = (data['d_exp'] - RSL_DEFAULT)
    #data['d_exp'] = RSL_DEFAULT + AMPLITUDE*np.sin(2.0*np.pi*FREQ*data['time'])
    data['d_error'] = 100.0*(data['d_sim'] - data['d_exp'])/data['d_exp']
    #Check spring length every 10s to ensure it's prescribed correctly
    maxTime = np.amax(data['time'])
    for idx in range(int(np.trunc(maxTime))):
        timeStart = float(idx)
        print(timeStart)
        oscData = data[data['time'] >= timeStart].copy()
        oscData = oscData[oscData['time'] < timeStart + PERIOD].copy()
        oscData = oscData.reset_index(drop=True)
        print('maxError = ',np.amax(oscData['d_error']))
        oscData['tau'] = (oscData['time'] - timeStart)*FREQ
        oscData['A_sim'] = (oscData['d_sim'] - RSL_DEFAULT)
        oscData['A_exp'] = (oscData['d_exp'] - RSL_DEFAULT)
        PlotSpringLength(cwd,oscData,Re,idx)       
    
    #Renormalize
    data['xL'] /= RSMALL
    data['yL'] /= RSMALL
    data['xS'] /= RSMALL
    data['yS'] /= RSMALL
    data['xCM'] /= RSMALL
    data['yCM'] /= RSMALL
    #Calculate vavg
    #Obtain period data
    lenData = len(data['time'])
    nPer = int(np.trunc(lenData/20.0))
    #print('nPer = ',nPer)
    perIdx = [15+20*idx for idx in range(nPer-1)]
    perData = data.iloc[perIdx]
    perData = perData.sort_values(by=['time','yCM','yL','yS'])
    perData = perData.reset_index(drop=True)
    #print(perData)
    #Shift yCM by init pos
    perData['yCM_shifted'] = perData['yCM'] - perData.loc[0,'yCM']
    #print(perData['yCM_shifted'])
    #Calc vavg for each period
    perData['vavg'] = 0.0
    for idx in range(1,len(perData['time'])):
        perData.loc[idx,'vavg'] = perData.loc[idx,'yCM_shifted'] - perData.loc[idx-1,'yCM_shifted']
    perData['vavg'] /= FREQ
    perData['tau'] = (perData['time'] - 0.075)*FREQ
    #print(perData['vavg'])
    print('Plot Re = '+str(Re))
    #Remove any outliers (CM jump)
    perData['zscore'] = np.abs(stats.zscore(perData['vavg']))
    if(Re == 40.0 or Re == 50.0):
        perData = perData[perData['zscore'] < 3]
        perData.reset_index(drop=True)
    perData['vsmooth'] = SmoothCurve(perData,'vavg')
    PlotVvsTime(cwd,perData,Re)
    #Extract Steady State Velocity
    #if Re < 40.0: avg last 3 periods
    #if Re >= 40.0, find maximum value
    if(Re < 40.0):
        vSS = np.mean(perData.loc[len(perData['vsmooth'])-3:,'vsmooth'])
    else:
        vSS = np.amax(perData['vsmooth'])
    return vSS

def PlotVvsTime(cwd,data,Re):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(4,3),dpi=200)
    ax.set_title(r'$v_{avg}$ vs time: Re = %.1f'%Re,**csfont,fontsize=16)
    ax.set_xlabel(r'time $\tau$',**csfont,fontsize=12)
    ax.set_ylabel(r'$v_{avg}$ (v/rf)',**csfont,fontsize=12)
    ax.plot(data['tau'],data['vsmooth'],lw=1.5,c='r')
    ax.scatter(data['tau'],data['vavg'],s=9,c='k')
    ax.set_ylim(-0.015,0.03)
    fig.tight_layout()
    fig.savefig(cwd+plotName(Re))
    fig.clf()
    plt.close()
    return

def SmoothCurve(data,varName):
    #After plotting the raw data for forces, the term mdu/dt was very noisy
    #We will apply a smoothing algorithm to all velocities in the y-direction
        
    #Use Savitsky-Golay filter from scipy
    yval = data[varName].tolist()
    varValue = savgol_filter(yval, 51, 3)
    
    return varValue

def PlotVvsRe(cwd,Re,vSS):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(4,3),dpi=200)
    ax.set_title(r'$v_{SS}$ vs Re',**csfont,fontsize=16)
    ax.set_xlabel('Re',**csfont,fontsize=12)
    ax.set_ylabel(r'$v_{SS}$ (v/rf)',**csfont,fontsize=12)
    ax.plot([-10.0,60.0],[0.0,0.0],c='k')
    ax.plot(Re,vSS,lw=1.5,c='r')
    ax.scatter(Re,vSS,s=9,c='k',zorder=3)
    ax.set_xlim(0.0,55.0)
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)

    fig.tight_layout()
    fig.savefig(cwd+VvsReName())
    fig.clf()
    plt.close()
    return

if __name__ == '__main__':
    #Current Directory
    cwd_PYTHON = os.getcwd()
    ReList = [0.5,1.0,2.0,3.0,4.0,5.0,10.0,12.0,
              14.0,16.0,18.0,20.0,30.0,40.0,50.0]
    SteadyState = np.zeros(len(ReList))
    idx = 0
    for Re in ReList:
        SteadyState[idx] = StoreData(cwd_PYTHON,Re)
        idx += 1
    PlotVvsRe(cwd_PYTHON,ReList,SteadyState)

