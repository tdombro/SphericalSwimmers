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
DENS = 1000.0
RADIUS_LARGE = 3.0e-1
RADIUS_SMALL = 1.5e-1
RSL_DEFAULT = 3.0*RADIUS_LARGE
LENGTH = 75.0
dt = 1.0e-4
FIGNUM=0
csfont = {'fontname':'Times New Roman'}

# constructs a filepath for the pos data of Re = $Re
def pname(eps,sval,Par):
    return "/../PosData/{0}/pd_e{1}_s{2}.txt".format(Par,eps,sval)


# constructs a filepath for the pos data of Re = $Re
def plotName(eps,sval,Par):
    strDir = "../Figures/{0}/VvsTime/".format(Par)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return "/../Figures/{0}/VvsTime/VvsTime_e{1}_s{2}.png".format(Par,eps,sval)

# constructs a filepath for the pos data of Re = $Re
def dispName(eps,sval,Par):
    strDir = "../Figures/{0}/MSD/".format(Par)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return "/../Figures/{0}/MSD/MSD_e{1}_s{2}.png".format(Par,eps,sval)


def oscName(eps,sval,Par):
    if(Par == 'epsilon'):
        strDir = "../Figures/{0}/TopDown/e{1}/".format(Par,eps)
    else:
        strDir = "../Figures/{0}/TopDown/s{1}/".format(Par,sval)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return strDir

def springName(eps,sval,Par,idx):
    if(Par == 'epsilon'):
        strDir = "../Figures/{0}/SpringLength/e{1}/".format(Par,eps)
    else:
        strDir = "../Figures/{0}/SpringLength/s{1}/".format(Par,sval)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return '/'+strDir+"SpringLengthCheck_e{0}_s{1}_{2}.png".format(eps,sval,idx)

def PlotSpringLength(cwd,data,eps,sval,par,idx):
    global RADIUS_LARGE
    AMPLITUDE = eps*RADIUS_LARGE
    fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(8,3),dpi=200)
    ax[0].set_title(r'%s Method: $\epsilon$ = %.1f: $s^2$ = %.1f: time = %.1f'%('CIB',eps,sval,PERIOD*idx),**csfont,fontsize=16)
    ax[0].set_xlabel(r'time $\tau$',**csfont,fontsize=12)
    ax[0].set_ylabel(r'A/R',**csfont,fontsize=12)
    ax[1].set_xlabel(r'time $\tau$',**csfont,fontsize=12)
    ax[1].set_ylabel(r'% error',**csfont,fontsize=12)
    #Amplitude Comparison
    ax[0].plot(data['tau'],data['A_exp']/RADIUS_LARGE,lw=1.5,c='k')
    ax[0].plot(data['tau'],data['A_sim']/RADIUS_LARGE,lw=1.5,c='r')
    #Percent Error
    ax[1].plot(data['tau'],data['d_error'],lw=1.5,c='k')
    ax[0].axis([0.0,1.0,-1.1*AMPLITUDE/RADIUS_LARGE,1.1*AMPLITUDE/RADIUS_LARGE])
    minE, maxE = np.amin(data['d_error']), np.amax(data['d_error'])
    ax[1].axis([0.0,1.0,minE - 0.02,maxE + 0.02])
    fig.tight_layout()
    fig.savefig(cwd+springName(str(eps),str(sval),par,str(idx)))
    fig.clf()
    plt.close()
    
    return

def StoreData(cwd,eps,sval,sweepName):
    global RADIUS_SMALL, RADIUS_LARGE, RSL_DEFAULT
    #Load position data
    #Columns
    #xL yL xS yS curr_spr des_spr time
    arrPos = np.loadtxt(cwd+pname(str(eps),str(sval),sweepName))
    posDict = {'x1':arrPos[:,0],'y1':arrPos[:,1],'x2':arrPos[:,2],'y2':arrPos[:,3],
               'd_sim':arrPos[:,4],'d_exp':arrPos[:,5],'time':arrPos[:,6]}
    data = pd.DataFrame(data=posDict)
    data = data.sort_values(by=['time'])
    data = data.reset_index(drop=True)
    #Add initial positioning
    initDict = {'x2':[0.0],'y2':[-4.5e-1],'x1':[0.0],'y1':[4.5e-1],
                'd_sim':[RSL_DEFAULT],'d_exp':[RSL_DEFAULT],'time':[0.0]}
    initData = pd.DataFrame(data=initDict)
    data = data.append(initData)
    data = data.sort_values(by=['time'])
    data = data.reset_index(drop=True)
    #print(data.head(100))
    #sys.exit(0)
    #Calculate CM
    data['xCM'] = (8.0/9.0)*data['x1'] + (1.0/9.0)*data['x2']
    data['yCM'] = (8.0/9.0)*data['y1'] + (1.0/9.0)*data['y2']
    #Calculate Spring Length
    data['A_sim'] = (data['d_sim'] - RSL_DEFAULT)
    data['A_exp'] = (data['d_exp'] - RSL_DEFAULT)
    #data['d_exp'] = RSL_DEFAULT + AMPLITUDE*np.sin(2.0*np.pi*FREQ*data['time'])
    data['d_error'] = 100.0*(data['d_sim'] - data['d_exp'])/data['d_exp']
    #Check spring length every 10s to ensure it's prescribed correctly
    maxTime = np.amax(data['time'])
    for idx in range(int(np.trunc(maxTime/PERIOD))):
        timeStart = PERIOD*idx
        print(timeStart)
        oscData = data[data['time'] >= timeStart].copy()
        oscData = oscData[oscData['time'] < timeStart + PERIOD].copy()
        oscData = oscData.reset_index(drop=True)
        print('maxError = ',np.amax(oscData['d_error']))
        oscData['tau'] = (oscData['time'] - timeStart)*FREQ
        oscData['A_sim'] = (oscData['d_sim'] - RSL_DEFAULT)
        oscData['A_exp'] = (oscData['d_exp'] - RSL_DEFAULT)
        PlotSpringLength(cwd,oscData,eps,sval,sweepName,idx) 
    #return
 
    #Renormalize
    data['x1'] /= RADIUS_LARGE
    data['y1'] /= RADIUS_LARGE
    data['x2'] /= RADIUS_LARGE
    data['y2'] /= RADIUS_LARGE
    data['xCM'] /= RADIUS_LARGE
    data['yCM'] /= RADIUS_LARGE
    #Calculate vavg
    #Obtain period data
    lenData = len(data['time'])
    timestep=1.0e-4
    if(eps == 1.0):
        if(sval == 10.0 or sval == 20.0 or sval == 30.0):
            timestep = 5.0e-5
    itPer = int(PERIOD/timestep)
    nPer = int(np.trunc(lenData/float(itPer)))
    if(nPer == 0):
        return []
    print('nPer = ',nPer)
    perIdx = [itPer*idx for idx in range(nPer)]
    perData = data.iloc[perIdx]
    perData = perData.sort_values(by=['time','yCM'])
    perData = perData.reset_index(drop=True)
    #Plot Average Spring Length
    #PlotAvgSpringLength(cwd,perData,Model,Nres,L)  
    #Shift yCM by init pos
    perData['yCM_shifted'] = perData['yCM'] - perData.loc[0,'yCM']
    #Calc vavg for each period
    perData['vavg_x'] = 0.0
    perData['vavg_y'] = 0.0
    perData['vavg'] = 0.0
    for idx in range(1,len(perData['time'])):
        perData.loc[idx,'vavg_x'] = (perData.loc[idx,'xCM'] - perData.loc[idx-1,'xCM'])/PERIOD
        perData.loc[idx,'vavg_y'] = (perData.loc[idx,'yCM_shifted'] - perData.loc[idx-1,'yCM_shifted'])/PERIOD
        #perData.loc[idx,'vavg'] = perData.loc[idx,'yCM_shifted'] - perData.loc[idx-1,'yCM_shifted']
    perData['vavg'] = np.hypot(perData['vavg_x'],perData['vavg_y'])
    perData['vavg_x'] /= FREQ
    perData['vavg_y'] /= FREQ
    perData['vavg'] /= FREQ
    perData['nPer'] = perData['time']*FREQ
    #print(perData['vavg'])
    
    #Visualize Swimmer Position in Lab Frame
    #PlotPositions(cwd,perData,eps,sval,sweepName)
    
    '''
    #Remove any outliers (CM jump)
    perData['zscore'] = np.abs(stats.zscore(perData['vavg']))
    perData = perData[perData['zscore'] < 3]
    perData.reset_index(drop=True)
    '''
    
    #perData['vsmooth_x'] = SmoothCurve(perData,'vavg_x')
    #perData['vsmooth_y'] = SmoothCurve(perData,'vavg_y')
    
    PlotVvsTime(cwd,perData,eps,sval,sweepName)
    PlotDispvsTime(cwd,perData,eps,sval,sweepName)
        
    print('Final Position X = ',data.loc[len(data['time'])-1,'xCM'])
    print('Final Position Y = ',data.loc[len(data['time'])-1,'yCM'])

    return perData
    

def PlotVvsTime(cwd,data,eps,sval,par):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(4,3),dpi=200)
    #ax[0].set_title(r'$\langle v_x \rangle$ vs time: $\epsilon$ = %.1f: $s^2$ = %.1f'%(eps,sval),**csfont,fontsize=16)
    #ax[0].set_xlabel(r'time (nPer)',**csfont,fontsize=12)
    #ax[0].set_ylabel(r'$\langle v_x \rangle$ (v/Rf)',**csfont,fontsize=12)
    ax.set_title(r'$\langle v_y \rangle$ vs time: $\epsilon$ = %.1f: $s^2$ = %.1f'%(eps,sval),**csfont,fontsize=16)
    ax.set_xlabel(r'time (nPer)',**csfont,fontsize=12)
    ax.set_ylabel(r'$\langle v_y \rangle$ (v/Rf)',**csfont,fontsize=12)
    ax.scatter(data['nPer'],data['vavg_y'],s=4,c='k')
    fig.tight_layout()
    fig.savefig(cwd+plotName(str(eps),str(sval),par))
    fig.clf()
    plt.close()
    return

def PlotDispvsTime(cwd,data,eps,sval,par):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(4,3),dpi=200)
    ax.set_title(r'$\Delta y_{CM}$ vs time',**csfont,fontsize=16)
    ax.set_xlabel(r'time (nPer)',**csfont,fontsize=12)
    ax.set_ylabel(r'$\Delta y_{CM}$ (y/R)',**csfont,fontsize=12)
    ax.plot(data['nPer'],data['yCM_shifted'],lw=1.5,c='k')
    fig.tight_layout()
    fig.savefig(cwd+dispName(str(eps),str(sval),par))
    fig.clf()
    plt.close()

def SmoothCurve(data,varName):
    #After plotting the raw data for forces, the term mdu/dt was very noisy
    #We will apply a smoothing algorithm to all velocities in the y-direction
        
    #Use Savitsky-Golay filter from scipy
    yval = data[varName].tolist()
    varValue = savgol_filter(yval, 51, 3)
    
    return varValue

def PlotPositions(cwd,data,eps,sval,par):
    #Get Variables
    time = data['time']
    x1 = data['x1']
    y1 = data['y1']
    x1_i, y1_i = x1.iloc[0],y1.iloc[0]
    x2 = data['x2']
    y2 = data['y2']
    x2_i, y2_i = x2.iloc[0],y2.iloc[0]
    xCM = data['xCM']
    yCM = data['yCM']
    xCM_i, yCM_i = xCM.iloc[0],yCM.iloc[0]
    for idx in range(len(time)):
        #Create Figure
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(4,4),dpi=200)
        ax.set_title(r'$\tau$ = %.1f'%(float(time[idx]/PERIOD)),**csfont,fontsize=16)
        #Plot Top Sphere
        Circle1 = plt.Circle((x1.iloc[idx],y1.iloc[idx]),2.0,color='k', clip_on=False)
        ax.add_artist(Circle1)
        ax.plot([-0.5*LENGTH,0.5*LENGTH],[y1_i,y1_i],'b')
        ax.plot([x1.iloc[idx]],[y1.iloc[idx]],'b',marker='o',markersize=5)
        #print('x1 = ',x1.iloc[idx])
        #print('y1 = ',y1.iloc[idx])
        #Plot Bottom Sphere
        Circle2 = plt.Circle((x2.iloc[idx],y2.iloc[idx]),1.0,color='k', clip_on=False)
        ax.add_artist(Circle2)
        ax.plot([-0.5*LENGTH,0.5*LENGTH],[y2_i,y2_i],'b')
        ax.plot([x2.iloc[idx]],[y2.iloc[idx]],'b',marker='o',markersize=5)
        #Plot COM using 'COM'
        ax.plot([-0.5*LENGTH,0.5*LENGTH],[xCM_i,yCM_i],'r')
        ax.plot([xCM.iloc[idx]],[yCM.iloc[idx]],'r',marker='o',markersize=5)
        #Set axes boundaries
        ax.axis([-8.0,8.0,-8.0,8.0])
        ax.set_aspect(aspect='equal')
        fig.tight_layout()
        plotName = oscName(str(eps),str(sval),par)
        fig.savefig(cwd+'/'+plotName+'TD_%04d.png'%idx)
        fig.clf()
        #sys.exit(0)
        plt.close()

                             
    return

if __name__ == '__main__':
    #Current Directory
    cwd_PYTHON = os.getcwd()
    Nres = '512'
    epsList = [0.1,0.2,0.4,1.0]
    
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(4,4),dpi=200)
    colorList = ['purple','blue','green','red']
    count = 0
    for eps in epsList:
        if(eps == 0.1):
            svalList = [1.0,2.0,3.0,4.0,5.0,7.5,10.0,50.0,100.0,200.0,300.0,400.0,500.0]
        elif(eps == 0.2):
            svalList = [1.0,2.0,3.0,4.0,5.0,7.5,10.0,50.0,100.0,150.0,300.0]
        elif(eps == 0.4):
            svalList = [0.5,1.0,2.0,3.0,4.0,5.0,7.5,10.0,50.0,100.0,150.0,200.0,300.0]
        elif(eps == 1.0):
            svalList = [1.0,2.0,3.0,4.0,5.0,10.0,20.0,30.0]
    
        #s^2 sweep
    
        felDict = {'avgU':[0.0],'eps2':[eps*eps],'s2':[0.0],'eps':[eps]}
        felData = pd.DataFrame(data=felDict)
    
        fig1, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(4,4),dpi=200)
        ax1.set_title(r'v vs. time: $s^2$ sweep')
        ax1.set_xlabel(r'time (nPer)',**csfont,fontsize=12)
        ax1.set_ylabel(r'$\langle v_y \rangle$ (v/Rf)',**csfont,fontsize=12)
        cm = matplotlib.cm.viridis
        norm.autoscale([0.0,max(svalList)])
        sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
        sm.set_array([])
    
        for sval in svalList:
            print('sval = ',sval)
            data = StoreData(cwd_PYTHON,eps,sval,'eps'+str(eps))
            if(len(data) != 0):
                data['eps'], data['sval'] = eps, sval
                ax1.plot(data['nPer'],data['vavg_y'],lw=2,c=cm(norm(sval)),label=sval)
                simDict = {'avgU':[data.loc[len(data)-1,'vavg_y']*RADIUS_LARGE*FREQ],'eps2':[eps*eps],'s2':[sval],'eps':[eps]}
                simData = pd.DataFrame(data=simDict)
                felData = pd.concat([felData,simData],ignore_index=True)
        ax1.legend(loc='best',fontsize='x-small')
        fig1.tight_layout()
        fig1.savefig(cwd_PYTHON+'/../Figures/eps'+str(eps)+'/VvsTime.png')
    
        #Plot Felderhof s^2 plot
        felData = felData.sort_values(by=['s2'])
        felData = felData.reset_index(drop=True)
        ax.set_title(r'$\langle U_{sw} \rangle$ vs. $s^2$')
        ax.set_xlabel(r'$\epsilon s^2$',**csfont,fontsize=12)
        ax.set_ylabel(r'$10^3\langle U_{sw} \rangle$/($\omega$R$\epsilon^2$)',**csfont,fontsize=12)
        #ax.plot(felData['s2']*felData['eps'],1.0e3*felData['avgU']/(OMEGA*RADIUS_LARGE*felData['eps2']),c=colorList[count],label='eps'+str(eps))
        #ax.scatter(felData['s2']*felData['eps'],1.0e3*felData['avgU']/(OMEGA*RADIUS_LARGE*felData['eps2']),c=colorList[count])
        ax.plot(felData['s2'],1.0e3*felData['avgU']/(OMEGA*RADIUS_LARGE*felData['eps2']),c=colorList[count],label='eps'+str(eps))
        ax.scatter(felData['s2'],1.0e3*felData['avgU']/(OMEGA*RADIUS_LARGE*felData['eps2']),c=colorList[count])
        #xval = np.linspace(0.0,30.0,100)
        #yval = 1000.0*(0.00234 - 6.0e-5*xval**2)/(OMEGA*RADIUS_LARGE)
        #ax.plot(xval,yval,c='r')
        #ax.axis([0.0,30.0,0.0,2.5])
        count += 1
    ax.plot([0.0,500.0],[0.0,0.0],c='k')
    ax.set_xlim(0.0,500.0)
    ax.set_ylim(-5.0,5.0)
    ax.legend(loc='upper left',fontsize='x-small')
    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'/../Figures/fel_epsSweep_v_vs_s2.png')
    ax.set_xlim(0.0,50.0)
    fig.tight_layout()
    fig.savefig(cwd_PYTHON+'/../Figures/fel_epsSweep_v_vs_s2_zoom_noeps.png')
    fig.clf()
    plt.close()
    
    #Save Dataframe into .csv file
    #felData.to_csv(cwd_PYTHON+'/../Figures/eps'+str(eps)+'/fel_sSweep.csv',sep=' ',index=False)
