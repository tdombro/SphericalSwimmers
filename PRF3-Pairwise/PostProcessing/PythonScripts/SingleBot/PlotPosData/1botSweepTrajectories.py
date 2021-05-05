#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:54:13 2019

@author: thomas
"""

#MODULES
import os,sys
import numpy as np
import pandas as pd
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)
from scipy import stats
from scipy.signal import savgol_filter
import pathlib

mpl.rcParams['axes.linewidth'] = 1.5 #set the value globally

#Gather position data located in ../PositionData/Re$ReIdx/pd.txt
#Create a pandas database from each of them
#Plot Trajectories of small sphere and large sphere
#Look at data and decide how to handle calculating velocity

#CONSTANTS
cwd_PYTHON = os.getcwd()
PERIOD = 0.1
DT = 5.0e-3
RADIUSLARGE = 0.002
RADIUSSMALL = 0.001
FREQ = 10.0

csfont = {'fontname':'Times New Roman'}

#Lists
ReList = ([0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0,5.0,5.5,6.0,6.5,7.0,7.5,10.0,12.5,15.0,17.5,20.0,25.0,30.0,35.0,40.0,50.0,60.0])

# constructs a filepath for the pos data of Re = $Re
def pname(Re):
    return "/../PosData/pd_Re{0}.txt".format(Re)

def StoreData(cwd,Re):
    #Reset position data every Re
    data = []
    #Load position data
    #xL yL xS yS time
    posData = pd.read_csv(cwd+pname(Re),delimiter=' ')
    #Split up individual sphere data by index given
    topData = posData[posData['idx'] == 6].copy()
    botData = posData[posData['idx'] == 19].copy()
    topData = topData.sort_values(by=['time'])
    botData = botData.sort_values(by=['time'])
    topData = topData.reset_index(drop=True)
    botData = botData.reset_index(drop=True)
    #print(topData.head(100))
    #print(botData.head(100))
    #Rename columns to previous data names
    topData = topData.rename(columns={"x": "x1", "y": "y1"})
    botData = botData.rename(columns={"x": "x2", "y": "y2"})
    #print(topData.head(100))
    #print(botData.head(100))
    #Combine top and bot data to form previous dataframe used
    splitDict = {'xUp':topData['x1'],'yUp':topData['y1'],
        'xLow':botData['x2'],'yLow':botData['y2'],'time':topData['time']}
    data = pd.DataFrame(data=splitDict)
    #Add initial positioning
    initDict = {'xUp':[data.loc[0,'xUp']],'yUp':[0.0],'xLow':[data.loc[0,'xLow']],'yLow':[-5.0e-3],'time':[0.0]}
    initData = pd.DataFrame(data=initDict)
    data = data.append(initData)
    data = data.sort_values(by=['time'])
    data = data.reset_index(drop=True)
    #Center of Mass
    data["xCM"] = 0.8*data.xUp + 0.2*data.xLow
    data["yCM"] = 0.8*data.yUp + 0.2*data.yLow
    #Calculate Instantaneous Velocity
    data.loc[0,"vxCM"] = 0.0
    data.loc[0,"vyCM"] = 0.0
    for idx in range(1,len(data.time)):
        data.loc[idx,"vxCM"] = (data.loc[idx,'xCM'] - data.loc[idx-1,'xCM'])/DT
        data.loc[idx,"vyCM"] = (data.loc[idx,'yCM'] - data.loc[idx-1,'yCM'])/DT
    
    #Renormalize
    data['xCM'] /= RADIUSLARGE
    data['yCM'] /= RADIUSLARGE
    #Calculate vavg
    #Obtain period data
    lenData = len(data['time'])
    timestep=DT
    itPer = int(PERIOD/timestep)
    print('itPer = ',itPer)
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
    perData['xCM_shifted'] = perData['xCM'] - perData.loc[0,'xCM']
    #Calc vavg for each period
    perData['vavg_x'] = 0.0
    perData['vavg_y'] = 0.0
    perData['vavg'] = 0.0
    perData['disp'] = 0.0
    for idx in range(1,len(perData['time'])):
        perData.loc[idx,'vavg_x'] = (perData.loc[idx,'xCM_shifted'] - perData.loc[idx-1,'xCM_shifted'])/PERIOD
        perData.loc[idx,'vavg_y'] = (perData.loc[idx,'yCM_shifted'] - perData.loc[idx-1,'yCM_shifted'])/PERIOD
        if(perData.loc[idx,'vavg_y'] < 0.0):
            perData.loc[idx,'vavg'] = -1.0*np.hypot(perData.loc[idx,'vavg_x'],perData.loc[idx,'vavg_y'])
        else:
            perData.loc[idx,'vavg'] = np.hypot(perData.loc[idx,'vavg_x'],perData.loc[idx,'vavg_y'])
        if(perData.loc[idx,'yCM_shifted'] < 0.0):
            perData.loc[idx,'disp'] = -1.0*np.hypot(perData.loc[idx,'yCM_shifted'],perData.loc[idx,'xCM_shifted'])
        else:
            perData.loc[idx,'disp'] = np.hypot(perData.loc[idx,'yCM_shifted'],perData.loc[idx,'xCM_shifted'])
    #perData['vavg'] = np.hypot(perData['vavg_x'],perData['vavg_y'])
    #print(perData.head(10))
    #sys.exit(0)
    perData['vavg_x'] /= FREQ
    perData['vavg_y'] /= FREQ
    perData['vavg'] /= FREQ
    perData['nPer'] = perData['time']*FREQ
    #print(perData['vavg'])
    
    #Visualize Swimmer Position in Lab Frame
    #PlotPositions(cwd,perData,eps,sval,sweepName)

    #Remove any outliers (CM jump)
    perData['zscore'] = np.abs(stats.zscore(perData['vavg']))
    perData = perData[perData['zscore'] < 3]
    perData.reset_index(drop=True)
    
    #perData['vsmooth_x'] = SmoothCurve(perData,'vavg_x')
    perData['vsmooth_y'] = SmoothCurve(perData,'vavg_y')
    perData['vsmooth'] = SmoothCurve(perData,'vavg')
    
    PlotVvsTime(cwd,perData,Re)
    PlotDispvsTime(cwd,perData,Re)
        
    print('Final Position X = ',data.loc[len(data['time'])-1,'xCM'])
    print('Final Position Y = ',data.loc[len(data['time'])-1,'yCM'])

    #Calculate Swimmer steady state velocity
    if(float(Re) <= 10.0):
        #Swims SSL
        tempData = perData[perData['time'] >= 5.0].copy()
        avgData = tempData[tempData['time'] <= 18.0].copy()
        avgData = avgData.reset_index(drop=True)
        mean = np.amin(avgData['vsmooth'])
    elif(float(Re) > 15.0):
        #Swims LSL
        tempData = perData[perData['time'] >= 5.0].copy()
        avgData = tempData[tempData['time'] <= 18.0].copy()
        avgData = avgData.reset_index(drop=True)
        mean = np.amax(avgData['vsmooth'])
    elif(float(Re) == 12.5):
        tempData = perData[perData['time'] >= 5.0].copy()
        avgData = tempData[tempData['time'] <= 18.0].copy()
        avgData = avgData.reset_index(drop=True)
        #mean = avgData['vsmooth'].mean()
        mean = np.amin(avgData['vsmooth'])
    elif(float(Re) == 15.0):
        tempData = perData[perData['time'] >= 5.0].copy()
        avgData = tempData[tempData['time'] <= 15.0].copy()
        avgData = avgData.reset_index(drop=True)
        mean = avgData['vsmooth'].mean()
    '''
    #Calculate mean and std for avg vel
    if(Re == '50.0'):
        tempData = perData[perData['time'] >= 10.0].copy()
        avgData = tempData[tempData['time'] <= 15.0].copy()
    else:
        avgData = perData[perData['time'] >= 10.0].copy()
    mean = avgData['vavg_y'].mean()
    std  = avgData['vavg_y'].std()
    '''
    print('v_mean = ',mean)
    #print('v_std = ',std)

    #PlotTrajectory(data,Re)
    return mean
    #return (mean, std)

def SmoothCurve(data,varName):
    #After plotting the raw data for forces, the term mdu/dt was very noisy
    #We will apply a smoothing algorithm to all velocities in the y-direction
    
    #Use Savitsky-Golay filter from scipy
    yval = data[varName].tolist()
    varValue = savgol_filter(yval, 51, 3)
    
    return varValue

# constructs a filepath for the pos data of Re = $Re
def plotName(Par,Re):
    strDir = "../Figures/{0}/Rd2/".format(Par)
    pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    return strDir+"{0}_Re{1}.png".format(Par,Re)

def PlotVvsTime(cwd,data,Re):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(4,3),dpi=200)
    #ax[0].set_title(r'$\langle v_x \rangle$ vs time: $\epsilon$ = %.1f: $s^2$ = %.1f'%(eps,sval),**csfont,fontsize=16)
    #ax[0].set_xlabel(r'time (nPer)',**csfont,fontsize=12)
    #ax[0].set_ylabel(r'$\langle v_x \rangle$ (v/Rf)',**csfont,fontsize=12)
    ax.set_title(r'$\langle v \rangle$ vs time: Re = %s'%(Re),**csfont,fontsize=16)
    ax.set_xlabel(r'time (nPer)',**csfont,fontsize=12)
    ax.set_ylabel(r'$\langle v \rangle$ (v/rf)',**csfont,fontsize=12)
    #v_x
    #ax[0].plot(data['nPer'],data['vsmooth_x'],lw=1.5,c='r',zorder=5)
    #ax[0].scatter(data['nPer'],data['vavg_x'],s=4,c='k')
    #ax[0].plot([FREQ*(t_end-50),FREQ*t_end],[0.0,0.0],c='k')
    #v_y
    #ax.plot(data['nPer'],data['vsmooth_y'],lw=1.5,c='r',zorder=5)
    #ax.scatter(data['nPer'],data['vavg_y'],s=4,c='k')
    ax.plot(data['nPer'],data['vsmooth'],lw=1.5,c='r',zorder=5)
    ax.scatter(data['nPer'],data['vavg'],s=4,c='k')
    #ax[1].plot([FREQ*(t_end-50),FREQ*t_end],[0.0,0.0],c='k')
    #Axes
    ax.set_xlim(0.0,200.0)
    ax.set_ylim(-0.1,0.3)
    #ax[0].axis([FREQ*(t_end-50),FREQ*t_end,-0.01,0.01])
    #ax[1].axis([FREQ*(t_end-50),FREQ*t_end,-0.01,0.01])
    fig.tight_layout()
    fig.savefig(cwd+plotName('VvsTime',Re))
    fig.clf()
    plt.close()
    return

def PlotDispvsTime(cwd,data,Re):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(4,3),dpi=200)
    #ax[0].set_title(r'$\Delta x_{CM}$ vs time: $\epsilon$ = %.1f: $s^2$ = %.1f'%(eps,sval),**csfont,fontsize=16)
    #ax[0].set_xlabel(r'time (nPer)',**csfont,fontsize=12)
    #ax[0].set_ylabel(r'$\Delta x_{CM}$ (x/r)',**csfont,fontsize=12)
    #ax[0].plot(data['nPer'],data['xCM'],lw=1.5,c='k')
    #ax[0].set_xlim(FREQ*(t_end-50),FREQ*t_end)
    #ax[0].set_ylim(-2.0,2.0)
    ax.set_title(r'$\Delta y_{CM}$ vs time: Re = %s'%Re,**csfont,fontsize=16)
    ax.set_xlabel(r'time (nPer)',**csfont,fontsize=12)
    ax.set_ylabel(r'$\Delta y_{CM}$ (y/R)',**csfont,fontsize=12)
    #ax.plot(data['nPer'],data['yCM_shifted'],lw=1.5,c='k')
    ax.plot(data['nPer'],data['disp'],lw=1.5,c='k')
    #ax[1].set_xlim(FREQ*(t_end-50),FREQ*t_end)
    #ax[1].set_ylim(-2.0,2.0)
    ax.set_xlim(0.0,200.0)
    #ax.axis([FREQ*(t_end-50),FREQ*t_end,-0.01,0.01])
    fig.tight_layout()
    fig.savefig(cwd+plotName('MSD',Re))
    fig.clf()
    plt.close()
    return

def PlotTrajectory(data,Re):
    #Create Folder for Images
    pathlib.Path('../TrajectoryImages/').mkdir(parents=True, exist_ok=True)
    cwd_TRAJ = cwd_PYTHON + '/../TrajectoryImages/'
    #GENERATE FIGURE
    csfont = {'fontname':'Times New Roman'}
    fig = plt.figure(num=0,figsize=(12,4),dpi=120)
    ax = fig.add_subplot(131)
    ax.set_title('Trajectory: Re = %s'%Re,fontsize=20,**csfont)
    ax.set_ylabel('y (m)',fontsize=14,**csfont)
    ax.set_xlabel('x (m)',fontsize=14,**csfont)
    im = ax.scatter(data.xCM,data.yCM,c=data.vyCM,cmap='rainbow',s=1)
    ax.scatter(data.loc[len(data.xCM) - 1,"xCM"],data.loc[len(data.yCM) - 1,"yCM"],color='k',s=3,zorder=5)
    #ax.plot(data.xLow,data.yLow,color='b')
    ax.axis([-0.025,0.025,-0.025,0.025])
    SetAxesParameters(ax,-0.025,0.025,0.01,0.005,0.005)
    fig.colorbar(im, ax=ax)
    
    ax2=fig.add_subplot(132)
    ax2.set_title('Velocity: Re = %s'%Re,fontsize=20,**csfont)
    ax2.set_ylabel(r'$v_y$ (m)',fontsize=14,**csfont)
    ax2.set_xlabel('t (s)',fontsize=14,**csfont)
    ax2.plot(data.time,data.vyCM)
    ax2.axis([0.0,10.0,np.amin(data.vyCM) - 0.01,np.amax(data.vyCM)+0.01])
    
    ax3=fig.add_subplot(133)
    ax3.set_title('Velocity: Re = %s'%Re,fontsize=20,**csfont)
    ax3.set_ylabel(r'$v_x$ (m)',fontsize=14,**csfont)
    ax3.set_xlabel('t (s)',fontsize=14,**csfont)
    ax3.plot(data.time,data.vxCM)
    ax3.axis([0.0,10.0,np.amin(data.vxCM) - 0.01,np.amax(data.vxCM)+0.01])
    
    fig.tight_layout()
    fig.savefig(cwd_TRAJ+'Traj_Re%s.png'%Re)
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

def VvsReName():
    return "/../Figures/VvsRe_new.png"

def PlotVvsRe(cwd,velData):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6),dpi=200)
    ax.set_title(r'$\langle v \rangle$ vs Re',**csfont,fontsize=16)
    ax.set_xlabel('Re',**csfont,fontsize=12)
    ax.set_ylabel(r'$\langle v \rangle$ (v/Rf)',**csfont,fontsize=12)
    ax.plot([-10.0,60.0],[0.0,0.0],c='k')
    #ax.errorbar(velData['Re'],velData['v_avg'],yerr=2.0*velData['v_std'],c='k')
    ax.plot(velData['Re'],velData['v_avg'],c='k')
    ax.scatter(velData['Re'],velData['v_avg'],s=12,c='k',zorder=3)
    ax.set_xlim(0.0,60.0)
    ax.set_ylim(-0.1,0.35)
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    
    fig.tight_layout()
    fig.savefig(cwd+VvsReName())
    fig.clf()
    plt.close()
    return

if __name__ == '__main__':
    cwd_PYTHON = os.getcwd() + '/'
    velStats = np.zeros((len(ReList),2))
    mean = np.zeros(len(ReList))
    count = 0
    for Re in ReList:
        print('Re = %.1f'%Re)
        #velStats[count,0], velStats[count,1] = StoreData(cwd_PYTHON,str(Re))
        mean[count] = StoreData(cwd_PYTHON,str(Re))
        count += 1

    #velDict = {'Re':ReList,'v_avg':velStats[:,0],'v_std':velStats[:,1]}
    #velData = pd.DataFrame(data=velDict)
    velDict = {'Re':ReList,'v_avg':mean}
    velData = pd.DataFrame(data=velDict)
    PlotVvsRe(cwd_PYTHON,velData)
    velData.to_csv('SingleSwimmer_VelData_new.csv',index=False)



