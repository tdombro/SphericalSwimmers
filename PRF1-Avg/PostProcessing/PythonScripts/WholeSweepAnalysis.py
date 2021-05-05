#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 12:01:38 2018

@author: thomas
"""

# Walks through the directories of a parameter sweep with each directory titled "ParameterA###/ParameterB###"
# For each simulation directory it:
# (1) extracts position data from the pd.txt files
# (2) Creates a dictionary with time: y1, y2 (large and small sphere positions)
# (3) Performs a linear regression to determine the velocity (via the slope of the regression)
#     Saves Maximum velocity during 10 oscillations
#     May be changed to record an average steady state velocity but unsure as of now
# (4) Plots the y data as a function of time, and the regression line (for visual confirmation)
# (5) Creates a master data.txt file that contains ParameterA, ParameterB, Velocity for the entire sweep

import collections
import os
import re
from shutil import copyfile

from time import perf_counter
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#FUTURE USE MAYBE?
def getRegressionBounds(time):
    #Default Regression bounds
    REGRESSION_START = 12 # start time of linear regression (s)
    REGRESSION_END = REGRESSION_START+REGRESSION_TIME
    REG_START_POSITION = int(REGRESSION_START/OUTPUT_STEP_SIZE)
    REG_END_POSITION = int(REGRESSION_END/OUTPUT_STEP_SIZE)
    #Record Last time value to make sure we 1)ran long enough for steady state
    LAST_TIME = time[-1]
        
    #Get Bounds of regression
    if(LAST_TIME < REGRESSION_END):
        if(LAST_TIME >= REGRESSION_TIME):
            START_TIME = LAST_TIME - REGRESSION_TIME
            REGRESSION_START = START_TIME
            REGRESSION_END = LAST_TIME
            REG_START_POSITION = int(REGRESSION_START/OUTPUT_STEP_SIZE)
            REG_END_POSITION = int(REGRESSION_END/OUTPUT_STEP_SIZE)
        else:
            REGRESSION_START = 0
            REGRESSION_END = LAST_TIME
            REG_START_POSITION = int(REGRESSION_START/OUTPUT_STEP_SIZE)
            REG_END_POSITION = int(REGRESSION_END/OUTPUT_STEP_SIZE)
            
    return (REGRESSION_START, REGRESSION_END, REG_START_POSITION, REG_END_POSITION)

#FUTURE USE MAYBE?
def linear_regression(time_to_xy):
    print('Starting linear_regression')
    start = perf_counter()
    
    just_times = list(time_to_xy.keys())
    just_y1s = [y1 for (x1, y1, x2, y2) in time_to_xy.values()]

    #Get Bounds of regression
    (REGRESSION_START,REGRESSION_END,REG_START_POSITION,REG_END_POSITION) = getRegressionBounds(just_times)

    regression_times=just_times[REG_START_POSITION:REG_END_POSITION]
    regression_y1s=just_y1s[REG_START_POSITION:REG_END_POSITION]
    
    #LINEAR REGRESSION
    slope, intercept, r_value, p_value, std_err = stats.linregress(regression_times,regression_y1s)

    end = perf_counter()
    print('Total linear regression time: {:.3f}s'.format(end - start))

    return(slope,intercept,r_value)

def MakeDirectory(directory,seed):
    if not os.path.exists(directory+'/'+str(seed)):
        os.makedirs(directory+'/'+str(seed))
    return

def main():
    start = perf_counter()
    
    #CONSTANTS
    PERIOD = 0.1
    FREQ = 10.0
    MAXAMP = 0.3
    RSL_DEFAULT = 0.75
    dt = 1.0e-4
    fignum=0

    #FILES
    file = open("WholeSweepData.txt","w")
    file.write('visc ampFrac rslFrac vel\n')
    
    #Get All RSL directories
    cwd_RSL = os.getcwd()
    dirsRSL = [d for d in os.listdir(cwd_RSL) if os.path.isdir(d)]
    #To test with a single directory (Uncomment)
    #dirsRSL = ['RSL1.0']
    #Loop over all RSL dir
    for dirRSL in dirsRSL:
        #Go into RSL#.#/Structures directory
        os.chdir(cwd_RSL+'/'+dirRSL)
        rsl = dirRSL
        rslString, rslValue = tuple(re.split('(\d.*)',rsl)[:2])
        os.chdir(cwd_RSL+'/'+dirRSL+'/Structures')
        #Save current working directory (Structures)
        cwd_STRUCT = os.getcwd()
        #Find all directories when using 'ls' in current directory (Viscosity)
        dirsVisc = [d for d in os.listdir(cwd_STRUCT) if os.path.isdir(d)]
        #print(dirsVisc)
        #to test with single directory (Uncomment)
        #dirsVisc = ['v0.04']
        for dirVisc in dirsVisc:
            #Go into 'v#.##' directory
            os.chdir(cwd_STRUCT+'/'+dirVisc)
            #Save current working directory (temporary)
            cwd_VISC = os.getcwd()
            #Find all directories when using 'ls' in current directory (resting spring length)
            dirsAmp = [d for d in os.listdir(cwd_VISC) if os.path.isdir(d)]
            #print(dirsAmp)
            #To test with single directory (Uncomment)
            #dirsAmp = ['amp0.8']
            for dirAmp in dirsAmp:
                dir_start = perf_counter()
                #Go into 'amp#.#' directory
                os.chdir(cwd_VISC+'/'+dirAmp)
                #Get Current Working Directory
                cwd_SIM = os.getcwd()
                #Obtain the amplitude and viscosity being used
                amp = os.path.basename(cwd_SIM)
                #Split 'rsl' into string 'rsl' and value 'value'
                ampString, ampValue = tuple(re.split('(\d.*)',amp)[:2])
                visc = os.path.split(os.path.dirname(cwd_SIM))[-1]
                #Split 'Re' dir into string 'Re' and value 'value'
                viscString, viscValue = tuple(re.split('(\d.*)',visc)[:2])
                #Calculate Reynolds Number of current simulation
                ReValue = 2.0*np.pi*FREQ*(float(ampValue)*MAXAMP)**2.0/float(viscValue)
                
                #Process pd.txt file. Store in pandas database
                #NEED TO ADD COLUMNS USING AWK BEFORE PROCESSING FILES (SEE add_header.sh)
                data = pd.read_csv('new.txt',delimiter = ' ')
                #print(data)
                #Save Values that occur when a period elapses
                savedValues = {}
                for i in range(0,len(data.time),int(PERIOD/dt)):
                    print('Full Period has elapsed. Time = %.3f'%data.time[i])
                    #Add to dictionary the COMy values and time for each elapsed period
                    if(float(ReValue) <= 10.0): 
                        if(data.time[i] <= 6.0): #Swimmer accelerates near wall and distorts data
                            savedValues[data.time[i]] = (data.upCOMy[i],data.lowCOMy[i])
                            print(savedValues[data.time[i]])
                    else:
                        savedValues[data.time[i]] = (data.upCOMy[i],data.lowCOMy[i])
                        print(savedValues[data.time[i]])
                
                #Perform Linear Regression for every 10 periods up till last timestep
                #ONLY IF max(time) > 1.0
                #Record largest velocity
                if(np.amax(data.time) > 1.0):
                    justTime = [y1 for (y1) in savedValues.keys()]
                    justUpCOMy = [y1 for (y1,y2) in savedValues.values()]
                    justLowCOMy = [y2 for (y1,y2) in savedValues.values()]
                    d = {'time':justTime,'upCOMy':justUpCOMy,'lowCOMy':justLowCOMy}
                    periodData = pd.DataFrame(data=d)
                    #print(periodData)
                    periodData = periodData.sort_values(by='time')
                    periodData = periodData.reset_index(drop=True)
                    print(periodData)
                    maxSlope = 0.0
                    j = 10
                    while(j < len(periodData.time)):
                        i = j-10
                        xData = periodData.time[i:j]
                        yData = periodData.upCOMy[i:j]
                        slope, intercept, r_value, p_value, std_err = stats.linregress(xData,yData)
                        print('Velocity = %.3f\ttmin = %.3f\ttmax = %.3f'%(slope,periodData.time[j-10],periodData.time[j]))
                        if(abs(maxSlope) < abs(slope) ):
                            maxSlope = slope
                            #maxSlope = max(maxSlope,abs(slope))
                        j+=1
                    '''for i in range(0,len(justTime)-10,1):
                        if(i+10 > len(justTime)):
                            j = -1
                        else:
                            j = i+10
                        #if(i+10 <= len(justTime)):
                        slope, intercept, r_value, p_value, std_err = stats.linregress(justTime[i:i+10],justUpCOMy[i:i+10])
                        print('Velocity = %.3f\ttmin = %.3f\ttmax = %.3f'%(slope,justTime[i],justTime[j]))
                        maxSlope = max(maxSlope,slope)'''
                    #ADD TO PLOT FUNCTION
                    #Plot Trajectory
                    #justTime = [y1 for (y1) in savedValues.keys()]
                    #justUpCOMy = [y1 for (y1,y2) in savedValues.values()]
                    #justLowCOMy = [y2 for (y1,y2) in savedValues.values()]
                    fig = plt.figure(num=fignum,figsize=(4,8),dpi=80)
                    ax = fig.add_subplot(211)
                    #Plot COM Trajectory
                    ax.set_title(visc+' '+amp+' '+rsl+': V = %.3f m/s'%maxSlope)
                    ax.plot(periodData.time,periodData.upCOMy)
                    ax.plot(periodData.time,periodData.lowCOMy)
                    #Plot Spring length difference
                    ax2 = fig.add_subplot(212)
                    ax2.set_title('Comparing Spring Lengths')
                    ax2.plot(data.time,data.currSL - float(rslValue)*RSL_DEFAULT,label='current')
                    ax2.plot(data.time,data.desSL - float(rslValue)*RSL_DEFAULT,label = 'desired')
                    ax2.legend(loc='upper left',fontsize='small')
                    ax2.axis([0.0,1.0,-0.3,0.3])
                    fig.savefig('Trajectory'+rsl+'.'+visc+'.'+amp+'.png')
                    fig.clf()
                    fignum+=1
                    #Write Parameters to data file
                    file.write('%.5e %.5e %.5e %.5e\n' %(float(viscValue),float(ampValue),float(rslValue),maxSlope))
                    dir_end = perf_counter()
                    print('Processed {} in {:.3f} seconds\n'.format(dirRSL+'Structures/'+dirVisc+'/'+dirAmp, dir_end - dir_start))
                    #Copy Image files to Image dir
                    #Include Subdirectories for organization
                    MakeDirectory(cwd_RSL,'Images/'+rsl+'/'+visc+'/'+amp)
                    cwd_IMAGE = cwd_RSL+'/Images/'+rsl+'/'+visc+'/'+amp+'/'
                    copyfile('Trajectory'+rsl+'.'+visc+'.'+amp+'.png', cwd_IMAGE+'Trajectory'+rsl+'.'+visc+'.'+amp+'.png')
                    #return
                else:
                    maxSlope = 'N/A'
                    file.write('%.5e %.5e %.5e %s\n' %(float(ReValue),float(ampValue),float(rslValue),maxSlope))
                    dir_end = perf_counter()
                    print('Processed {} in {:.3f} seconds\n'.format(dirRSL+'/Structures/'+dirVisc+'/'+dirAmp, dir_end - dir_start))
     
                #Processing and analysis for 1 amplitude complete. Return to viscosity directory (cwd_VISC)
                os.chdir(cwd_VISC)
            #Processing and analysis complete for all amplitudes and 1 viscosity. Return to cwd_STRUCT to change Viscosity
            os.chdir(cwd_STRUCT)
        #Processing and analysis complete for all viscosities and amplitudes. Return to cwd_RSL to change the resting spring length
        os.chdir(cwd_RSL)
    #All Processing and analysis complete! 
    #CLOSE FILES   
    file.close()
    end = perf_counter()
    print('Processed {:d} dirs in {:.3f} seconds'.format(len(dirsRSL), end - start))
    return

#--------------------__END MAIN__------------------------#
main()
