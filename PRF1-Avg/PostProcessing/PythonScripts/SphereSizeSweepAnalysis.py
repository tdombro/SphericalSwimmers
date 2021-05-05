#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 11:32:38 2018

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
    radiusLarge = 0.3
    DENS = 2.0

    #FILES
    file = open("NewSphereSizeData.txt","w")
    file.write('Rlarge Rsmall Nu ampFrac rsl amp vel\n')
    
    #Get All RSL directories
    cwd_STRUCT = os.getcwd()
    #dirsRlarge = [d for d in os.listdir(cwd_STRUCT) if os.path.isdir(d)]
    dirsRlarge = ['R0.2','R0.225','R0.375','R0.75']
    #To test with a single directory (Uncomment)
    #dirsRlarge = ['v2.0']
    #Loop over all Visc dir
    for dirRlarge in dirsRlarge:
        #Go into RSL#.#/Structures directory
        os.chdir(cwd_STRUCT+'/'+dirRlarge)
        Rlarge = dirRlarge
        RlargeString, RlargeValue = tuple(re.split('(\d.*)',Rlarge)[:2])        
        #Save current working directory (Viscosity)
        cwd_Rlarge = os.getcwd()
        #Find all directories when using 'ls' in current directory (Viscosity)
        dirsRsmall = [d for d in os.listdir(cwd_Rlarge) if os.path.isdir(d)]
        #print(dirsRsmall)
        #to test with single directory (Uncomment)
        #dirsRsmall = ['amp0.4']
        for dirRsmall in dirsRsmall:
            #Go into 'amp#.#' directory
            os.chdir(cwd_Rlarge+'/'+dirRsmall)
            #Save current working directory (temporary)
            cwd_Rsmall = os.getcwd()
            Rsmall = dirRsmall
            RsmallString, RsmallValue = tuple(re.split('(\d.*)',Rsmall)[:2])
            #Find all directories when using 'ls' in current directory (resting spring length)
            dirsVisc = [d for d in os.listdir(cwd_Rsmall) if os.path.isdir(d)]
            #print(dirsVisc)
            #To test with single directory (Uncomment)
            #dirsVisc = ['0.02']
            for dirVisc in dirsVisc:
                dir_start = perf_counter()
                #Go into 'amp#.#' directory
                os.chdir(cwd_Rsmall+'/'+dirVisc)
                #Get Current Working Directory
                cwd_Visc = os.getcwd()
                #Obtain the amplitude and viscosity being used
                visc = dirVisc
                #Split 'rsl' into string 'rsl' and value 'value'
                viscString, viscValue = tuple(re.split('(\d.*)',visc)[:2])
                
                #Find all directories when using 'ls' in current directory (resting spring length)
                dirsAmp = [d for d in os.listdir(cwd_Visc) if os.path.isdir(d)]
                #print(dirsAmp)
                #To test with single directory (Uncomment)
                #dirsAmp = ['0.7']
                for dirAmp in dirsAmp:
                    dir_start = perf_counter()
                    #Go into 'amp#.#' directory
                    os.chdir(cwd_Visc+'/'+dirAmp)
                    #Get Current Working Directory
                    cwd_SIM = os.getcwd()
                    #Obtain the amplitude and viscosity being used
                    Amp = dirAmp
                    #Split 'rsl' into string 'rsl' and value 'value'
                    AmpString, AmpValue = tuple(re.split('(\d.*)',Amp)[:2])
                
                    #Calculate Reynolds Number of current simulation
                    ReValue = 2.0*np.pi*FREQ*(float(AmpValue)*MAXAMP)*float(RlargeValue)/(float(viscValue)/DENS)
                    print('Re = %.3e'%ReValue)
                    
                    #Process pd.txt file. Store in pandas database
                    #NEED TO ADD COLUMNS USING AWK BEFORE PROCESSING FILES (SEE add_header.sh)
                    data = pd.read_csv('new.txt',delimiter = ' ')
                    #print(data)
                    #Save Values that occur when a period elapses
                    savedValues = {}
                    
                    #Conditional for R=0.75 atm. Timestep change 
                    if(float(RlargeValue) == 0.75):
                        savedIndices = {}
                        k = 1
                        for j in range(0,len(data.time)):
                            if(data.lowCOMy[j] >= -3.0 and data.upCOMy[j] <= 3.0):
                                if(data.time[j] >= PERIOD*k):
                                    print('2. Full Period has elapsed. Time = %.3f'%data.time[j])
                                    #Add to dictionary the COMy values and time and CurrSLfor each elapsed period
                                    savedValues[data.time[j]] = (data.upCOMy[j],data.lowCOMy[j],data.currSL[j])
                                    savedIndices[k] = j
                                    print(savedValues[data.time[j]])
                                    k += 1

                        #Print N/A as the velocity atm
                        #maxSlope = 'N/A'
                        #file.write('%.5e %.5e %.5e %.5e %.5e %.5e %s\n' %(float(RlargeValue),float(RsmallValue),float(viscValue)/DENS,float(AmpValue),1.3*RSL_DEFAULT,MAXAMP*0.7,maxSlope))
                        #dir_end = perf_counter()
                        #print('Processed {} in {:.3f} seconds\n'.format(dirRlarge+'/'+dirRsmall+'/'+dirVisc+'/'+dirAmp, dir_end - dir_start))
     
                    else: #All other Rlarge values
                        #Try something related to position
                        for j in range(int(PERIOD/dt),len(data.time),int(PERIOD/dt)):
                            if(data.lowCOMy[j] >= -3.0 and data.upCOMy[j] <= 3.0):
                                print('2. Full Period has elapsed. Time = %.3f'%data.time[j])
                                #Add to dictionary the COMy values and time and CurrSLfor each elapsed period
                                savedValues[data.time[j]] = (data.upCOMy[j],data.lowCOMy[j],data.currSL[j])
                                print(savedValues[data.time[j]])
                        
                    #Perform Linear Regression for every 10 periods up till last timestep
                    #ONLY IF max(time) > 1.0
                    #Record largest velocity
                    if(np.amax(data.time) > 1.0):
                        justTime = [y1 for (y1) in savedValues.keys()]
                        justUpCOMy = [y1 for (y1,y2,y3) in savedValues.values()]
                        justLowCOMy = [y2 for (y1,y2,y3) in savedValues.values()]
                        justCurrSL = [y3 for (y1,y2,y3) in savedValues.values()] 
                        d = {'time':justTime,'upCOMy':justUpCOMy,'lowCOMy':justLowCOMy,'rsl':justCurrSL}
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
                                #Calculate Spring Length
                                springLength = np.mean(periodData.rsl[j-10:j])
                                print('rsl = %.3e'%springLength)
                                savedJ = j
                            j+=1

                        #Conditional for R=0.75: Timestep change
                        if(float(RlargeValue) == 0.75):
                            print(savedIndices)
                            #We know time range PERIOD[j-10,j]
                            #Translate this to a time. j = k previous. 
                            #So savedIndices[k=j] = index!
                            minIndex = savedIndices[savedJ-10]
                            maxIndex = savedIndices[savedJ]
                            
                        
                        else: #All other Rlarge values
                            #Time to calculate the actual amplitude
                            minIndex = int((savedJ-10)*PERIOD/dt)
                            maxIndex = int(savedJ*PERIOD/dt)
                            
                        maxSL = np.amax(data.currSL[minIndex:maxIndex])
                        minSL = np.amin(data.currSL[minIndex:maxIndex])
                        Amplitude = (maxSL - minSL)*0.5
                            
                        #ADD TO PLOT FUNCTION
                        #Plot Trajectory
                        fig = plt.figure(num=fignum,figsize=(8,24),dpi=120)
                        ax = fig.add_subplot(311)
                        #Plot COM Trajectory
                        ax.set_title(Rlarge+' '+Rsmall+' '+visc+' '+Amp+' : V = %.3f m/s'%maxSlope)
                        ax.plot(periodData.time,periodData.upCOMy)
                        ax.plot(periodData.time,periodData.lowCOMy)
                        #Plot Spring length difference
                        ax2 = fig.add_subplot(312)
                        ax2.set_title('Comparing Spring Lengths')
                        ax2.plot(data.time,data.currSL,label='current')
                        ax2.plot(data.time,data.desSL,label = 'desired')
                        ax2.legend(loc='upper left',fontsize='small')
                        ax2.axis([0.0,1.0,0.7,1.3]) #RSLFrac = 1.3
                        #Plot Spring Length Percent Difference
                        percentDifference = 100.0*(data.currSL - data.desSL)/data.desSL
                        ax3 = fig.add_subplot(313)
                        ax3.set_title('Comparing Spring Length: Percent Difference')
                        ax3.plot(data.time,percentDifference)
                        fig.savefig('Trajectory'+Rlarge+'.'+Rsmall+'.'+visc+'.'+Amp+'.png')
                        fig.clf()
                        fignum+=1
                        #Write Parameters to data file
                        file.write('%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n'
                                   %(float(RlargeValue),float(RsmallValue),float(viscValue)/DENS,float(AmpValue),springLength, Amplitude,maxSlope))
                        dir_end = perf_counter()
                        print('Processed {} in {:.3f} seconds\n'.format(dirRlarge+'/'+dirRsmall+'/'+dirVisc+'/'+dirAmp, dir_end - dir_start))
                        #Copy Image files to Image dir
                        #Include Subdirectories for organization
                        MakeDirectory(cwd_STRUCT,'AllImages/'+Rlarge+'/'+Rsmall+'/'+visc+'/'+Amp)
                        cwd_IMAGE = cwd_STRUCT+'/AllImages/'+Rlarge+'/'+Rsmall+'/'+visc+'/'+Amp+'/'
                        copyfile('Trajectory'+Rlarge+'.'+Rsmall+'.'+visc+'.'+Amp+'.png', cwd_IMAGE+'Trajectory'+Rlarge+'.'+Rsmall+'.'+visc+'.'+Amp+'.png')
                        #return
                    else:
                        maxSlope = 'N/A'
                        file.write('%.5e %.5e %.5e %.5e %.5e %.5e %s\n' %(float(RlargeValue),float(RsmallValue),float(viscValue)/DENS,float(AmpValue),RSL_DEFAULT*1.3,MAXAMP*0.7,maxSlope))
                        dir_end = perf_counter()
                        print('Processed {} in {:.3f} seconds\n'.format(dirRlarge+'/'+dirRsmall+'/'+dirVisc+'/'+dirAmp, dir_end - dir_start))
     
                    #Processing and analysis for 1 Amp complete. Return to Visc directory (cwd_VISC)
                    os.chdir(cwd_Visc)
                #Processing and analysis for 1 Visc and all Amp complete. Return to Rsmall directory (cwd_Rsmall)
                os.chdir(cwd_Rsmall)
            #Processing and analysis complete for all Visc,Amp and 1 Rsmall. Return to cwd_Rlarge to change Rsmall
            os.chdir(cwd_Rlarge)
        #Processing and analysis complete for all Visc,Amp,Rsmall and 1 Rlarge. Return to cwd_STRUCT to change the Rlarge
        os.chdir(cwd_STRUCT)
    #All Processing and analysis complete! 
    #CLOSE FILES   
    file.close()
    end = perf_counter()
    print('Processed {:d} dirs in {:.3f} seconds'.format(len(dirsRlarge), end - start))
    return

#--------------------__END MAIN__------------------------#
main()
