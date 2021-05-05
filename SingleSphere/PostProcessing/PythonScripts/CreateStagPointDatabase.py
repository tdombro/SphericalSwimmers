#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 11:19:50 2018

@author: thomas
"""

import pandas as pd
import os
import re
from pathlib import Path
from shutil import copyfile
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

#Goal:  1) read all .curve files
#       2) Identify Hori zontal/Vertical Stagnation point locations
#       3) Where the sign(V) changes, the stagnation point is in between those 2 points
#       4) Store all points in database
#       5) Plots stagnation pt vs delta and color code by St

def readFile(direction,Re,St,minDist,Amp):
    #Check to see if file exists
    filePath = Path(os.getcwd()+'/'+direction+'.'+Re+'.'+St+'.curve')
    if(filePath.is_file()):
        print('='*40+'\n')
        print('direction = %s\tRe = %s\tSt = %s' %(direction,Re,St))
        #Name of file
        file = direction+'.'+Re+'.'+St+".curve"
        
        df = pd.read_csv(file,delimiter=' ')
        df = df.rename(index=str, columns={"#": "pos", "curve": "vel"})
        velOld = 0.0
        #Loop over all database entries
        for idx in range(len(df['pos'])):
            #Only continue if pos >= radius
            if(df['pos'][idx] >= minDist):
                velCurrent = df['vel'][idx]
                #print('Current = %.3f\tOld = %.3f'%(velCurrent,velOld))
                if(velCurrent*velOld < 0.0):
                    print('Velocity direction changed!')
                    if(direction == 'Horizontal'):
                        stagPos = (df['pos'][idx] + df['pos'][idx-1])/2.0
                    if(direction == 'Vertical'):
                        stagPos = (df['pos'][idx] + df['pos'][idx-1])/2.0 - Amp
                        #stagPos = df['pos'][idx]
                    print('stagPos = ',stagPos)
                    print('-'*40+'\n')
                    return (stagPos,True)
                velOld = velCurrent
        print('-'*40+'\n')
    return (0.0,False)

def main():
    
    #Parameters for Reading file and extracting info
    LineDirection = ['Horizontal','Vertical']
    discRadius = 0.15
    radius = 0.2
    FREQ = 10.0
    minDist = radius
    figNum = 0
    
    #Parameters for Plotting
    #[Re,St,direction] = (stagPos,velChange?)
    stagPos = np.zeros(2)
    velChange = np.zeros(2)
    horData = pd.DataFrame({'direction': [],'Re': [],'St': [],'stagPos': [],'velChange': []})
    vertData = pd.DataFrame({'direction': [],'Re': [],'St': [],'stagPos': [],'velChange': []})
    LineoutData = [horData, vertData]
    ReynoldsList = np.arange(5.0,16.0,1.0)
    StrouhalList = np.arange(0.25,1.75,0.25)
    colorList = ['red','orange','green','blue','purple','black']
    print(ReynoldsList)
    print(StrouhalList)
    
    #obtain current directory (Python)
    cwd_PYTHON = os.getcwd()
    
    '''This part is for Sweep2 Re = [5.0,15.0,1.0]'''
    #Find all directories when using 'ls' in current directory (Reynolds)
    dirsRe = [d for d in os.listdir(cwd_PYTHON) if os.path.isdir(d)]
    #Single Directory (For testing)
    #dirsRe = ['Re14.0']
    for dirRe in dirsRe:
        #Go into Re directory
        os.chdir(cwd_PYTHON+'/'+dirRe)
        #Save Current Working Directory (Re)
        cwd_Re = os.getcwd()
        Re = dirRe
        #Split 'Re' into string 'Re' and value 'Re'
        ReString, ReValue = tuple(re.split('(\d.*)',Re)[:2])
        #Find all directories when using 'ls' in current directory (St)
        dirsSt = [d for d in os.listdir(cwd_Re) if os.path.isdir(d)]
        #Single Directory (For Testing)
        #dirsSt = ['St1.5']
        for dirSt in dirsSt:
            #Go into St directory
            dirPath = Path(cwd_Re+'/'+dirSt)
            if(dirPath.is_dir()):
                os.chdir(cwd_Re+'/'+dirSt)
            #Save Current Working Directory
            cwd_CURVE = os.getcwd()
            St = dirSt
            #Split 'Re' into string 'Re' and value 'Re'
            StString, StValue = tuple(re.split('(\d.*)',St)[:2])
            #Open CURVE files
            #Obtain Amplitude
            Amp = radius/float(StValue)
            #Loop over 2 directions
            idx = 0
            for direction in LineDirection:
                #Read the curve file for direction 'direction'
                print('direction = %s\tminDist = %.3f'%(direction,minDist))
                stagPos[idx], velChange[idx] = readFile(direction,Re,St,minDist,Amp)
                if(velChange[idx] == True):
                    tempData = pd.DataFrame({'direction': [direction],'Re': [float(ReValue)],
                                             'St': [float(StValue)],'stagPos': [stagPos[idx]],
                                             'velChange': [velChange[idx]]})
                    if(direction == 'Horizontal'):
                        #Add to hor database
                        frames = [horData,tempData]
                        horData = pd.concat(frames)
                    if(direction == 'Vertical'):
                        #Add to vert database
                        frames = [vertData,tempData]
                        vertData = pd.concat(frames)    
                idx += 1
            #return
            #Return to dir 'cwd_Re'
            os.chdir(cwd_Re)
        #Return to dir 'cwd_PYTHON'
        os.chdir(cwd_PYTHON)
        
    '''This part is for Sweep1 Re = [14.0,58.0,4.0]'''
    os.chdir(cwd_PYTHON+'/../VisitFiles/')
    cwd_Sweep1 = os.getcwd()
    print(cwd_Sweep1)
    #Find all directories when using 'ls' in current directory (Reynolds)
    dirsRe = [d for d in os.listdir(cwd_Sweep1) if os.path.isdir(d)]
    #Single Directory (For testing)
    #dirsRe = ['Re14.0']
    for dirRe in dirsRe:
        #Go into Re directory
        os.chdir(cwd_Sweep1+'/'+dirRe)
        #Save Current Working Directory (Re)
        cwd_Re = os.getcwd()
        Re = dirRe
        #Split 'Re' into string 'Re' and value 'Re'
        ReString, ReValue = tuple(re.split('(\d.*)',Re)[:2])
        #Find all directories when using 'ls' in current directory (St)
        dirsSt = [d for d in os.listdir(cwd_Re) if os.path.isdir(d)]
        #Single Directory (For Testing)
        #dirsSt = ['St1.5']
        for dirSt in dirsSt:
            #Go into St directory
            dirPath = Path(cwd_Re+'/'+dirSt)
            if(dirPath.is_dir()):
                os.chdir(cwd_Re+'/'+dirSt)
            #Save Current Working Directory
            cwd_CURVE = os.getcwd()
            St = dirSt
            #Split 'Re' into string 'Re' and value 'Re'
            StString, StValue = tuple(re.split('(\d.*)',St)[:2])
            #Open CURVE files
            #Obtain Amplitude
            Amp = radius/float(StValue)
            #Loop over 2 directions
            idx = 0
            for direction in LineDirection:
                #Read the curve file for direction 'direction'
                stagPos[idx], velChange[idx] = readFile(direction,Re,St,minDist,Amp)
                if(velChange[idx] == True):
                    tempData = pd.DataFrame({'direction': [direction],'Re': [2.0*float(ReValue)],
                                             'St': [float(StValue)],'stagPos': [stagPos[idx]],
                                             'velChange': [velChange[idx]]})
                    if(direction == 'Horizontal'):
                        #Add to hor database
                        frames = [horData,tempData]
                        horData = pd.concat(frames)
                    if(direction == 'Vertical'):
                        #Add to vert database
                        frames = [vertData,tempData]
                        vertData = pd.concat(frames)    
                idx += 1
            #return
            #Return to dir 'cwd_Re'
            os.chdir(cwd_Re)
        #Return to dir 'cwd_PYTHON'
        os.chdir(cwd_PYTHON)
    
    print(horData)
    print('='*40+'\n')
    print(vertData)    
    print('='*40+'\n')
    
    LineoutData = [horData, vertData]
    for data in LineoutData:
        data = data.sort_values(by=['Re','St'])
        data = data.rename(index=str,columns={'Re':'Re_s'})
        data = data.reset_index(drop=True)
        #Adding Variables to Data
        data['Amp'] = discRadius/data.St
        data['Nu'] = 2.0*np.pi*FREQ*data.Amp**2/data.Re_s
        data['delta'] = np.sqrt(data.Nu/(2.0*np.pi*FREQ))
        data['Re_Paper'] = data.Amp*discRadius/data.delta**2
        data['epsilon'] = data.Amp/discRadius
        data['Tvar'] = data.stagPos - discRadius
        data['Msquared'] = data.Re_Paper/data.epsilon
        data['Re_M'] = discRadius**2/data.delta**2
        print(data[data.St == 0.25])
        directionList = data.direction.tolist()
        #Create Database of 'data'
        data.to_csv('StagPoint_r0.15m_'+directionList[0]+'.csv')
        
    return
    
#------------------------__END_MAIN__-------------------------------------
main()