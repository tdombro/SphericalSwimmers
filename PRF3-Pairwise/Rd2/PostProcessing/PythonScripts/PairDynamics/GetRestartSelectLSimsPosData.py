#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 21:59:46 2020

@author: thomas
"""

import numpy as np
import pandas as pd
import os, sys
import time as t
import subprocess
from shutil import copyfile
import pathlib

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'

# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    return cwd+"startLData_Re2_.csv"

def GetLData(cwd):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    data['parHy'] *= 2.0
    data['parHx'] *= 2.0
    data = data[data['parThetaBW'] <= 180.0].copy()
    data = data.sort_values(by=['parThetaBW','parHx','parHy'])
    data = data.reset_index(drop=True)
    return data

def GetRestartIndex(cwd):
    file = pathlib.Path(cwd+'pd.txt')
    if file.exists():
        #print('file_exists')
        data = pd.read_csv(cwd+'pd.txt',delimiter=' ')
        UAdata = data[data['idx'] == 6].copy()
        #Sort data by time and reset indices
        UAdata = UAdata.sort_values(by=['time'])
        UAdata = UAdata.reset_index(drop=True)
        lastData = UAdata.tail(1)
        lastData = lastData.reset_index(drop=True)
        endTime = lastData.loc[0,'time']
        endTime = int(np.trunc(endTime))
        return int(endTime*1e5)
    else:
        return 0

if __name__ == '__main__':
    
    #READ NOTAList.txt to get all sims that did not complete
    #Whether through error or through no end state
    #Pull Hx, Hy, Theta parameters for each
    #Change directory to Theta$Theta/Hx$Hx/Hy$Hy
    #Modify 'script_restart.sh and copy to specified directory
    #Copy input2D_restart into directory
    #Submit with subprocess the command "sbatch script_restart.sh"
    cwd_PYTHON = os.getcwd() + '/'
    data = GetLData(cwd_PYTHON)
    #Restart simulation where it left off. Some at 40s. Some at 20s.
    for idx in range(len(data['endTime'])):
        parTheta = np.round(data.loc[idx,'parThetaBW'],1)
        parHx = np.round(data.loc[idx,'parHx'],1)
        parHy = int(np.round(data.loc[idx,'parHy'],1))
        #Find restart interval from pd.txt data
        cwd_POS = cwd_PYTHON+'Theta{0}/Hx{1}/Hy{2}/'.format(parTheta,parHx,parHy)
        restartIndex = GetRestartIndex(cwd_POS)
        print('Re2: Theta={0}: Hx={1}: Hy={2}: time = {3}'.format(parTheta,parHx,parHy,restartIndex/1e5))
        #Copy pd.txt as pd3.txt over to another directory
        strDir = cwd_PYTHON+'PosData/Re2/Theta{0}/Hx{1}/Hy{2}/'.format(parTheta,parHx,parHy)
        pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
        posFile = cwd_POS+'pd.txt'
        newposFile = strDir+'pd3.txt'
        os.system('cp '+posFile+' '+newposFile)
