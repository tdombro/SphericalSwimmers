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
def filename(cwd):
    return cwd+"startLData_Re2_.csv"

def GetLData(cwd):
    data = pd.read_csv(filename(cwd),delimiter=' ')
    data['parHy'] *= 2.0
    data['parHx'] *= 2.0
    data = data[data['parThetaBW'] == 112.5].copy()
    data = data.sort_values(by=['parThetaBW','parHx','parHy'])
    data = data.reset_index(drop=True)
    return data

if __name__ == '__main__':
    
    #READ NOTAList.txt to get all sims that did not complete
    #Whether through error or through no end state
    #Pull Hx, Hy, Theta parameters for each
    #Change directory to Theta$Theta/Hx$Hx/Hy$Hy
    #Modify 'script_restart.sh and copy to specified directory
    #Copy input2D_restart into directory
    #Submit with subprocess the command "sbatch script_restart.sh"
    cwd_PYTHON = os.getcwd() + '/'
    data = GetLData(cwd_PYTHON+'../')
    #Restart simulation where it left off. Some at 40s. Some at 20s.
    for idx in range(len(data['endTime'])):
        parTheta = np.round(data.loc[idx,'parThetaBW'],1)
        parHx = np.round(data.loc[idx,'parHx'],1)
        parHy = int(np.round(data.loc[idx,'parHy'],1))
        cwd_VTK = cwd_PYTHON+"../Fluid/Theta{0}/Hx{1}/Hy{2}/VTK/AVG".format(parTheta,parHx,parHy)
        strSBATCH = "sbatch -J AvgVTK_T{0}_Hx{1}_Hy{2}_ -t 1-0 -n 1 -p general -o %x.out --mem-per-cpu=25000 scriptAvg.sh {3} {4} {5}".format(parTheta,parHx,parHy,parTheta,parHx,parHy)
        print(strSBATCH)
        os.system(strSBATCH)
            
