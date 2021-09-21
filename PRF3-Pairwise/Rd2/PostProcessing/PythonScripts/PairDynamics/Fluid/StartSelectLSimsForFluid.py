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
        #Copy files from PD Sim directory to a new fluid directory
        cwd_FLUID = cwd_PYTHON+"Fluid/Theta{0}/Hx{1}/Hy{2}/".format(parTheta,parHx,parHy)
        pathlib.Path(cwd_FLUID).mkdir(parents=True, exist_ok=True)
        strDir = 'Theta{0}/Hx{1}/Hy{2}/'.format(parTheta,parHx,parHy)
        cwd_SIM = cwd_PYTHON+strDir
        strCopy = 'cp '+cwd_SIM+'bot* '+cwd_FLUID
        os.system(strCopy)
        strCopy = 'cp '+cwd_SIM+'skel* '+cwd_FLUID
        os.system(strCopy)
        copyfile('input2D_fluid',cwd_FLUID+'input2D_fluid')
        
        #Change to the Theta$Theta/Hx$Hx/Hy$Hy directory

        #Create Executable
        strTemplate = cwd_PYTHON+'Template/main2d'
        strSIM = cwd_FLUID+'main2d.exe'
        strLink = 'ln -s ' + strTemplate+' ' + strSIM
        file = pathlib.Path(strSIM)
        if not file.exists():
            os.system(strLink)
            print('Re2: Linked: Theta={0}: Hx={1}: Hy={2}'.format(parTheta,parHx,parHy))
        #Create sbatch script for fluid sims
        f = open('script_fluid.sh','w')
        f.write('#!/bin/sh\n')
        strSBATCH1 = '#SBATCH --job-name=LFluid_Theta{0}_Hx{1}_Hy{2}_Re2\n'.format(parTheta,parHx,parHy)
        f.write(strSBATCH1)
        f.write('#SBATCH --ntasks=1\n')
        f.write('#SBATCH --time=7-0\n')
        f.write('#SBATCH --partition=general\n\n')
        strSBATCH2 = './main2d.exe input2D_fluid'
        f.write(strSBATCH2)
        f.close()
        
        copyfile('script_fluid.sh',cwd_FLUID+'script_fluid.sh')
        print('Re2: Submitting Fluid Simulation: Theta={0}: Hx={1}: Hy={2}'.format(parTheta,parHx,parHy))
        os.chdir(cwd_FLUID)
        subprocess.call(["sbatch","script_fluid.sh"])
        os.chdir(cwd_PYTHON)
            
