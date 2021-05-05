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

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'

# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    return cwd+"/NOTA_List.txt"

def GetNOTAData(cwd):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    data = data.sort_values(by=['Theta','parHx','parHy'])
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
    data = GetNOTAData(cwd_PYTHON)
    #2 choices: Either ran for full 20s seconds or stopped due to error
    #1: Ran full 20s
    for idx in range(len(data['20s'])):
    #for idx in range(1):
        #Get Parameters Values
        parTheta = data.loc[idx,'Theta']
        parHx = data.loc[idx,'parHx']
        parHy = data.loc[idx,'parHy']
        #Change to the Theta$Theta/Hx$Hx/Hy$Hy directory
        strDir = 'PosData/Theta{0}/Hx{1}/Hy{2}/'.format(parTheta,parHx,parHy)
        cwd_SIM = cwd_PYTHON+'../'+strDir
        if(data.loc[idx,'20s'] == 'y'):
            #Copy pd.txt to the Theta$Theta/Hx$Hx/Hy$Hy directory on laptop as pd2.txt
            os.chdir(cwd_SIM)
            strLL = "tdombro@longleaf.unc.edu:/pine/scr/t/d/tdombro/Spherobot/Pairwise/PairDynamics/LongTime/Re10/Theta{0}/Hx{1}/Hy{2}/pd.txt".format(parTheta,parHx,parHy)
            subprocess.call(["scp",strLL,"pd2.txt"])
            os.chdir(cwd_PYTHON)
        elif(data.loc[idx,'20s'] == 'n'):
            os.chdir(cwd_SIM)
            strLL = "tdombro@longleaf.unc.edu:/pine/scr/t/d/tdombro/Spherobot/Pairwise/PairDynamics/LongTime/Re10/Theta{0}/Hx{1}/Hy{2}/pd.txt".format(parTheta,parHx,parHy)
            subprocess.call(["scp",strLL,"pd2.txt"])
            os.chdir(cwd_PYTHON)

 
