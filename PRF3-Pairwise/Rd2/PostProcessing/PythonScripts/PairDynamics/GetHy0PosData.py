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
ReList = ['2','7','10']
ThetaList = [0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5]
HxList = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5]

if __name__ == '__main__':
    
    #Copy pd.txt for Hy0 sims into separate PosData directory
    cwd_PYTHON = os.getcwd() + '/'
    for Re in ReList:
        cwd_RE = cwd_PYTHON+'Re{0}/'.format(Re)
        for Theta in ThetaList:
            cwd_THETA = cwd_RE+'Theta{0}/'.format(Theta)
            HxDirs = [ name for name in os.listdir(cwd_THETA)
                      if os.path.isdir(os.path.join(cwd_THETA, name)) ]
            for dirHx in HxDirs:
                cwd_POS = cwd_THETA+'{0}/Hy0/'.format(dirHx)
                print('Re{0}: Theta={1}: {2}: Hy0'.format(Re,Theta,dirHx))
                posFile = cwd_POS+'pd.txt'
                strDir = cwd_PYTHON+'PosData/Re{0}/Theta{1}/{2}/Hy0/'.format(Re,Theta,dirHx)
                pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
                newposFile = strDir+'pd.txt'
                os.system('cp '+posFile+' '+newposFile)
