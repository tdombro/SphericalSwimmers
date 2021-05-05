#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 13:49:41 2020

@author: thomas
"""

import os, sys
import pandas as pd
import time as t
import pathlib

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
ReList=['0.5','0.6','0.7','0.8','0.9','1.0','2.0','3.0','4.0','5.0','5.5','6.0','6.5','7.0','7.5',
          '10.0','12.5','15.0','17.5','20.0','25.0','30.0','35.0','40.0','50.0','60.0']

# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    #return cwd+"/pd.txt"
    #cwd = cwd_PYTHON
    return cwd+"/pd.txt"

def GetPosDataLength(cwd):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    return int((len(data['time'])-1)/2)

#Find missing VTK files for V-shape sims
print('Finding V-shape missing VTKs')
fp = open(cwd_PYTHON+'missingVTK.txt','w')
for Re in ReList:
    cwd_Re = cwd_PYTHON+'../SweepRe/Re'+Re+'/'
    nDumps = GetPosDataLength(cwd_Re)
    cwd_DATA = cwd_Re+'/VTK/'
    missingVTK = []
    for idx in range(nDumps):
        file = pathlib.Path(cwd_DATA+'DATA%05d.vtk'%idx)
        if not file.exists ():
            print ("File %i does not exist"%idx)
            missingVTK.append(idx)
            fp.write('%i '%idx)
    print('Missing VTK: Re = '+Re)
    print(missingVTK)
    fp.write('\n')
fp.close()
