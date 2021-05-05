#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 11:26:03 2018

@author: thomas
"""

import pandas as pd
import numpy as np
import os
import re
from pathlib import Path
from shutil import copyfile
import matplotlib.pyplot as plt
from scipy import stats

#CONSTANTS
radiusLarge = 0.30 #Radius of large sphere
radiusSmall = 0.15 #Radius of small sphere
RSL = 0.75 #Default resting spring length
FREQ = 10.0 #frequency of oscillation
AMPFRAC = 0.8 #Fraction of max amplitude
MAXAMP = 0.30 #maximum amplitude of oscillation
DENSITY = 2.0 #Fluid density (kg/m-3)
L = 0.3 #Length Scale (Currently just radius of large sphere)
ssFrac = 2.0
ALPHA = radiusLarge/radiusSmall

def main():
    
    #Read and import WholeSweep data file
    #Columns used by data.'NameofColumnHere'
    data = pd.read_csv('WholeSweep/FinalSweepData.txt',delimiter = ' ')
    data = data.sort_values(by=['visc','ampFrac','rslFrac'])
    #Remove data (Would go here below)
    #BAD DATA POINTS
    
    #RSL1.0: amp=0.8, visc <= 0.05
    #RSL>=1.1: amp=0.8, visc <= 0.02
    #RSL2.0: amp=0.7, visc = 0.01
    
    #Tom Data Removal
    data = data.ix[~((data.rslFrac == 1.0) & (data.ampFrac == 0.8) & (data.visc <= 0.05))]
    data = data.ix[~((data.rslFrac >= 1.1) & (data.ampFrac == 0.8) & (data.visc <= 0.02))]
    data = data.ix[~((data.ampFrac == 0.7) & (data.rslFrac >= 1.6) & (data.visc <= 0.01))]
    #Regime II
    #data = data[(data.visc > 0.01) & (data.visc >= 0.02) & (data.visc <=0.06)]
    #Regime I
    
    #Viscosity range
    data = data[(data.visc >= 0.01) & (data.visc < 2000.08)]# & (data.visc < 10.0)]
    #Amplitude range
    data = data[data.ampFrac >= 0.4]
    #RSL range
    #data = data[data.rslFrac <= 1.6]
    #Velocity range
    #data = data[(data.vel < 0.15)]# & (data.vel >= -0.03)]
    
    print(data)
    #Add new column to data by the following: data['NewColumnName'] = Expression
    data['rsl'] = data.rslFrac*RSL
    data['amp'] = data.ampFrac*MAXAMP
    data['visc'] = data.visc/DENSITY #Converts mu to nu 
    data['cd'] = data.rsl - data.amp - radiusLarge - radiusSmall #Closest distance the 2 spheres get to each other
    
    #Create dimensionless parameters
    data['rsldivamp'] = data.rsl/data.amp
    data['delta'] = np.sqrt(data.visc/(2.0*np.pi*FREQ))
    data['VdivFrsl'] = data.vel/(FREQ*data.rslFrac)
    data['VdivFamp'] = data.vel/(FREQ*data.ampFrac)
    data['VdivFdelta'] = data.vel/(FREQ*data.delta)
    data['VdivF'] = data.vel/FREQ
    data['ampdivdelta'] = data.ampFrac/data.delta
    data['rsldivdelta'] = data.rslFrac/data.delta #Prob not gonna be used
    data['surf'] = data.rsl - radiusSmall - radiusLarge
    data['length'] = data.rsl + radiusSmall + radiusLarge 
    data['AmpSmall'] = data.amp*0.8
    data['AmpLarge'] = data.amp*0.2
    data['Rsmall'] = radiusSmall
    data['Rlarge'] = radiusLarge
    
    #Dimensionless Quantities
    data['Re'] = data.AmpSmall*data.Rsmall/data.delta**2
    data['Re_s'] = data.AmpSmall**2/data.delta**2
    print('minRe = %.4f\tmaxRe = %.4f'%(np.amin(data.Re),np.amax(data.Re)))
    data['M_R'] = data.Rlarge/data.delta
    data['M_r'] = data.Rsmall/data.delta
    data['A_RdivR'] = data.AmpLarge/data.Rlarge
    data['A_rdivr'] = data.AmpSmall/data.Rsmall
    data['eps_R'] = data.AmpLarge/data.Rlarge
    data['eps_r'] = data.AmpSmall/data.Rsmall
    data['A_RdivRSL'] = data.AmpLarge/data.rsl
    data['A_rdivRSL'] = data.AmpSmall/data.rsl
    
    #Import SphereSize Ratio Data
    dataSphere = pd.read_csv('SphereSize/Shannon/NewSphereSizeData.txt',delimiter = ' ') #Large
    dataSphere = dataSphere.sort_values(by=['Nu','Rlarge','Rsmall']) #Large
    dataSphere = dataSphere[dataSphere.Nu > 0.005] #Large
    dataSphere['delta'] = np.sqrt(dataSphere.Nu/(2.0*np.pi*FREQ)) #Large
    dataSphere['rsl'] = 1.3*RSL
    dataSphere['length'] = dataSphere.rsl + dataSphere.Rlarge + dataSphere.Rsmall
    dataSphere['sphereRatio'] = dataSphere.Rlarge/dataSphere.Rsmall #Large
    dataSphere['massRatio'] = dataSphere.sphereRatio**2
    dataSphere['amp'] = dataSphere.ampFrac*MAXAMP
    dataSphere['AmpSmall'] = dataSphere.amp*(dataSphere.massRatio/(dataSphere.massRatio + 1.0))
    dataSphere['AmpLarge'] = dataSphere.amp*(1.0/(dataSphere.massRatio + 1.0))
    
    #Dimensionless Quantities
    dataSphere['Re'] = dataSphere.AmpSmall*dataSphere.Rsmall/data.delta**2
    dataSphere['Re_s'] = dataSphere.AmpSmall**2/dataSphere.delta**2
    print('minRe = %.4f\tmaxRe = %.4f'%(np.amin(dataSphere.Re),np.amax(dataSphere.Re)))
    dataSphere['M_R'] = dataSphere.Rlarge/dataSphere.delta
    dataSphere['M_r'] = dataSphere.Rsmall/dataSphere.delta
    dataSphere['A_RdivR'] = dataSphere.AmpLarge/dataSphere.Rlarge
    dataSphere['A_rdivr'] = dataSphere.AmpSmall/dataSphere.Rsmall
    dataSphere['eps_R'] = dataSphere.AmpLarge/dataSphere.Rlarge
    dataSphere['eps_r'] = dataSphere.AmpSmall/dataSphere.Rsmall
    dataSphere['A_RdivRSL'] = dataSphere.AmpLarge/dataSphere.rsl
    dataSphere['A_rdivRSL'] = dataSphere.AmpSmall/dataSphere.rsl
    
    #Calculate Minimum and Maximum values of dimensional parameters
    #Dimensional Variables: delta, Rsmall, Rlarge, AmpSmall, AmpLarge, RSL
    #Dimensionless Variables: Re, Re_s, M_R, M_r, A_R/R, A_r/r, eps_R, eps_r, A_R/RSL, A_r/RSL
    
    #Large Sweep
    #Viscosity range
    data = data[(data.visc >= 0.01) & (data.visc < 2000.08)]
    #Amplitude range
    data = data[data.ampFrac >= 0.4]
    #RSL range
    #data = data[data.rslFrac <= 1.6]
    #Velocity range
    #data = data[(data.vel < 0.15)]# & (data.vel >= -0.03)]
    
    #Sphere Sweep
    #Viscosity range
    #dataSphere = dataSphere[(dataSphere.Nu >= 0.01) & (dataSphere.Nu < 2000.08)]
    #Amplitude range
    #dataSphere = dataSphere[dataSphere.ampFrac >= 0.4]
    #RSL range
    #dataSphere = dataSphere[dataSphere.rsl <= 0.75]
    #Velocity range
    #dataSphere = dataSphere[(dataSphere.vel < 0.15)]# & (dataSphere.vel >= -0.03)]
    
    f = open('RangeOfSweeps.txt','w')
    f.write('='*40+'\n')
    f.write('Start of Large Sweep\n')
    f.write('='*40+'\n\n')
    varNames = ['delta','Rsmall','Rlarge','AmpSmall','AmpLarge','rsl','Re','Re_s','M_R',
                'M_r','A_RdivR','A_rdivr','eps_R','eps_r','A_RdivRSL','A_rdivRSL']
    #LargeSweep
    for Name in varNames:
        minVar = np.amin(data[Name])
        maxVar = np.amax(data[Name])
        varString = Name+': ['+str(round(minVar,3))+', '+str(round(maxVar,3))+']'
        f.write(varString+'\n')
    f.write('\n')
    f.write('='*40+'\n')
    f.write('End of Large Sweep\n')
    f.write('='*40+'\n')
    f.write('Start of Sphere Sweep\n')
    f.write('='*40+'\n\n')
    #SphereSweep
    for Name in varNames:
        minVar = np.amin(dataSphere[Name])
        maxVar = np.amax(dataSphere[Name])
        varString = Name+': ['+str(round(minVar,3))+', '+str(round(maxVar,3))+']'
        f.write(varString+'\n')
    f.write('\n')
    f.write('='*40+'\n')
    f.write('End of Sphere Sweep\n')
    f.write('='*40+'\n')
    f.close()
    return

#------------------__END MAIN__-----------------------------------
main()