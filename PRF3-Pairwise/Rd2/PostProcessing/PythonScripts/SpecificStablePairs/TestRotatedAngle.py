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
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
import pathlib
from matplotlib.colors import Normalize
from scipy import interpolate
norm = Normalize()
from resource import getrusage, RUSAGE_SELF

def Rotate(xy, theta):
    # https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions
    #First Rotate based on Theta
    #Allocate Arrays
    rotationMatrix = np.zeros((2,2))
    #Calculate rotation matrix
    rotationMatrix[0,0] = np.cos(theta)
    rotationMatrix[0,1] = -1.0*np.sin(theta)
    rotationMatrix[1,0] = np.sin(theta)
    rotationMatrix[1,1] = np.cos(theta)
    return rotationMatrix.dot(xy)

def CalcLabAngle(pos):
    #Find swimming axis (normal y-axis)
    xU, xL = pos.loc[0,'aXU'], pos.loc[0,'aXL']
    yU, yL = pos.loc[0,'aYU'], pos.loc[0,'aYL']
    labX   = xU - xL
    labY   = yU - yL
    length = np.hypot(labX,labY)
    normY = labY/length
    #2) Calculate Theta
    if(normY >= 0.0):
        theta = np.arccos(normY)
    else:
        theta = -1.0*np.arccos(normY)+2.0*np.pi
    print('theta = ',theta*180.0/np.pi)
    #return 2.0*np.pi - theta
    return theta

#Test Angle
xy = np.array([1,-1])
norm = xy/np.hypot(xy[0],xy[1])
if(norm[0] <= 0.0):
    theta = np.arccos(norm[1]) #- np.pi/2.0
else:
    theta = -1.0*np.arccos(norm[1])+2.0*np.pi #- np.pi/2.0
print('theta = ',theta*180.0/np.pi)
theta_rotate = 2.0*np.pi - theta
#theta_rotate = theta
print('theta_rot = ',theta_rotate*180.0/np.pi)
print('xy = ',xy)
#theta = 5.0*np.pi/4.0
xy_rot = Rotate(xy,theta_rotate)
print('xy_rot = ',xy_rot)

#Test Field Interpolation
Ux = np.random.rand(1024,1024)
Uy = np.random.rand(1024,1024)
x = np.linspace(-0.05,0.05,1024)
y = np.linspace(-0.05,0.05,1024)
mx,my = np.meshgrid(x,y)
theta_rotate = np.pi/8.0
mxy = np.array([mx.flatten(),my.flatten()])
mxy_rot = np.zeros((2,1024*1024))
for jdx in range(1024*1024):
    mxy_rot[:,jdx] = Rotate(mxy[:,jdx],theta_rotate)
mx_rot = mxy_rot[0,:].reshape((1024,1024))
my_rot = mxy_rot[1,:].reshape((1024,1024))
x_new = np.linspace(-0.025,0.025,512)
y_new = np.linspace(-0.025,0.025,512)
mx_new, my_new = np.meshgrid(x_new,y_new)

arrayUx_new=interpolate.griddata((mx_rot.flatten(),my_rot.flatten()),Ux.flatten() , (mx_new,my_new),method='cubic')
print('X transformation complete')
print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)
arrayUy_new=interpolate.griddata((mx_rot.flatten(),my_rot.flatten()),Uy.flatten() , (mx_new,my_new),method='cubic')
print('Y tranformation complete')
print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)




