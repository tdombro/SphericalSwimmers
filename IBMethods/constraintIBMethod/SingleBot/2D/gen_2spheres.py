#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#LIB_PATH = 'c:/home/python'
LIB_PATH = '/nas/longleaf/home/hongnt/python'
# /mnt/c/home/python
import os, sys, random, os.path, shutil
sys.path.append(LIB_PATH)
import numpy as np
import utilities as ut

infile = 'input2d'

rLarge = ut.read_from_input(infile, 'RADIUS_LARGE', 'f')
rSmall = ut.read_from_input(infile, 'RADIUS_SMALL', 'f')
cc_distance = ut.read_from_input(infile, 'CCDIST', 'f')
Xlo     = ut.read_from_input(infile, 'Xlo', 'f')
Xhi     = ut.read_from_input(infile, 'Xhi', 'f')
Ylo     = ut.read_from_input(infile, 'Ylo', 'f')
Yhi     = ut.read_from_input(infile, 'Yhi', 'f')
ratio   = ut.read_from_input(infile, 'REF_RATIO', 'f')
level   = ut.read_from_input(infile, 'MAX_LEVELS', 'f')
Nx      = ut.read_from_input(infile, 'Nx', 'f')
Ny      = ut.read_from_input(infile, 'Ny', 'f')
dX      = (Xhi - Xlo)/(Nx * ratio**(level - 1))
dY      = (Yhi - Ylo)/(Ny * ratio**(level - 1))

vertices = []
centerId = []
idx = -1

LargeNx = int(np.floor(rLarge/dX))
SmallNx = int(np.floor(rSmall/dX))
print('SmallNx, LargeNx: ', SmallNx, LargeNx)

for i in range(-LargeNx, LargeNx + 1):
    x = i*dX
    NumPtsY = int(np.floor(np.sqrt(rLarge*rLarge - x*x)/dY))
    for j in range(-NumPtsY, NumPtsY + 1):
        idx += 1
        y = j*dY
        if i == 0 and j == 0:
            centerId.append(idx)
        vertices.append([x,y])

for i in range(-SmallNx, SmallNx + 1):
    x = i*dX
    NumPtsY = int(np.floor(np.sqrt(rSmall*rSmall - x*x)/dY))
    for j in range(-NumPtsY, NumPtsY + 1):
        idx += 1
        y = j*dY
        if i == 0 and j == 0:
            centerId.append(idx)
        vertices.append([x + cc_distance, y])
print('number of vertices:', idx + 1)
print('indices of centers: ', centerId)

f = open('centers.txt','w')
f.write('sphere centers lag ids: %d %d' % (centerId[0], centerId[1]))
f.close()
#ut.plot_vertices_and_bonds(vertices,[],'','')
ut.save_points_file(vertices,'gg','', '.', 'spheres.vertex')
