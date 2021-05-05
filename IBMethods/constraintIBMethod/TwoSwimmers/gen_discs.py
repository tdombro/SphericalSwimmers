#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#LIB_PATH = 'c:/home/python'
LIB_PATH = '/nas/longleaf/home/hongnt/python'
#LIB_PATH = '/mnt/c/home/python'
import os, sys, random, os.path, shutil
sys.path.append(LIB_PATH)
import numpy as np
import utilities as ut
from CircleClass import *
from InputClass import *
from CentersClass import *

infile  = 'input2d'
rLarge  = ut.read_from_input(infile, 'RADIUS_LARGE', 'f')
rSmall  = ut.read_from_input(infile, 'RADIUS_SMALL', 'f')
Xlo     = ut.read_from_input(infile, 'Xlo', 'f')
Xhi     = ut.read_from_input(infile, 'Xhi', 'f')
Ylo     = ut.read_from_input(infile, 'Ylo', 'f')
Yhi     = ut.read_from_input(infile, 'Yhi', 'f')
ratio   = ut.read_from_input(infile, 'REF_RATIO', 'f')
level   = ut.read_from_input(infile, 'MAX_LEVELS', 'f')
Nx      = ut.read_from_input(infile, 'Nx', 'f')
Ny      = ut.read_from_input(infile, 'Ny', 'f')
DS      = ut.read_from_input(infile, 'DS', 'f')
dX      = (Xhi - Xlo)/(Nx * ratio**(level - 1))
dY      = (Yhi - Ylo)/(Ny * ratio**(level - 1))
PI      = ut.read_from_input(infile, 'PI', 'f')
nbot    = ut.read_from_input(infile, 'NBOT', 'i')
radii   = [rLarge, rSmall]

if (dX != dY):
    print("ERROR: dX, dY not the same!!")
    quit()

gen_bot_input = ut.read_from_input(infile, 'gen_bot_input', 's')
trans_mom  = ut.read_from_input(infile, 'TRANS_MOM', 'af')
rot_mom    = ut.read_from_input(infile, 'ROT_MOM', 'af')
ib_update  = ut.read_from_input(infile, 'IB_UPDATE', 's')


centers = []
allVertices  = []
struct_names = []

def twospheres_square_mesh(centers, radii, ds):
    # input: center of two spheres [center1, center2]
    # radii = [radius1, radius2]
    # ds = mesh spacing
    if len(centers) != 2 or len(radii) != 2:
        print('Error::length of centers or radii not equal 2', centers)
        return
    print('ds = %g, radii = [%g, %g]' % (ds, radii[0], radii[1]))

    large = centers[0]
    small = centers[1]

    lo = ut.distance(large, small)
    cc = [small[0] - large[0], small[1] - large[1]]
    angle = ut.angleRad(cc, [1, 0])
    if cc[1] < 0: angle = 2*PI - angle
    print('bot length = %g, angle (small to large) wrt to x-axis = %g' % (lo, angle * 180/PI))

    Nx = []
    for radius in radii:
        Nx.append(int(np.round(radius/ds)))

    #print('SmallNx, LargeNx: ', SmallNx, LargeNx)

    idx = -1
    vertices = []
    npoints = []
    centerId = []
    for discId in range(2):
        npts = 0
        nx = Nx[discId]
        radius = radii[discId]

        for i in range(-nx, nx + 1):
            x = i*ds
            NumPtsY = int(np.floor(np.sqrt(abs(nx*nx - i*i))))

            for j in range(-NumPtsY, NumPtsY + 1):
                idx += 1
                npts += 1
                y = j*ds

                if i == 0 and j == 0:
                    centerId.append(idx)
                vertices.append([x + discId*lo,y])
        npoints.append(npts)

    print('accepted radii = [%g, %g]' % (Nx[0]*ds, Nx[1]*ds))
    print('npoints along diameter = [%d, %d] ' % (2*Nx[0] + 1, 2*Nx[1] + 1))
    print('total vertices = %d, [large = %d, small = %d]' % (len(vertices), npoints[0], npoints[1]))

    ut.rotate2d(vertices, angle)
    ut.displace(vertices, large)
    return Nx, centerId, vertices

def twospheres_polar_mesh(centers, radii, ds):
    if len(centers) != 2 or len(radii) != 2:
        print('Error::length of centers or radii not equal 2', centers)
        return
    print('ds = %g, radii = [%g, %g]' % (ds, radii[0], radii[1]))

    large = centers[0]
    small = centers[1]

    lo = ut.distance(large, small)
    cc = [small[0] - large[0], small[1] - large[1]]
    angle = ut.angleRad(cc, [1, 0])
    if cc[1] < 0: angle = 2*PI - angle
    print('bot length = %g, angle (small to large) wrt to x-axis = %g' % (lo, angle * 180/PI))

    Nx = []
    for radius in radii:
        Nx.append(int(np.round(radius/ds)))

    #print('SmallNx, LargeNx: ', SmallNx, LargeNx)

    idx = 0
    vertices = []
    npoints = []
    centerId = []
    CENTERS = [large, small]
    for discId in range(2):
        center = CENTERS[discId]
        nx = Nx[discId]
        rValues = [ds * i for i in range(1, nx + 1)]
        toDEG = 180 / PI

        # the first vertex is the center
        vertices.append([0.0 + lo*discId, 0.0])
        centerId.append(idx)
        npts = 1
        idx += 1

        for r in rValues:
            nPhi = int(np.round(0.5*PI * r / ds))
            dPhi = 0.5*PI / nPhi
            phiValues = [dPhi * i for i in range(-2*nPhi, 2*nPhi)]
            #print('\tr = %10g, nPhi = %4d,  phi = [%5.1f, %5.1f]' %
            #      (r, nPhi, min(phiValues)*180/PI, max(phiValues)*180/PI))
            for phi in phiValues:
                idx += 1
                npts += 1
                x = r * np.cos(phi)
                y = r * np.sin(phi)
                vertices.append([x + lo*discId, y])

        npoints.append(npts)

    print('accepted radii = [%g, %g]' % (Nx[0]*ds, Nx[1]*ds))
    print('npoints along diameter = [%d, %d] ' % (2*Nx[0] + 1, 2*Nx[1] + 1))
    print('total vertices = %d, [large = %d, small = %d]' % (len(vertices), npoints[0], npoints[1]))

    ut.rotate2d(vertices, angle)
    ut.displace(vertices, large)
    return Nx, centerId, vertices

def get_centers(*argc):
    gen_bot_input = argc[0]
    nbot          = argc[1]
    radii         = argc[2]

    inp = InputDatabase(gen_bot_input)
    gen_method = inp.getString('gen_method')
    subdb = inp.getDatabase(gen_method)
    origin = inp.getDoubleArray('origin')

    if gen_method == 'RANDOM':
        limit_box = subdb.getDoubleArray('limit_box')
        # create list of bot length
        list_length = []
        for botId in range(nbot):
            list_length.append(subdb.getDouble('rand_cc_dist'))
        # create a random position in the limit box
        c = Centers(origin, nbot, radii, list_length)
        c.genRandomRectangle(limit_box)
        list_centers = c.getTwoCenters()
        del c
    elif gen_method == 'MANUAL':
        list_centers = []
        list_length = []
        for botId in range(nbot):
            lpos = subdb.getDoubleArray(str(botId + 1) + 'L')
            spos = subdb.getDoubleArray(str(botId + 1) + 'S')
            list_centers.append([lpos, spos])
            list_length.append(ut.distance(lpos, spos))

    elif gen_method == 'LATTICE':
        limit_box = subdb.getDoubleArray('limit_box')
        nrow = subdb.getInt('num_row')
        ncol = subdb.getInt('num_col')
        rand_updown = subdb.getBool('rand_updown')

        if (nrow*ncol) != nbot:
            print('gen lattice: invalid number of row(%d) or col(%d)'%(nrow, ncol))
            quit()
        # create list of bot length
        list_length = []
        for botId in range(nbot):
            list_length.append(subdb.getDouble('latt_cc_dist'))
        # create a random position in the limit box
        c = Centers(origin, nbot, radii, list_length)
        c.genUniformRectangle(limit_box, nrow, ncol, rand_updown)
        list_centers = c.getTwoCenters()
        del c
    else:
        print('gen_method %s not supported' %(gen_method))
        quit()
    return list_centers

def main(*argc):
    centers = argc[0]
    radii   = argc[1]
    DS      = argc[2]
    dX      = argc[3]
    allVertices     = argc[4]
    struct_names    = argc[5]

    lag_idx = 0
    for bid in range(len(centers)):
        print("bot %d, mesh POLAR" % (bid + 1))
        botname = 'disc_' + str(bid + 1)
        struct_names.append(botname)

        #nx, centId, vertices = twospheres_square_mesh(centers[bid], [rLarge, rSmall], DS*dX)
        nx, centId, vertices = twospheres_polar_mesh(centers[bid], radii, DS*dX)

        print("Lag center idx = [%d, %d]" % (lag_idx + centId[0], lag_idx + centId[1]))
        allVertices += vertices
        ut.save_points_file(vertices, 'gg', '', '.', botname + '.vertex')
        lag_idx += len(vertices)
        print("="*60)


    for bid in range(len(centers)):
        [l, s] = [centers[bid][0], centers[bid][1]]

        print('bot %d: [%g, %g, %g, %g]'% (bid, l[0], l[1], s[0], s[1]))

    print("="*60)
    return # main

# get centers of all spherobots
centers = get_centers(gen_bot_input, nbot, radii)

main(centers, radii, DS, dX, allVertices, struct_names)

ut.plot_vertices_list([allVertices], 2)
ut.constraintib_database('struct_info.txt', struct_names, trans_mom, rot_mom, ib_update, centers)
