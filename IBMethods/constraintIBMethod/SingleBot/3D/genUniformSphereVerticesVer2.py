#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 11:33:33 2018

@author: thomas
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import os, sys
import zipfile
from shutil import copyfile
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import pathlib

def MakeDirectory(directory,seed):
    if not os.path.exists(directory+'/'+str(seed)):
        os.makedirs(directory+'/'+str(seed))
    return

def zipFiles(src,dst):
    zf = zipfile.ZipFile('%s.zip' % (dst), 'w', zipfile.ZIP_DEFLATED)
    abs_src = os.path.abspath(src)
    for dirname, subdirs, files in os.walk(src):
        for filename in files:
            absname = os.path.abspath(os.path.join(dirname, filename))
            arcname = absname[len(abs_src) + 1:]
            print('zipping %s as %s' % (os.path.join(dirname, filename),
                                        arcname))
            zf.write(absname, arcname)
    zf.close()
    return

def GenPolarMesh(nLayers,dS,vertices,origin):
    count = 1
    for jdx in range(1,nLayers):
        #Ring Parameters
        radius = jdx*dS
        bandLength = 2.0*np.pi*radius                   #Circumference of band (m)
        npts = int(np.ceil(bandLength/dS))             #Number of vertex points along band
        print('i = %i\tnpts = %i\tradius = %.5e' %(jdx,npts,radius))
        dtheta = 2.0*np.pi/npts                         #Change in angle between vertex points
        for kdx in range(npts):
            #Vertex coordinate (including shift)
            x = radius*np.cos(kdx*dtheta) + origin[0]
            y = radius*np.sin(kdx*dtheta) + origin[1]
            vertices.append([x,y])
            count += 1
    return (vertices, count)

def genPolarSpheres(nShells,dS,vertices,origin):
    count = 1
    for i in range(nShells):
        #Ring Parameters
        radius = (i+1)*dS
        nTheta = m.ceil(np.pi*radius/dS)+1
        for j in range(nTheta+1):
            Theta = (j/nTheta)*np.pi
            rho = radius*np.sin(Theta)
            bandLength = 2.0*np.pi*rho  #Circumference of x,y plane (m)
            nPhi = m.ceil(bandLength/dS)  #Number of vertex points along ring
            #npts_shell += nPhi
            print('i = %i\tnPhi = %i' %(i,nPhi))
            #print('i = %i\tnTheta = %i'%(i,nTheta))
            #print('i = %i\tnpts_shell = %i' %(i,npts_shell))
            for k in range(nPhi):
                Phi = (k/nPhi)*2.0*np.pi
                #Vertex coordinate (including shift)
                x = rho*np.cos(Phi) + origin[0]
                y = rho*np.sin(Phi) + origin[1]
                z = radius*np.cos(Theta) + origin[2]
                vertices.append([x,y,z])
                count += 1
    print('count = ',count)
    return (vertices, count)

def genDiscMesh(Radii,centers,dX, dirName,FACTOR):
    #INPUTS:
    #dX = fluid grid spacing.
    #dS = 0.5dX (Standard); dS = 2.0dX (CIB)
    #definedRadius = Radius of sphere
    #origin = center of sphere
    #structureName = name of sphere
    #dirName = Where files are saved
    
    #Between Spheres
    maxRadius = Radii
    vertices = []
    dS_sphere = dX
    count = 0
    nShells = (maxRadius - dS_sphere)/dS_sphere + 1
    print(nShells)
    dX_new = dX*nShells/int(np.ceil(nShells))
    nShells = int(np.ceil(nShells))
    print(nShells)
    dS = FACTOR*dX_new
    nLayers = int(np.ceil(nShells/(dS/dS_sphere)))
    #Create center vertices
    origin = centers
    vertices.append(origin)
    count = 1
    vertices, count = genPolarSpheres(nLayers,dS,vertices,origin)
    nVert = count
    
    print(nVert)
    print('-'*40)
    return (nVert,vertices)
    
def WriteVertexData(dirName,structureName,data,nUpper,nLower):
    nPoints = len(data[:,0])
    f = open(dirName+structureName+'.vertex','w')
    f.write('%i #%i large %i small\n' % (nPoints,nUpper,nLower))
    for idx in range(nPoints):
        f.write('%.5e %.5e %.5e\n' %(data[idx,0],data[idx,1],data[idx,2]))
    f.close()
    return

def PlotVertexData(dirName,structureName,data,nUpper):
    #Plot parameters
    fig1 = plt.figure(num=1,figsize=(10,5),dpi=120)
    ax1 = fig1.add_subplot(111,projection='3d')
    ax1.scatter(data[:,0], data[:,1], data[:,2] ,s=1,c='r')
    ax1.scatter(data[0,0],data[0,1], data[0,2], s=1,c='k',zorder=5)
    ax1.scatter(data[nUpper,0],data[nUpper,1],data[nUpper,2], s=1,c='k',zorder=5)
    #ax1.axis([-0.8,0.8,0.0,0.4,-0.4,0.4])
    ax1.set_xlim(-0.8,0.8)
    ax1.set_ylim(-0.4,0.4)
    ax1.set_zlim(-0.4,0.4)
    ax1.set_aspect('equal')
    fig1.savefig(dirName+structureName+'.png',bbox_inches='tight',pad_inches=0.5)
    fig1.clf()
    #plt.show()
    plt.close()
    return

def main():
    #Goal: Generate a uniform 2D disc mesh from a loop of rings
    
    #Sim Parameters
    boxLength = 8.0            #length of computational domain (m)
    Nres = 512                  #NFINEST
    dX = boxLength/(1.0*Nres)   #grid mesh width (m)
    
    #Sphere Parameters
    listRadii = np.array([0.30,0.15])
    structureName = 'spheres'
    restSpringLength = 3.0*listRadii[0]
    discOrigin = [[-0.5*restSpringLength,0.0,0.0],[0.5*restSpringLength,0.0,0.0]]
    vertexFactor = 1.0
    #Vertices
    upper_vertices = []
    lower_vertices = []

    #dirName = os.getcwd()+'/../Structures/Nick3D/constraintIB/N'+str(Nres)+'/'
    dirName = os.getcwd()+'/'
    pathlib.Path(dirName).mkdir(parents=True, exist_ok=True)
    #Upper
    #Vertices
    nUpper,upper_vertices = genDiscMesh(listRadii[0],discOrigin[0],dX, dirName,vertexFactor)
    #Lower
    #Vertices
    nLower,lower_vertices = genDiscMesh(listRadii[1],discOrigin[1],dX, dirName,vertexFactor)
    #Combine Vertices (plot and write)
    all_vertices = upper_vertices+lower_vertices
    print(len(all_vertices))
    PlotVertexData(dirName,structureName,np.array(all_vertices),nUpper)
    WriteVertexData(dirName,structureName,np.array(all_vertices),nUpper,nLower)
    
    return
    
    '''#Grid Parameters
    boxLength = 2.0             #length of computational domain (m)
    Nres = int(sys.argv[1])                  #NFINEST
    dX = boxLength/(1.0*Nres)   #grid mesh width (m)
    Reynolds = [2.0]
    
    #Disc Parameters
    listRadii = np.array([0.3,0.15])
    structureNames = ['large','small']
    restSpringLength = 6.0*listRadii[1]
    discOrigin = [[-0.5*restSpringLength,0.0,0.0],[0.5*restSpringLength,0.0,0.0]]
    vertexFactor = float(sys.argv[2])

    for Re in Reynolds:
        dirName = '../Structures/Nick3D/constraintIB/N'+str(Nres)+'/'
        pathlib.Path(dirName).mkdir(parents=True, exist_ok=True)
        for i in range(len(listRadii)):
            genDiscMesh(dX,listRadii[i],discOrigin[i],structureNames[i],dirName,vertexFactor)
    
    return'''
    
#-----------------END MAIN-----------------------
main()
