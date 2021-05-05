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
#from genSkeleton import SkeletonGenerator
import os, sys
import zipfile
from shutil import copyfile
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

def genDiscMesh(dX,definedRadius,origin,vertexFactor):
    #Mesh Parameters
    #dS = 2.0*dX                                      #CIB
    dS = vertexFactor*dX                                     #Standard IB
    minRadius = dS                                  #Minimum radius of disc (m)
    maxRadius = definedRadius - 1.25*dX                           #Maximum radius of disc (m)
    print('nRings = ',(maxRadius - minRadius)/dS+1)
    #Make new dX value for Ring separation
    nRings = (maxRadius - minRadius)/dS + 1       #Number of rings in disc
    print('b4: dX = ',dX)
    #dX = 2.0*(maxRadius-minRadius)/nRings
    dX *= nRings/int(np.ceil(nRings))
    print('a4: dX = ',dX)
    nRings = int(np.ceil(nRings))                       #Number of rings in disc
    print('nRings = ',nRings)
    dS = vertexFactor*dX
    totalPoints = 1                                     #Total number of points in mesh (only origin at first)
    npts = totalPoints
    #Obtain total number of points in current disc to create array of that size
    for i in range(nRings):
        #Ring Parameters
        radius = (i+1)*dS
        #radius = minRadius + i*0.5*dX                   #Radius of band (m)
        bandLength = 2.0*np.pi*radius                   #Circumference of band (m)
        npts = m.ceil(bandLength/dS)  #Number of vertex points along band
        print('i = %i\tnpts = %i' %(i,npts))
        totalPoints += npts
    print('totalPoints = ',totalPoints)
    
    #Vertex parameters
    arrDiscVertex = np.zeros((totalPoints,2))
    arrDiscVertex[0,0], arrDiscVertex[0,1] = origin[0], origin[1]   #First vertex point is at the origin
    #Record origin as first vertex coordinate
    indexNumber = 1
    
    for i in range(nRings):
        #Ring Parameters
        radius = (i+1)*dS
        #radius = minRadius + i*0.5*dX                   #Radius of band (m)
        bandLength = 2.0*np.pi*radius                   #Circumference of band (m)
        npts = m.ceil(bandLength/dS)              #Number of vertex points along band
        print('i = %i\tnpts = %i\tradius = %.5e' %(i,npts,radius))
        dtheta = 2.0*np.pi/npts                         #Change in angle between vertex points

        #Allocate Arrays
        arrVertex = np.zeros((npts,2))                  #Allocate array for vertex coord information

        for j in range(npts):
            #Vertex coordinate (including shift)
            arrVertex[j,0] = radius*np.cos(j*dtheta) + origin[0]
            arrVertex[j,1] = radius*np.sin(j*dtheta) + origin[1]
            arrDiscVertex[indexNumber,0] = arrVertex[j,0]
            arrDiscVertex[indexNumber,1] = arrVertex[j,1]
            indexNumber += 1
    
    #Generate the springs that hold the mesh together!
    #arrDiscSpring = genDiscSprings(arrDiscVertex)
    
    #return (arrDiscVertex, arrDiscSpring)
    return arrDiscVertex

def PlotVertexData(dirName,structureName,data):
    #Plot parameters
    fig1 = plt.figure(num=1,figsize=(5,5),dpi=120)
    ax1 = fig1.add_subplot(111)
    ax1.scatter(data[:,0], data[:,1] ,s=1,c='r')
    ax1.scatter(data[0,0],data[0,1],s=1,c='k',zorder=5)
    ax1.axis('equal')
    fig1.savefig(dirName+structureName+'.png',bbox_inches='tight',pad_inches=0.5)
    fig1.clf()
    #plt.show()
    plt.close()
    return

def WriteVertexData(dirName,structureName,data):
    nPoints = len(data[:,0])
    f = open(dirName+structureName+'.vertex','w')
    f.write('%i\n' % nPoints)
    for idx in range(nPoints):
        f.write('%.5e %.5e\n' %(data[idx,0],data[idx,1]))
    f.close()
    return

#Distance between 2 identified vertices
def getSpringLength(xpos,ypos):
    return np.sqrt(xpos*xpos + ypos*ypos)

def genDiscSprings(pointArr):
    #Generate Springs!
    #Kstiff = 5.0e4
    #springArr = np.zeros((100000,4))
    nSprings = 0

    #2D Delaunay Triangulation
    tri = Delaunay(pointArr)

    pairs = np.zeros((tri.simplices.size,2))
    pairCount = 0
    #Generate Pairs of vertices
    for simplex in tri.simplices:
        #Store Vertex pairs
        pairs[pairCount,:] = [simplex[0],simplex[1]]
        pairCount += 1
        pairs[pairCount,:] = [simplex[1],simplex[2]]
        pairCount += 1
        pairs[pairCount,:] = [simplex[0],simplex[2]]
        pairCount += 1
    
    print('pairCount = ',pairCount)
    
    savedPairs = np.zeros((int(pairCount),2))
    boolPair = 1
    sPairCount = 0
    
    print('Size of tri.simplices = ',tri.simplices.size)
    
    #Report all non-repeating pairs
    for i in range(pairCount):
        for j in range(i+1, pairCount):
            if(j != i):
                if((pairs[i,0] == pairs[j,0] and pairs[i,1] == pairs[j,1]) 
                or (pairs[i,0] == pairs[j,1] and pairs[i,1] == pairs[j,0])):
                    boolPair = 0
                    break
            else:
                boolPair = 0
                break
        if(boolPair == 1):
            #convert index value used to find dist b/w 2 points
            savedPairs[sPairCount,:] = [pairs[i,0],pairs[i,1]]
            #print("sPC[] = [%f, %f]\tsPC = %i" %(savedPairs[sPairCount,0], savedPairs[sPairCount,1],sPairCount))
            sPairCount += 1
        boolPair = 1
        
    sPx, sPy = savedPairs.shape
    print('Shape of savedPairs: %i, %i' %(sPx,sPy))
    
    springArr = np.zeros((sPairCount,3))
    #Generate spring file in a new file called disc_test.spring
    for i in range(sPairCount):
        #convert index value used to find vertex coordinates
        vert1 = [pointArr[int(savedPairs[i,0]),0], pointArr[int(savedPairs[i,0]),1]]
        vert2 = [pointArr[int(savedPairs[i,1]),0], pointArr[int(savedPairs[i,1]),1]]
        
        #Find distance b/w the 2 vertices
        dist = getSpringLength(vert1[0] - vert2[0],vert1[1] - vert2[1])
        springArr[nSprings,:] = [int(savedPairs[i,0]),int(savedPairs[i,1]),dist]
        nSprings += 1
    print('nSprings = ',nSprings)

    return springArr

def PlotSpringData(dirName,structureName,vData,sData):
    #Plot vertices and springs b/w each pair of vertices
    fig2 = plt.figure(num=2,figsize=(10,10),dpi=120)
    ax2 = fig2.add_subplot(111)
    nSprings = len(sData[:,0])
    for i in range(nSprings):
        ax2.plot([vData[int(sData[i,0]),0],vData[int(sData[i,1]),0]], [vData[int(sData[i,0]),1],vData[int(sData[i,1]),1]], color = 'g')
    ax2.axis('equal')
    fig2.savefig(dirName+structureName+'Springs.png')
    fig2.clf()
    plt.close()
    return

def WriteSpringData(dirName,structureName,data,Kstiff):
    f = open(dirName+structureName+'.spring','w')
    #All springs have been generated
    nSprings = len(data[:,0])
    f.write('%i\n' % nSprings)
    #Store springArr in .spring file
    for i in range(nSprings):
        f.write('%i %i %.5e %.5e\n' %(int(data[i,0]),int(data[i,1]),Kstiff,data[i,2]))
    f.close()
    return
    
def main():
    #Goal: Generate a uniform 2D disc mesh from a loop of rings
    
    #Grid Parameters
    boxLength = 16.0            #length of computational domain (m)
    Nres = 400                  #NFINEST
    dX = boxLength/(1.0*Nres)   #grid mesh width (m)
    vertexFactor = 0.5 #StandardIB
    
    #Disc Parameters
    listRadii = 0.5
    structureNames = 'cylinder'
    discOrigin = [0.0,0.0]

    #allData = [[[] for a in range(2)] for b in range(2)] #a = #rows = vertex,spring, b = #cols = #structures
    dirName = os.getcwd()+'/'
    #pathlib.Path(dirName).mkdir(parents=True, exist_ok=True)
    vertexData = genDiscMesh(dX,listRadii,discOrigin,vertexFactor)
    PlotVertexData(dirName,structureNames,vertexData)
    WriteVertexData(dirName,structureNames,vertexData)
    
    '''for kdx in range(2,11,2):
        K_osc = 1.0e4*kdx
        K_body = 5.0*K_osc
        dirName = os.getcwd()+'/../Structures/N512/K%.1f/'%(round(K_osc/1.0e4,1))
        pathlib.Path(dirName).mkdir(parents=True, exist_ok=True)
        for i in range(len(listRadii)):
            PlotVertexData(dirName,structureNames[i],allData[i][0])
            WriteVertexData(dirName,structureNames[i],allData[i][0])
            PlotSpringData(dirName,structureNames[i],allData[i][0],allData[i][1])
            WriteSpringData(dirName,structureNames[i],allData[i][1],K_body)
            #genDiscMesh(dX,listRadii[i],discOrigin[i],structureNames[i],dirName,K_body)
        SkeletonGenerator(listRadii,restSpringLength,dX,dirName,K_body,K_osc)
    '''
    
    return
    
#-----------------END MAIN-----------------------
main()
