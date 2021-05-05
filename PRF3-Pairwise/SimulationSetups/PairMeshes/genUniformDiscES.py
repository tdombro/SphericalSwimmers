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
from genSkeletonES import SkeletonGenerator
import os
import zipfile
from shutil import copyfile

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

def genDiscMesh(dX,definedRadius,origin,structureName,dirName):
    #FILES
    file = open(dirName+structureName+'.vertex','w')
    
    #Plot parameters
    fig1 = plt.figure(num=1,figsize=(5,5),dpi=120)
    ax1 = fig1.add_subplot(111)
    
    #Mesh Parameters
    minRadius = 0.5*dX                                  #Minimum radius of disc (m)
    maxRadius = definedRadius - 1.25*dX                 #Maximum radius of disc (m)
    print((maxRadius - minRadius)/(0.5*dX)+1)
    #Make new dX value for Ring separation
    nRings = (maxRadius - minRadius)/(0.5*dX) + 1       #Number of rings in disc
    print('b4: dX = ',dX)
    #dX = 2.0*(maxRadius-minRadius)/nRings
    dX *= nRings/int(np.ceil(nRings))
    print('a4: dX = ',dX)
    nRings = int(np.ceil(nRings))                       #Number of rings in disc
    print('nRings = ',nRings)
    totalPoints = 1                                     #Total number of points in mesh (only origin at first)
    npts = totalPoints
    #Obtain total number of points in current disc to create array of that size
    for i in range(nRings):
        #Ring Parameters
        radius = (i+1)*0.5*dX
        #radius = minRadius + i*0.5*dX                   #Radius of band (m)
        bandLength = 2.0*np.pi*radius                   #Circumference of band (m)
        npts = m.ceil((2.0*bandLength/dX))  #Number of vertex points along band
        print('i = %i\tnpts = %i' %(i,npts))
        totalPoints += npts
    print('totalPoints = ',totalPoints)
    #Record total number of vertices
    file.write('%i\n' % totalPoints)
    
    #Vertex parameters
    arrDiscVertex = np.zeros((totalPoints,2))
    arrDiscVertex[0,0], arrDiscVertex[0,1] = origin[0], origin[1]   #First vertex point is at the origin
    file.write('%.5e %.5e\n' %(arrDiscVertex[0,0],arrDiscVertex[0,1]))
    #Record origin as first vertex coordinate
    ax1.scatter(origin[0],origin[1],s=1,c='r')
    indexNumber = 1
    
    for i in range(nRings):
        #Ring Parameters
        radius = (i+1)*0.5*dX
        #radius = minRadius + i*0.5*dX                   #Radius of band (m)
        bandLength = 2.0*np.pi*radius                   #Circumference of band (m)
        npts = m.ceil((2.0*bandLength/dX))              #Number of vertex points along band
        print('i = %i\tnpts = %i\tradius = %.5e' %(i,npts,radius))
        ds = bandLength/npts                            #distance between vertex points
        dtheta = 2.0*np.pi/npts                         #Change in angle between vertex points

        #Allocate Arrays
        arrVertex = np.zeros((npts,2))                  #Allocate array for vertex coord information

        for j in range(npts):
            #Vertex coordinate (including shift)
            arrVertex[j,0] = radius*np.cos(j*dtheta) + origin[0]
            arrVertex[j,1] = radius*np.sin(j*dtheta) + origin[1]
            arrDiscVertex[indexNumber,0] = arrVertex[j,0]
            arrDiscVertex[indexNumber,1] = arrVertex[j,1] 
            file.write('%.5e %.5e\n' %(arrDiscVertex[indexNumber,0],arrDiscVertex[indexNumber,1]))
            indexNumber += 1

        ax1.scatter(arrVertex[:,0], arrVertex[:,1] ,s=1,c='r')
    
    #Close text file
    file.close()
    #even axes and save fig
    ax1.axis('equal')
    fig1.savefig(dirName+structureName+'.png',bbox_inches='tight',pad_inches=0.5)
    fig1.clf()
    #plt.show()
    plt.close()
    
    #Generate the springs that hold the mesh together!
    #genDiscSprings(arrDiscVertex,structureName,dirName)
    
    return

#Distance between 2 identified vertices
def getSpringLength(xpos,ypos):
    return np.sqrt(xpos*xpos + ypos*ypos)

def genDiscSprings(pointArr,structureName,dirName):    
    #Generate Springs!
    Kstiff = 5.0e4
    springArr = np.zeros((100000,4))
    nSprings = 0
    #FILES
    f = open(dirName+structureName+'.spring','w')
    
    #Plot vertices and springs b/w each pair of vertices
    fig2 = plt.figure(num=2,figsize=(10,10),dpi=120)
    ax2 = fig2.add_subplot(111)

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
    
    #Generate spring file in a new file called disc_test.spring
    for i in range(sPairCount):
        #convert index value used to find vertex coordinates
        vert1 = [pointArr[int(savedPairs[i,0]),0], pointArr[int(savedPairs[i,0]),1]]
        vert2 = [pointArr[int(savedPairs[i,1]),0], pointArr[int(savedPairs[i,1]),1]]
        
        #Find distance b/w the 2 vertices
        dist = getSpringLength(vert1[0] - vert2[0],vert1[1] - vert2[1])
        springArr[nSprings,:] = [int(savedPairs[i,0]),int(savedPairs[i,1]),Kstiff,dist]
        nSprings += 1
    print('nSprings = ',nSprings)
        
    #All springs have been generated
    f.write('%i\n' % nSprings)
    #Store springArr in .spring file
    for i in range(nSprings):
        f.write('%i %i %.5e %.5e\n' %(int(springArr[i,0]),int(springArr[i,1]),springArr[i,2],springArr[i,3]))
        ax2.plot([pointArr[int(springArr[i,0]),0],pointArr[int(springArr[i,1]),0]], [pointArr[int(springArr[i,0]),1],pointArr[int(springArr[i,1]),1]], color = 'g')     
    f.close()
    #ax1.axis([-0.40,-0.30,-0.05,0.05])
    ax2.axis('equal')   
    fig2.savefig(dirName+structureName+'Springs.png')
    fig2.clf()
    plt.close()
    return
    
def main():
    #Goal: Generate a uniform 2D disc mesh from a loop of rings
    
    #Grid Parameters
    boxLength = 0.05             #length of computational domain (m)
    Nres = 512                  #NFINEST
    dX = boxLength/(1.0*Nres)   #grid mesh width (m)
    Reynolds = [2.0,80.0]
    
    #Disc Parameters
    listRadii = np.array([0.0015,0.0015])
    structureNames = ['botupES','botlowES']
    restSpringLength = 1.0*(np.sum(listRadii) + 0.002)
    discOrigin = [[0.0,0.0],[0.0,-1.0*restSpringLength]]

    dirName = os.getcwd()+'/'
    for i in range(len(listRadii)):
        genDiscMesh(dX,listRadii[i],discOrigin[i],structureNames[i],dirName)
    SkeletonGenerator(listRadii,restSpringLength,dX,dirName)

    '''for Re in Reynolds:
        MakeDirectory('../Structures','Re'+str(Re))
        dirName = '../Structures/Re'+str(Re)+'/'
        for i in range(len(listRadii)):
            genDiscMesh(dX,listRadii[i],discOrigin[i],structureNames[i],dirName)
        SkeletonGenerator(listRadii,restSpringLength,dX,dirName)'''
        
        #Copy files needed for ibamr
        #copyfile('1botSIB.bsub',dirName+'1botSIB.bsub')
        #copyfile('input2d',dirName+'input2d')
        #copyfile('main.C',dirName+'main.C')
        #copyfile('Makefile',dirName+'Makefile')
        #copyfile('update_springs.C',dirName+'update_springs.C')
        #copyfile('update_springs.h',dirName+'update_springs.h')
    
    #zipFiles('../Structures/','UniformMesh')
    
    return
    
#-----------------END MAIN-----------------------
main()
