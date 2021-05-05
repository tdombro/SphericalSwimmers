#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:46:45 2019

@author: thomas
"""

import os, shutil
import pandas as pd
import numpy as np
import random
import zipfile
import matplotlib
import matplotlib.pyplot as plt
from shutil import copyfile
import pathlib

#Script will generate randomly placed circles in a 2D plane
#Script will then CheckCircleBounds
#If circles overlap, then they will move apart from one another in opposite directions
#Continue above (2,3) until no circle overlaps
#When circles collide, they will try to push off like a collision

cwd_PYTHON = os.getcwd()

#CONSTANTS GRID PLACEMENT
RADIUSLARGE = 0.002

#CONSTANTS SPHEROBOT
RList = [3.0,5.0,6.0,7.5]
structureNames = ['skeleton','botup','botlow']

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
    
def StoreVertexInfo():
    
    #BOTUP
    #Read vertex file line by line and save the values in a list using the \n delimiter
    linesup = [line.strip() for line in open('botup.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    linesup = [line.split() for line in linesup]
    nvertup = int(linesup[0][0])
    #Allocate Array for Large Sphere Vertex Positions
    vertUp = np.zeros((2,nvertup))
    #Store Vertices
    for i in range(1,nvertup+1):
        vertUp[0,i-1] = float(linesup[i][0])
        vertUp[1,i-1] = float(linesup[i][1])
    #BOTLOW
    #Read vertex file line by line and save the values in a list using the \n delimiter
    lineslow = [line.strip() for line in open('botlow.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    lineslow = [line.split() for line in lineslow]
    nvertlow = int(lineslow[0][0])
    #Allocate Array for Small Sphere Vertex Positions
    vertLow = np.zeros((2,nvertlow))
    #Store Vertices
    for i in range(1,nvertlow+1):
        vertLow[0,i-1] = float(lineslow[i][0])
        vertLow[1,i-1] = float(lineslow[i][1])
    #Read vertex file line by line and save the values in a list using the \n delimiter
    linesskel = [line.strip() for line in open('skeleton.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    linesskel = [line.split() for line in linesskel]
    nvertskel = int(linesskel[0][0])
    #Allocate Array for Skeleton Vertex Positions
    vertSkel = np.zeros((2,nvertskel))
    #Store Vertices
    for i in range(1,nvertskel+1):
        vertSkel[0,i-1] = float(linesskel[i][0])
        vertSkel[1,i-1] = float(linesskel[i][1])
    
    nvert = [nvertskel,nvertup,nvertlow]
    vertList = [vertSkel,vertUp,vertLow]
    
    return (vertList,nvert)

def DisplaceSpherobots(vertList, nvert, structureNames, R, Theta, idxT, idxConfig):
    
    #First Rotate based on idxConfig
    #Allocate Arrays
    rotationMatrix = np.zeros((2,2))
    if(idxConfig == 1):
        theta = np.pi
    elif(idxConfig == 2):
        theta = np.pi/2.0
    elif(idxConfig == 3):
        theta = -np.pi/2.0
    else:
        theta = 0.0
    #rotatedPosition = np.zeros((Nbots,2,nvert))
    
    #Generate Random Angles
    #print('theta[%i] = %.3e' %(i,theta[i]))
    rotationMatrix[0,0] = np.cos(theta)
    rotationMatrix[0,1] = -1.0*np.sin(theta)
    rotationMatrix[1,0] = np.sin(theta)
    rotationMatrix[1,1] = np.cos(theta)  
    
    #Displaces the spheres where they are a distance R apart at an angle Theta
    #x1 and x2
    x1Arr = np.zeros(2)
    x1Arr[0] = -0.5*R*RADIUSLARGE*np.cos(Theta)
    #print(np.sin(Theta))
    x2Arr = np.zeros(2)
    x2Arr[0], x2Arr[1] = 0.5*R*RADIUSLARGE*np.cos(Theta), R*RADIUSLARGE*np.sin(Theta)
    xList = [x1Arr, x2Arr]
    
    if(idxConfig == 0):
        pathlib.Path('../Structures/Periodic/PhiPI/'+str(int(R))+'/PI'+str(idxT)+'/Parallel/').mkdir(parents=True, exist_ok=True)
        cwd_PARALLEL = cwd_PYTHON + '/../Structures/Periodic/PhiPI/'+str(int(R))+'/PI'+str(idxT)+'/Parallel/'
    elif(idxConfig == 1):
        pathlib.Path('../Structures/Periodic/PhiPI/'+str(int(R))+'/PI'+str(idxT)+'/Anti/').mkdir(parents=True, exist_ok=True)
        cwd_ANTI = cwd_PYTHON + '/../Structures/Periodic/PhiPI/'+str(int(R))+'/PI'+str(idxT)+'/Anti/'
    elif(idxConfig == 2):
        pathlib.Path('../Structures/Periodic/PhiPI/'+str(int(R))+'/PI'+str(idxT)+'/PerpL/').mkdir(parents=True, exist_ok=True)
        cwd_PERPL = cwd_PYTHON + '/../Structures/Periodic/PhiPI/'+str(int(R))+'/PI'+str(idxT)+'/PerpL/'
    else:
        pathlib.Path('../Structures/Periodic/PhiPI/'+str(int(R))+'/PI'+str(idxT)+'/PerpS/').mkdir(parents=True, exist_ok=True)
        cwd_PERPS = cwd_PYTHON + '/../Structures/Periodic/PhiPI/'+str(int(R))+'/PI'+str(idxT)+'/PerpS/'
        
    #Generate Figure to show Pairwise Placement
    fig = plt.figure(num=0,figsize=(4,4),dpi=120)
    ax = fig.add_subplot(111)
    ax.set_title('Pairwise Initial Configuration: \nR = %.4f m Theta = PI*%.2f m'%(R*0.002,Theta/np.pi))
    ax.axis([-0.05,0.05,-0.05,0.05])
    
    #Displace Spherobots
    for idxBot in range(2):
        dispArr = xList[idxBot]
        #print(dispArr)
        for idxName in range(len(vertList)):
            name = structureNames[idxName]
            vertPos = vertList[idxName].copy()
            if(idxConfig == 0):
                f = open(cwd_PARALLEL+name+str(idxBot+1)+'.vertex','w')
                #Copy spring/beam files for 'name'
                copyfile(name+'.spring',cwd_PARALLEL+name+str(idxBot+1)+'.spring')
            elif(idxConfig == 1):
                f = open(cwd_ANTI+name+str(idxBot+1)+'.vertex','w')
                #Copy spring/beam files for a'name'
                copyfile(name+'.spring',cwd_ANTI+name+str(idxBot+1)+'.spring')
            elif(idxConfig == 2):
                f = open(cwd_PERPL+name+str(idxBot+1)+'.vertex','w')
                #Copy spring/beam files for a'name'
                copyfile(name+'.spring',cwd_PERPL+name+str(idxBot+1)+'.spring')
            else:
                f = open(cwd_PERPS+name+str(idxBot+1)+'.vertex','w')
                #Copy spring/beam files for a'name'
                copyfile(name+'.spring',cwd_PERPS+name+str(idxBot+1)+'.spring')
                
            f.write('%i\n'%nvert[idxName])

            for idxVert in range(nvert[idxName]):
                #Rotate Skeleton2 by Theta given idxConfig
                if(idxName == 0 and idxBot == 1):
                    #print('b4: idxVert = %i: xPos = %.5e: yPos = %.5e'%(idxVert,vertPos[0,idxVert],vertPos[1,idxVert]))
                    if(idxVert <= 12):
                        vertPos[:,idxVert] = rotationMatrix.dot(vertPos[:,idxVert])
                    else:
                        CM = np.array([0.0,-0.005])
                        vertPos[:,idxVert] = rotationMatrix.dot(vertPos[:,idxVert].copy() - CM)
                        vertPos[:,idxVert] += CM
                    #print('a4: idxVert = %i: xPos = %.5e: yPos = %.5e'%(idxVert,vertPos[0,idxVert],vertPos[1,idxVert]))
                #Displace Spherobot by xList[idxBot]
                vertPos[:,idxVert] += dispArr[:]
                #Flip Spherobot 2 if AntiParallel
                if(idxConfig == 1 and idxBot == 1):
                    if(idxName == 0):
                        #Skeleton
                        if(idxVert <= 12):
                            vertPos[1,idxVert] -= 0.005
                        else:
                            vertPos[1,idxVert] += 0.005
                    elif(idxName == 1):
                        #Large Sphere
                        vertPos[1,idxVert] -= 0.005
                    elif(idxName == 2):
                        #Small Sphere
                        vertPos[1,idxVert] += 0.005
                #Rotate 90 degrees if PerpL
                if(idxConfig == 2 and idxBot == 1):
                    if(idxName == 0):
                        #Skeleton
                        if(idxVert > 12):
                            vertPos[0,idxVert] += 0.005
                            vertPos[1,idxVert] += 0.005
                    elif(idxName == 2):
                        #Small Sphere
                        vertPos[0,idxVert] += 0.005
                        vertPos[1,idxVert] += 0.005
                #Flip Spherobot 2 and rotate 90 degrees if PerpS
                if(idxConfig == 3 and idxBot == 1):
                    if(idxName == 0):
                        #Skeleton
                        if(idxVert <= 12):
                            vertPos[0,idxVert] += 0.005
                        else:
                            vertPos[1,idxVert] += 0.005
                    elif(idxName == 1):
                        #Large Sphere
                        vertPos[0,idxVert] += 0.005
                    elif(idxName == 2):
                        #Small Sphere
                        vertPos[1,idxVert] += 0.005
                #Write vertex coordinates down in .vertex file
                if(idxVert == nvert[idxName] - 1):
                    f.write('%.5e %.5e'%(vertPos[0,idxVert],vertPos[1,idxVert]))
                else:
                    f.write('%.5e %.5e\n' %(vertPos[0,idxVert],vertPos[1,idxVert]))
                
                                
            f.close()
            
            #Plot Displaced Spherobots
            if(idxName == 0):
                #Skeleton
                ax.plot(vertPos[0,:],vertPos[1,:],'ro',zorder=5,markersize=2)
                ax.plot(vertPos[0,13],vertPos[1,13],'bo',zorder=6,markersize=2)
                ax.plot(vertPos[0,0],vertPos[1,0],'bo',zorder=6,markersize=2)
            else:
                #Large and Small Spheres
                ax.plot(vertPos[0,:],vertPos[1,:],'ko',zorder=1,markersize=2)
    
    fig.tight_layout()
    if(idxConfig == 0):
        #Parallel Configuration
        fig.savefig(cwd_PARALLEL+'InitConfig.png')
    elif(idxConfig == 1):
        #AntiParallel Configuration
        fig.savefig(cwd_ANTI+'InitConfig.png')
    elif(idxConfig == 2):
        #PerpL Configuration
        fig.savefig(cwd_PERPL+'InitConfig.png')
    else:
        #PerpS Configuration
        fig.savefig(cwd_PERPS+'InitConfig.png')
    fig.clf()
    plt.close()
    
    return
    

if __name__ == '__main__':
    #Generate Placement for Pairwise Configurations: Parallel and Anti-Parallel
    #1) Read in .vertex files
    #2) Store Vertices in array (vertexPos array)
    #3) Loop over R and Theta
    #4) Displace spherobot 1 by x1 and spherobot2 by x2
    #5) x1 = (-1.0*Lcos(theta),0.0); x2 = (Lcos(theta),Lsin(theta))
    #6) If Antiparallel, switch ind sphere locations: LS -= 0.005; SS += 0.005
    #7) Write new vertex positions to new .vertex file
    #8) copy .beam and .spring files over to same dir as new .vertex files
    #9) Zip newly created files
    
    #1)Read in .vertex files
    #2)Store Vertices in array
    vertList, nvert = StoreVertexInfo()
    nvertskel, nvertup, nvertlow = nvert[0], nvert[1], nvert[2]
    vertSkel, vertUp, vertLow = vertList[0], vertList[1], vertList[2]
    
    #3)Loop over R and Theta (Parallel, Anti-Parallel, PerpL, and PerpS)
    for idxConfig in range(0,4):
            #for idxR in range(0,3):
        for idxR in range(4):
            #R = 5.0 + 2.5*idxR
            #R = 5.0 + 1.0*idxR
            R = RList[idxR]
            for idxT in range(5):
                Theta = -1.0*np.pi/2.0 + idxT*np.pi/4.0
                #4) Displace spherobot by x1 and x2
                DisplaceSpherobots(vertList, nvert, structureNames, R, Theta, idxT,idxConfig)
        
