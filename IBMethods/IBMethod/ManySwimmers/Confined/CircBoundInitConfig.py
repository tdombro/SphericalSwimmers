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

#Script will generate randomly placed circles in a 2D Circle
#Script will then CheckCircleBounds
#If circles overlap, then they will move apart from one another in opposite directions
#Continue above (2,3) until no circle overlaps
#When circles collide, they will try to push off like a collision

#CONSTANTS GRID PLACEMENT
Radius = 0.005625
UpperBoundR = Radius*8.0
LowerBoundX, UpperBoundX = 0.0, Radius*16.0
LowerBoundY, UpperBoundY = 0.0, Radius*16.0
FIGNUM = 0

#CONSTANTS SPHEROBOT
structureNames = ['skeleton','botup','botlow']
seed = 27510

#Given Area Fraction, Calculate Ncircles
AreaFraction = 0.6
AreaCircle = np.pi*UpperBoundR**2
Ncircles = np.int(AreaFraction*AreaCircle/(np.pi*Radius**2))
#Ncircles = 9
print('Ncircles = ',Ncircles)
AreaNumCircle = np.pi*Radius**2*Ncircles
AreaFractionActual = AreaNumCircle/AreaCircle
print('AreaFractionActual = ',AreaFractionActual)

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
    #BOTLOW
    #Read vertex file line by line and save the values in a list using the \n delimiter
    lineslow = [line.strip() for line in open('botlow.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    lineslow = [line.split() for line in lineslow]
    nvertlow = int(lineslow[0][0])
    #Read vertex file line by line and save the values in a list using the \n delimiter
    linesskel = [line.strip() for line in open('skeleton.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    linesskel = [line.split() for line in linesskel]
    nvertskel = int(linesskel[0][0])
    
    nvert = [nvertskel,nvertup,nvertlow]
    
    return (linesup,lineslow,nvert)
    
def FindCM(linesup,lineslow,nvertup,nvertlow):
    #Determine COM for structure
    CM = np.zeros((Ncircles,2))
    #Array of combined structure positions
    combinedPosition = np.zeros((2,nvertup+nvertlow))
    for i in range(1,nvertup+1):
        combinedPosition[0,i-1] = float(linesup[i][0])
        combinedPosition[1,i-1] = float(linesup[i][1])
    print('Xpos = ',combinedPosition[0,nvertup-1])
    print('Ypos = ',combinedPosition[1,nvertup-1])
    for i in range(nvertup+1,nvertup+nvertlow+1):
        combinedPosition[0,i-1] = float(lineslow[i-nvertup][0])
        combinedPosition[1,i-1] = float(lineslow[i-nvertup][1])
        if(i == nvertup+1):
            print('Xpos = ',combinedPosition[0,i-1])
            print('Ypos = ',combinedPosition[1,i-1])
    
    #Find the COM of our structure
    CM[:,0] = np.sum(combinedPosition[0,:])/float(nvertup+nvertlow)
    print('COMx = ',CM[0,0])
    CM[:,1] = np.sum(combinedPosition[1,:])/float(nvertup+nvertlow)
    print('COMy = ',CM[0,1])
    
    return CM

def GenerateRandomCenters(center,Ncircles):
    #global LowerBoundX, UpperBoundX, LowerBoundY, UpperBoundY
    global UpperBoundR
    #In this function, we will determine the position of the centers of our circles
    randomR = np.random.uniform(low=0.0,high=UpperBoundR-Radius,size=Ncircles)
    randomTheta = np.random.uniform(low=0.0,high=2.0*np.pi,size=Ncircles)
    center[0,:] = randomR*np.cos(randomTheta)
    center[1,:] = randomR*np.sin(randomTheta)
    #center[0,:] = np.random.uniform(low=LowerBoundX+Radius,high=UpperBoundX-Radius,size=Ncircles)
    #center[1,:] = np.random.uniform(low=LowerBoundY+Radius,high=UpperBoundY-Radius,size=Ncircles)
    #print(center)
    return
    
def CalculateDistanceBwCenters(circle1,circle2):
    #This function calculates the distance between specified circle centers
    distance = np.sqrt((circle1[0] - circle2[0])**2 + (circle1[1] - circle2[1])**2)
    return distance

def CalculateAngleBwCenters(circle1,circle2,dist):
    #This function calculates the angle between specified circle centers
    #It will be used to prevent overlap
    Xdist = circle1[0] - circle2[0]
    Angle = np.arccos(Xdist/dist)
    return Angle

def CheckCircleBounds(circle):
    global LowerBoundX, UpperBoundX, LowerBoundY, UpperBoundY
    #This function checks to see if the circle is within the bounds of the plane
    circleRadius = np.sqrt(circle[0]*circle[0] + circle[1]*circle[1])
    circleAngle = np.tan(circle[1]/circle[0])
    if(circleRadius + Radius > UpperBoundR):
        circle[0] = (UpperBoundR - 1.5*Radius)*np.cos(circleAngle)
        circle[1] = (UpperBoundR - 1.5*Radius)*np.sin(circleAngle)
    '''if(circle[0] - Radius < LowerBoundX):
        circle[0] = LowerBoundX + 1.5*Radius
    elif(circle[0] + Radius > UpperBoundX):
        circle[0] = UpperBoundX - 1.5*Radius
    if(circle[1] - Radius < LowerBoundY):
        circle[1] = LowerBoundY + 1.5*Radius
    elif(circle[1] + Radius > UpperBoundY):
        circle[1] = UpperBoundY - 1.5*Radius'''
    return circle

def PlotCircleLocations(center,Ncircles):
    global FIGNUM
    #Plot Random Centers on 2D Plane
    circles = [0]*Ncircles
    circlesBorder = [0]*Ncircles
    #Plot Initialized Random Circle Locations
    fig = plt.figure(num=FIGNUM,figsize=(8,8),dpi=500)
    ax = fig.add_subplot(111)
    
    ax.set_title('Final Configuration: N = %i R = %.3f'%(Ncircles,Radius))
    #ax.scatter(center[0],center[1],s=4)
    #Plot Circles generated from centers
    for idx in range(Ncircles):
        circles[idx] = plt.Circle((center[0,idx],center[1,idx]),Radius,color='r',clip_on=False)
        circlesBorder[idx] = plt.Circle((center[0,idx],center[1,idx]),Radius,color='k',clip_on=False,fill=False)
        ax.add_artist(circles[idx])
        ax.add_artist(circlesBorder[idx])
    BoundaryCircle = plt.Circle((0.0,0.0),UpperBoundR,color='k',fill=False,clip_on=False)
    ax.add_artist(BoundaryCircle)
    #ax.axis([LowerBoundX,UpperBoundX,LowerBoundY,UpperBoundY])
    ax.axis([-UpperBoundR,UpperBoundR,-UpperBoundR,UpperBoundR])
    fig.tight_layout()
    fig.savefig('FinalConfig.png')
    #fig.clf()
    #plt.close()
    FIGNUM += 1
    return (fig,ax)
    
def GenerateInitialPositions(ax,center,Ncircles):
    #Determine if any of the circles are overlapping
    #Calculate the distance between each center
    #If distance < 2*Radius, then return true for overlap
    
    #Allocate Array of overlapping indices
    #The array will be 1 dimensional. each 2 elements is a pair that overlap
    overlap = []
    boolOverlap = True
    
    while(boolOverlap == True):
        #print('boolOverlap = ',boolOverlap)
        #print('Checking Again!')
        count = 0
        for idx1 in range(Ncircles):
            for idx2 in range(idx1+1,Ncircles):
                dist = CalculateDistanceBwCenters(center[:,idx1],center[:,idx2])
                if(dist < 2.0*Radius):
                    count += 1
                    boolOverlap = True
                    #print('Circles %i and %i overlap!'%(idx1,idx2))
                    #overlap = np.append(overlap,[idx1,idx2,dist])
                    overlap = [idx1,idx2,dist]
                    #Center 1 is higher than center 2
                    if(center[1,int(overlap[0])] > center[1,int(overlap[1])]):
                        i1 = int(overlap[0])
                        i2 = int(overlap[1])
                    #Center 2 is higher than center 1
                    else:
                        i1 = int(overlap[1])
                        i2 = int(overlap[0])
                    #Angle Between Centers
                    Angle = CalculateAngleBwCenters(center[:,i1],center[:,i2],dist)
                    #Move centers of overlapping circles
                    #Center 1 needs to move opposite the angle
                    #(2.0*Radius - dist)
                    center[0,i1] -= 0.5*(2.1*Radius - dist)*np.cos(np.pi-Angle)
                    center[1,i1] += 0.5*(2.1*Radius - dist)*np.sin(np.pi-Angle)
                    center[0,i2] -= 0.5*(2.1*Radius - dist)*np.cos(Angle)
                    center[1,i2] -= 0.5*(2.1*Radius - dist)*np.sin(Angle)
                    #Check to see if circle is with boundaries again
                    center[:,idx1] = CheckCircleBounds(center[:,idx1])
                    center[:,idx2] = CheckCircleBounds(center[:,idx2])
                    #Plot circles that overlapped on initial plot to see changes
                    circleOverlap1 = plt.Circle((center[0,i1],center[1,i1]),0.2*Radius,color='b',clip_on=False)
                    ax.add_artist(circleOverlap1)
                    circleOverlap2 = plt.Circle((center[0,i2],center[1,i2]),0.2*Radius,color='g',clip_on=False)
                    ax.add_artist(circleOverlap2)
                    
                
        if(count == 0):
            print('No Overlap!')
            boolOverlap = False
    
    return center

def RotateVertices(Nbots,nvert,structureName,gridPosition,CM,center,ax1):  
    #Rotates the set up grid by a randomly generated angle about the COM
    
    #Allocate Arrays
    rotationMatrix = np.zeros((Nbots,2,2))
    theta = np.zeros(Nbots)
    rotatedPosition = np.zeros((Nbots,2,nvert))
    #COM = np.zeros((Nbots,2))
    
    #Generate Random Angles
    np.random.seed(seed)
    for i in range(Nbots):
        theta[i] = 2.0*np.pi*np.random.rand(1)
        #print('theta[%i] = %.3e' %(i,theta[i]))
        rotationMatrix[i,0,0] = np.cos(theta[i])
        rotationMatrix[i,0,1] = -1.0*np.sin(theta[i])
        rotationMatrix[i,1,0] = np.sin(theta[i])
        rotationMatrix[i,1,1] = np.cos(theta[i])     
        
    #Create Directory of random seed
    MakeDirectory('RotatedDisc',str(Nbots))
    MakeDirectory('RotatedDisc/'+str(Nbots),seed)
    #MakeDirectory('../'+str(Nbots)+'bot',str(seed))

    #Rotate Grid Structures around COM based off of above angle 
    for i in range(int(Nbots)):
        #open newly rotated vertex file
        f = open('RotatedDisc/'+str(Nbots)+'/'+str(seed)+'/'+structureName+str(i+1)+'.vertex','w')
        f.write('%i\n'%nvert)
        
        #Rotate coord values about center of mass
        for j in range(nvert):
            rotatedPosition[i,:,j] = rotationMatrix[i].dot(gridPosition[i,:,j] - CM[i,:])
            rotatedPosition[i,:,j] += center[:,i]
            if(j == nvert - 1):
                f.write('%.5e %.5e'%(rotatedPosition[i,0,j],rotatedPosition[i,1,j]))
            else:
                f.write('%.5e %.5e\n' %(rotatedPosition[i,0,j],rotatedPosition[i,1,j]))
        f.close()
        
        #Plot Rotated Skeleton
        ax1.plot(rotatedPosition[i,0,:],rotatedPosition[i,1,:],'bo')
    
    return

if __name__ == '__main__':
    #1)Given Area Fraction of circles, calculate Ncircles (Number of swimmers)
    #2)Read in .vertex files
    #3)Calculate CoM from .vertex files
    #4)Store vertices in array (vertexPos array)
    #5)Calculate Random CoM Positions (center array)
    #6)Rotate Randomly each swimmer
    #7)Displace CoM by center array (Add center array value to each vertexPos)
    #8)Write new vertex positions to new .vertex file
    #9) copy .beam and .spring files over to same dir as new .vertex files
    #10) Zip newly created files
    
    #1) Done above all code where CONSTANTS are given
    #2)Read in .vertex files
    linesup, lineslow, nvert = StoreVertexInfo()
    nvertup, nvertlow = nvert[1], nvert[2]
    
    #3) Calculate CM from .vertex files
    CM = FindCM(linesup,lineslow,nvertup,nvertlow)
    
    '''INITIAL PLACEMENT'''
    #5)Calculate Random CM Positions (center array)
    #Allocate Array for circle centers
    center = np.zeros((2,Ncircles))
    
    #Generate Random Centers
    GenerateRandomCenters(center,Ncircles)
    
    #Plot Random Centers on 2D Plane
    circles = [0]*Ncircles
    circlesBorder = [0]*Ncircles
    
    #Plot Initialized Random Circle Locations
    fig0 = plt.figure(num=FIGNUM,figsize=(4,4),dpi=120)
    ax0 = fig0.add_subplot(111)
    ax0.set_title('Initial Setup With Changes Overlapped\nN = %i R = %.1f'%(Ncircles,Radius))
    #Plot Circles generated from centers
    for idx in range(Ncircles):
        circles[idx] = plt.Circle((center[0,idx],center[1,idx]),Radius,color='r',clip_on=False)
        circlesBorder[idx] = plt.Circle((center[0,idx],center[1,idx]),Radius,color='k',clip_on=False,fill=False)
        ax0.add_artist(circles[idx])
        ax0.add_artist(circlesBorder[idx])
    BoundaryCircle = plt.Circle((0.0,0.0),UpperBoundR,color='k',fill=False,clip_on=False)
    ax0.add_artist(BoundaryCircle)                       
    #ax.scatter(center[0],center[1],s=4)
    #ax0.axis([LowerBoundX,UpperBoundX,LowerBoundY,UpperBoundY])
    ax0.axis([-UpperBoundR,UpperBoundR,-UpperBoundR,UpperBoundR])
    fig0.tight_layout()
    fig0.savefig('InitConfig.png')
    FIGNUM += 1
    
    #Determine if any of the circles are overlapping
    #Calculate the distance between each center
    #If distance < 2*Radius, then return true for overlap
    GenerateInitialPositions(ax0,center,Ncircles)
    
    fig0.savefig('ConfigChanges.png')
    fig0.clf()
    plt.close()
    fig1, ax1 = PlotCircleLocations(center,Ncircles)
    print(center)
    
    
    '''ROTATION'''
    #4)Store vertices in array (vertexPos array)
    for i in range(len(structureNames)):
        vertexPos = np.zeros((Ncircles,2,nvert[i]))
        #Read vertex file line by line and save the values in a list using the \n delimiter
        lines = [line.strip() for line in open(structureNames[i]+'.vertex')]
        #Break each list element into an array of numbers using the space delimiter
        lines = [line.split() for line in lines]
        for j in range(Ncircles):
            for k in range(1,nvert[i]+1):
                vertexPos[j,0,k-1] = float(lines[k][0])
                vertexPos[j,1,k-1] = float(lines[k][1])
        #6)Rotate Randomly each swimmer
        #7)Displace CoM by center array (Add center array value to each vertexPos)
        #8)Write new vertex positions to new .vertex file
        RotateVertices(Ncircles,nvert[i],structureNames[i],vertexPos,CM,center, ax1)
    
    ax1.plot(center[0,:],center[1,:],'ko')
    #ax1.axis([0.0,UpperBoundX,0.0,UpperBoundY])
    ax1.axis([-UpperBoundR,UpperBoundR,-UpperBoundR,UpperBoundR])
    #fig1.savefig('Rotated'+str(seed)+'.png')
    fig1.savefig('InitialConfiguration.png')
    fig1.clf()
    plt.close()
    #9) copy .beam and .spring files over to same dir as new .vertex files
    for i in range(len(structureNames)):
        for j in range(1,Ncircles+1):
            copyfile(structureNames[i]+'.spring','RotatedDisc/'+str(Ncircles)+'/'+str(seed)+'/'+structureNames[i]+str(j)+'.spring')
            '''if(structureNames[i] == 'skeleton'):
                copyfile(structureNames[i]+'.beam','RotatedDisc/'+str(Ncircles)+'/'+str(seed)+'/'+structureNames[i]+str(j)+'.beam')'''
    #10) Zip newly created files
        
