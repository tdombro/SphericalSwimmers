#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:05:05 2018

@author: thomas
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

#Distance from center of Structure
def getVertexDistance(mSL,i,xpos,ypos):
    return np.sqrt(xpos*xpos + (-mSL*i - ypos)*(-mSL*i -ypos))

#Distance between 2 identified vertices
def getSpringLength(xpos,ypos):
    return np.sqrt(xpos*xpos + ypos*ypos)

def SkeletonGenerator(radius, restSpringLength, dX, dirName):
    
    #Constants
    Ndiscs = 2 #Number of discs in swimmer
    #restSpringLength = 0.375 #Longest spring length
    radius -= 1.25*dX                        #Radius of discs
    dS = min(radius)/2.0                     #min spacing b/w vertices in skeleton
    print('dS = ',dS)
    #VERTEX
    epsilon = 1.0e-5 #Correction factor for identifying vertices
    nPoints = 0 #number of vertices
    vertexArr = np.zeros((100,2)) #vertex Array
    #SPRINGS
    nSprings = 0 #Number of springs
    Kosc = 1.0e4 #Spring constant for oscillating spring
    Kstiff = 5.0e4 #Spring constant for stiff springs
    springArr = np.zeros((200,4)) #spring Array
    #BEAMS
    nBeams = 0 #number of beams
    Kbeam = 5.0e1 #Elastic spring constant
    beamArr = np.zeros((100,4)) #beam Array
    
    
    #FILES
    f = open(dirName+'skeletonES.vertex','w')
    #Generate Vertices in Loop
    #Will be in order: disc# xpos ascending ypos ascending
    for i in range(Ndiscs):
        #minimum x-coord and y-coord
        minX = -1.0*radius[i]
        minY = -1.0*(restSpringLength*i + radius[i])
        print('maxIndex = ',2*int(radius[i]/dS)+1)
        for j in range(5):
            #x position
            xpos = minX + 0.5*radius[i]*j
            #xpos = minX + dS*j
            if(j == 0):
                print('xpos = ',xpos)
            for k in range(5):
                #y position
                ypos = minY + 0.5*radius[i]*k
                #ypos = minY + dS*k
                if(j == 1 and k ==1):
                    print('xpos = ',xpos)
                    print('ypos = ',ypos)
                    print('radius = ',getVertexDistance(restSpringLength,i,xpos,ypos))
                    
                #Distance from center
                dist = getVertexDistance(restSpringLength,i,xpos,ypos)
                #Exclude points outside circle of radius r
                if(dist <= radius[i] + epsilon):
                    vertexArr[nPoints,0] = xpos
                    vertexArr[nPoints,1] = ypos
                    #print(vertexArr[nPoints])
                    nPoints += 1
    
    print('nPoints = ',nPoints)
    plt.figure(1)
    #Record Total number of vertices
    f.write('%i\n' % nPoints)
    #Loop over total number of vertex points
    for i in range(nPoints):
        #Plot Skeleton Vertices
        plt.plot(vertexArr[i,0],vertexArr[i,1],'ro')
        f.write('%.5e %.5e\n' %(vertexArr[i,0],vertexArr[i,1]))
        if(i == nPoints - 1):
            f.write('%.5e %.5e' % (vertexArr[i,0],vertexArr[i,1]))
    
    f.close()
    plt.axis([-1.0,1.0,-0.5,0.5])
    plt.axis('equal')
    plt.savefig(dirName+'skeleton.png')
    plt.gcf().clear()
    plt.close()
    
    #Generate Springs!
    #First record the long oscillating springs
    #Indices of interest: [4,6,8,17,19,21,30,32,34]
    #FILES
    f = open(dirName+'skeletonES.spring','w')
    #Going Left to Right, Down to Up, Diagonals last
    oscillatingIndexArr = [[2,13],[6,19],[10,25],[2,25],[10,13]]
    #oscillatingIndexArr = [[4,17],[6,19],[8,21],[4,21],[8,17],
                           #[17,30],[19,32],[21,34],[17,34],[21,30]]
    for index in oscillatingIndexArr:
        xpos = vertexArr[index[1],0] - vertexArr[index[0],0]
        ypos = vertexArr[index[1],1] - vertexArr[index[0],1]
        springLength = getSpringLength(xpos,ypos)
        springArr[nSprings,:] = [index[0],index[1],Kosc,springLength]
        #f.write('%i %i %.5e\n' %(index[0],index[1],springLength))
        nSprings += 1
    print('nSprings = ',nSprings)
    print(springArr[0:nSprings])
    #Make a new Array that includes all vertex points
    pointArr = np.zeros((nPoints,2))
    for i in range(nPoints):
        pointArr[i] = np.copy(vertexArr[i])
    
    #Plot vertices and springs b/w each pair of vertices
    fig3 = plt.figure(num=3,figsize=(10,10),dpi=120)
    ax3 = fig3.add_subplot(111)
    
    nVert = int(nPoints/Ndiscs)
    for k in range(Ndiscs):
        maxIndex = nVert*(k+1)
        print('maxIndex = ',maxIndex)
        #2D Delaunay Triangulation
        tri = Delaunay(pointArr[nVert*k:maxIndex,:])
        print(pointArr[nVert*k:maxIndex,:])
        xPA, yPA = pointArr[nVert*k:maxIndex,:].shape
        print('shape of pointArr: %i %i' %(xPA,yPA))
    
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
            vert1 = [pointArr[int(savedPairs[i,0])+nVert*k,0], pointArr[int(savedPairs[i,0])+nVert*k,1]]
            vert2 = [pointArr[int(savedPairs[i,1])+nVert*k,0], pointArr[int(savedPairs[i,1])+nVert*k,1]]
            
            #Find distance b/w the 2 vertices
            dist = getSpringLength(vert1[0] - vert2[0],vert1[1] - vert2[1])
            springArr[nSprings,:] = [int(savedPairs[i,0])+nVert*k,int(savedPairs[i,1])+nVert*k,Kstiff,dist]
            nSprings += 1
        print(k)
        print('nSprings = ',nSprings)
        
    #All springs have been generated
    f.write('%i\n' % nSprings)
    #Store springArr in .spring file
    for i in range(nSprings):
        f.write('%i %i %.5e %.5e\n' %(int(springArr[i,0]),int(springArr[i,1]),springArr[i,2],springArr[i,3]))
        ax3.plot([pointArr[int(springArr[i,0]),0],pointArr[int(springArr[i,1]),0]], [pointArr[int(springArr[i,0]),1],pointArr[int(springArr[i,1]),1]], color = 'g')     
    f.close()
    #ax1.axis([-0.40,-0.30,-0.05,0.05])
    ax3.axis('equal')   
    fig3.savefig(dirName+'skeletonSprings.png')
    fig3.clf()
    plt.close()
    
    
    #Generate Beams!
    f = open(dirName+'skeletonES.beam','w')
    #List of beams for first disc (One on the left)
        #[[0,1,4],[0,2,6],[0,3,8],
        # [1,2,3],[1,5,9],[1,6,11],[2,6,10],[3,6,9],[3,7,11]
        # [4,5,6],[4,9,12],[5,6,7],[6,7,8],[6,10,12],[8,11,12],
        # [9,10,11]]
    beamList = [[0,1,4],[0,2,6],[0,3,8],
                [1,2,3],[1,5,9],[1,6,11],[2,6,10],[3,6,9],[3,7,11],
                [4,5,6],[4,9,12],[5,6,7],[6,7,8],[6,10,12],[8,11,12],
                [9,10,11]]
    for i in range(Ndiscs):
        for bL in beamList:
            beamArr[nBeams] = [bL[0]+i*nVert,bL[1]+i*nVert,bL[2]+i*nVert,Kbeam]
            nBeams += 1
    
    #Record beamArray in .beam file
    print('nBeams = ',nBeams)
    f.write('%i\n' % nBeams)
    for i in range(nBeams):
        f.write('%i %i %i %.5e\n' % (int(beamArr[i,0]), int(beamArr[i,1]), int(beamArr[i,2]), beamArr[i,3]))
    f.close()
    
    
#------------------END OF MAIN----------------------
#SkeletonGenerator(0.375)
