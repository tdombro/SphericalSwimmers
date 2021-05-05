#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 12:58:48 2020

@author: thomas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:46:45 2019

@author: thomas
"""

import os, shutil, sys
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
RADIUSSMALL = 0.001

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
    
    #CM
    CM = np.array([0.0,-1.0e-3])
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
        vertUp[0,i-1] = float(linesup[i][0]) - CM[0]
        vertUp[1,i-1] = float(linesup[i][1]) - CM[1]
    print('VertUp[1,0] = ',vertUp[1,0])
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
        vertLow[0,i-1] = float(lineslow[i][0]) - CM[0]
        vertLow[1,i-1] = float(lineslow[i][1]) - CM[1]
    print('VertLow[1,0] = ',vertLow[1,0])
    #Read vertex file line by line and save the values in a list using the \n delimiter
    linesskel = [line.strip() for line in open('skeleton.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    linesskel = [line.split() for line in linesskel]
    nvertskel = int(linesskel[0][0])
    #Allocate Array for Skeleton Vertex Positions
    vertSkel = np.zeros((2,nvertskel))
    #Store Vertices
    for i in range(1,nvertskel+1):
        vertSkel[0,i-1] = float(linesskel[i][0]) - CM[0]
        vertSkel[1,i-1] = float(linesskel[i][1]) - CM[1]
    print('VertSkel[1,0] = ',vertSkel[1,0])
    
    nvert = [nvertskel,nvertup,nvertlow]
    vertList = [vertSkel,vertUp,vertLow]
    
    return (vertList,nvert)

def DisplaceSpherobots(vertList, nvert, structureNames, Hx, Hy, Theta):
    
    #Create filePath and working directory
    #strDir = '../Structures/PairDynamics/Hx'+str(round(Hx,1))+'/Hy'+str(int(Hy))+'/Theta'+str(round(Theta/np.pi,3))+'/'
    #pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
    #cwd_PATH = cwd_PYTHON + strDir
    
    #Generate Figure to show Pairwise Placement
    #fig = plt.figure(num=0,figsize=(4,4),dpi=120)
    #ax = fig.add_subplot(111)
    #ax.set_title('Pairwise Initial Configuration: \nHx = %.1f Hy = %.1f Theta = PI*%.3f m'%(Hx,Hy,Theta/np.pi))
    #ax.axis([-0.025,0.025,-0.025,0.025])
    
    #First Rotate based on Theta
    #Allocate Arrays
    #rotationMatrix = np.zeros((2,2))
    #Calculate rotation matrix
    #rotationMatrix[0,0] = np.cos(Theta)
    #rotationMatrix[0,1] = -1.0*np.sin(Theta)
    #rotationMatrix[1,0] = np.sin(Theta)
    #rotationMatrix[1,1] = np.cos(Theta)  
    
    #Displace Spherobot 2 where they are a distance (Hx, Hy) apart
    x1Arr = np.array([0.0,0.0])
    x2Arr = np.array([Hx*RADIUSSMALL,Hy*RADIUSSMALL])
    xList = [x1Arr,x2Arr]
    
    #Check if Rectangles overlap before creating files
    # Create the square relative to (0, 0)
    isOverlap = CheckRectangleOverlap(Theta,xList[1])
    
    '''if(not(isOverlap)):
        #Create filePath and working directory
        strDir = '../Structures/PairDynamics/Hx'+str(round(Hx,1))+'/Hy'+str(int(Hy))+'/Theta'+str(round(Theta/np.pi,3))+'/'
        pathlib.Path(strDir).mkdir(parents=True, exist_ok=True)
        cwd_PATH = cwd_PYTHON + strDir
        #Generate Figure to show Pairwise Placement
        fig = plt.figure(num=0,figsize=(4,4),dpi=120)
        ax = fig.add_subplot(111)
        ax.set_title('Pairwise Initial Configuration: \nHx = %.1f Hy = %.1f Theta = PI*%.3f m'%(Hx,Hy,Theta/np.pi))
        ax.axis([-0.025,0.025,-0.025,0.025])
        #Displace and Rotate Spherobots
        for idxBot in range(2):
            dispArr = xList[idxBot]
            #print(dispArr)
            for idxName in range(len(vertList)):
                name = structureNames[idxName]
                vertPos = vertList[idxName].copy()
                #Generate new .vertex files
                f = open(cwd_PATH+name+str(idxBot+1)+'.vertex','w')
                copyfile(name+'.spring',cwd_PATH+name+str(idxBot+1)+'.spring')                
                f.write('%i\n'%nvert[idxName]) #Add #vertices
                #Rotate Spherobot 2 by Theta
                if(idxBot == 1):
                    #Spherobot 2
                    for idxVert in range(nvert[idxName]):
                        vertPos[:,idxVert] = Translate(Rotate(vertPos[:,idx],Theta),dispArr)
                        #vertPos[:,idxVert] = rotationMatrix.dot(vertPos[:,idxVert]) + dispArr[:]
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
    fig.savefig(cwd_PATH+'InitConfig.png')
    fig.clf()
    plt.close()'''
    
    return isOverlap

def CheckRectangleOverlap(Theta,disp):
    #Rectangle 1
    points1 = np.array([
        [-2.0*RADIUSSMALL, -5.0*RADIUSSMALL],
        [-2.0*RADIUSSMALL, 3.0*RADIUSSMALL],
        [2.0*RADIUSSMALL, 3.0*RADIUSSMALL],
        [2.0*RADIUSSMALL, -5.0*RADIUSSMALL]
    ])
    print('points1 = ',points1)
    #Rectangle 2
    #Rotate and Translate rectangle
    points2 = np.zeros((4,2))
    for idx in range(4):
        points2[idx,:] = Translate(Rotate(points1[idx],Theta),disp)
        #points2[idx,:] = rotationMatrix.dot(points1[idx,:]) + xList[1]
    print('points2 = ',points2)
    
    #Centers of Rectangles
    center1 = [0.0,0.0]
    center2 = (0.5*(np.amin(points2[:,0])+np.amax(points2[:,0])),
               0.5*(np.amin(points2[:,1])+np.amax(points2[:,1])))
    centers = np.array([center1,center2])
    print('centers = ',centers)
    
    #Plot to visually check overlap
    PlotRectangles(points1,points2,centers)
    
    Rec1 = Rectangle(0.0,'Rec1')
    Rec2 = Rectangle(Theta,'Rec2')
    Rec1.Diagonals(points1)
    Rec2.Diagonals(points2)
    Rec1.Projection()
    Rec2.Projection()
    isOverlap = Rec1.Intersects(Rec2)
    
    return isOverlap

def PlotRectangles(points1,points2,centers):
    #Plot Rectangles to visually check overlap
    figRec = plt.figure(num=1,figsize=(4,4),dpi=120)
    axRec = figRec.add_subplot(111)
    #Rectangle 1
    axRec.plot([points1[0,0],points1[1,0]],[points1[0,1],points1[1,1]],c='orange')
    axRec.plot([points1[1,0],points1[2,0]],[points1[1,1],points1[2,1]],c='r')
    axRec.plot([points1[2,0],points1[3,0]],[points1[2,1],points1[3,1]],c='b')
    axRec.plot([points1[3,0],points1[0,0]],[points1[3,1],points1[0,1]],c='g')
    #Rectangle 2    
    axRec.plot([points2[0,0],points2[1,0]],[points2[0,1],points2[1,1]],c='orange')
    axRec.plot([points2[1,0],points2[2,0]],[points2[1,1],points2[2,1]],c='r')
    axRec.plot([points2[2,0],points2[3,0]],[points2[2,1],points2[3,1]],c='b')
    axRec.plot([points2[3,0],points2[0,0]],[points2[3,1],points2[0,1]],c='g')
    #Centers
    axRec.scatter(centers[:,0],centers[:,1],c='k')
    #Add Rectangle 2 axes
    x2 = np.linspace(-5.0,10.0,100)
    y2x = centers[1,1] -1.0*(x2 - centers[1,0])
    y2y = centers[1,1] +1.0*(x2 - centers[1,0])
    axRec.plot(x2,y2x,c='k',ls='--')
    axRec.plot(x2,y2y,c='k',ls='--')
    axRec.axis([-0.005,0.015,-0.006,0.015])
    
    figRec.tight_layout()
    plt.show()
    plt.close()
    
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

def Translate(xy, offset):
    return xy + offset
    
class Rectangle:
    def __init__(self,angle,name):
        self.angle = angle
        self.name = name
        self.normX = np.array([np.cos(Theta),np.sin(Theta)])
        self.normY = np.array([np.cos(Theta+np.pi/2.0),np.sin(Theta+np.pi/2.0)])
        self.diag = np.zeros((4,2)) #diagonal coordinates on Rec1 axes
        self.proj = np.zeros((4,2)) #Diagonal coordinates on Rec2 axes
        print('normX = ',self.normX)
        print('normY = ',self.normY)

    def Diagonals(self,points):
        #Here we find diagonal coordinates for the rectangle in normal x-y space
        for idx in range(4):
            #self.diag[idx] = points[idx] - self.center
            self.diag[idx] = points[idx]
        print(self.name)
        print(self.diag)
        self.diag_xmin, self.diag_xmax = np.amin(self.diag[:,0]), np.amax(self.diag[:,0])
        self.diag_ymin, self.diag_ymax = np.amin(self.diag[:,1]), np.amax(self.diag[:,1])
        
    def Projection(self):
        #Here we find the diagonal coordinates for rectangle 2's axes
        for idx in range(4):
            self.proj[idx,0] = np.dot(self.normX,self.diag[idx])
            self.proj[idx,1] = np.dot(self.normY,self.diag[idx])
        print(self.name)
        print(self.proj)
        self.proj_xmin, self.proj_xmax = np.amin(self.proj[:,0]), np.amax(self.proj[:,0])
        self.proj_ymin, self.proj_ymax = np.amin(self.proj[:,1]), np.amax(self.proj[:,1])
        
    def Intersects(self,other):
        #If all 4 are not separated, then there is overlap
        #If one of these is separated, then there is no overlap
        #x1 proj
        bool_x1 = other.diag_xmax < self.diag_xmin or self.diag_xmax < other.diag_xmin
        print('bool_x1 = ',bool_x1)
        #y1 proj
        bool_y1 = other.diag_ymax < self.diag_ymin or self.diag_ymax < other.diag_ymin
        print('bool_y1 = ',bool_y1)
        #x2 proj
        bool_x2 = other.proj_xmax < self.proj_xmin or self.proj_xmax < other.proj_xmin
        print('bool_x2 = ',bool_x2)
        #y2 proj
        bool_y2 = other.proj_ymax < self.proj_ymin or self.proj_ymax < other.proj_ymin
        print('bool_y2 = ',bool_y2)
        
        self.overlap = not(bool_x1 or bool_x2 or bool_y1 or bool_y2)
        if(self.overlap):
            #No side is separated. They are overlapping and intersecting
            print('Intersection and Overlap!')
            print(self.overlap)
        else:
            print('They are Separated!')
            print(self.overlap)
        return self.overlap
    
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
    
    #Mesh has been stored. Now we need to position them appropriately 
    #based on three parameters
    #1)H_x: The distance between the 2 CM in the x-dir
    #2)H_y: The dist b/w the 2 CM in the y-dir
    #3)Theta: The angle b/w the 2 swimmers (arrow points in dir of large sphere)
    
    arrH_x = np.linspace(2.5,12.5,11)
    arrH_y = np.linspace(-13.0,9.0,12)
    arrTheta = np.linspace(0.0,2.0*np.pi,17)
    arrTheta = arrTheta[:-1].copy()
    print(arrH_x)
    print(arrH_y)
    print(arrTheta/np.pi)
    #Store all viable positions and orientations
    viableDict = {'Hx':[],'Hy':[],'Theta':[]}
    viableConfigs = pd.DataFrame(data=viableDict)
    '''#Test Rectangle Overlap Check
    intHx = random.randint(0,10)
    intHy = random.randint(0,11)
    intTheta = random.randint(0,15)
    Hx = arrH_x[0]
    Hy = arrH_y[7]
    Theta = arrTheta[10]
    print('Hx = %.3f\tHy = %.3f\t Theta/PI = %.3f'%(Hx,Hy,180*Theta/np.pi))
    DisplaceSpherobots(vertList, nvert, structureNames, Hx, Hy, Theta)'''
    #Loop over all Configurations
    for idx in range(len(arrH_x)):
        Hx = arrH_x[idx]
        for jdx in range(len(arrH_y)):
            Hy = arrH_y[jdx]
            for kdx in range(len(arrTheta)):
                Theta = arrTheta[kdx]
                print('Hx = %.3f\tHy = %.3f\t Theta/PI = %.3f'%(Hx/RADIUSSMALL,Hy/RADIUSSMALL,Theta/np.pi))
                isOverlap = DisplaceSpherobots(vertList, nvert, structureNames, Hx, Hy,Theta)
                if(not(isOverlap)):
                    #They are separated! Store in dataframe
                    #config = np.array([Hx,Hy,Theta/np.pi])
                    config = {'Hx':[Hx],'Hy':[Hy],'Theta':[Theta/np.pi]}
                    data = pd.DataFrame(data=config)
                    viableConfigs = pd.concat([viableConfigs,data],ignore_index=True)
                    print(len(viableConfigs.Hx))
    
    
    
        
