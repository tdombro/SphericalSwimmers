#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 21:34:38 2020

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

cwd_PYTHON = os.getcwd()

#CONSTANTS GRID PLACEMENT
RADIUSLARGE = 0.002
RADIUSSMALL = 0.001

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

# Python program to check if rectangles overlap 
class Point: 
    def __init__(self, x, y): 
        self.x = x 
        self.y = y 
        
class Rectangle:
    def __init__(self,x,y,angle,name):
        #self.x = x #center in x axis
        #self.y = y #center in y axis
        #self.center = np.array([self.x,self.y])
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
        #self.width = np.amax(self.diag[:,0]) - np.amin(self.diag[:,0])
        #self.height = np.amax(self.diag[:,1]) - np.amin(self.diag[:,1])
        
    def Projection(self):
        #Here we find the diagonal coordinates for rectangle 2's axes
        for idx in range(4):
            self.proj[idx,0] = np.dot(self.normX,self.diag[idx])
            self.proj[idx,1] = np.dot(self.normY,self.diag[idx])
        print(self.name)
        print(self.proj)
        self.proj_xmin, self.proj_xmax = np.amin(self.proj[:,0]), np.amax(self.proj[:,0])
        self.proj_ymin, self.proj_ymax = np.amin(self.proj[:,1]), np.amax(self.proj[:,1])
        #self.width_proj = np.amax(self.proj[:,0]) - np.amin(self.proj[:,0])
        #self.height_proj = np.amax(self.proj[:,1]) - np.amin(self.proj[:,1])
        
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
        
        '''if(abs(self.x - other.x) < abs(width_self + width_other)/2.0 and
           abs(self.y - other.y) < abs(height_self + height_other)/2.0):
            print('They collide and intersect!')
        else:
            print('They donot intersect')'''
            
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
    axRec.axis([-5,10,-5,10])
    
    figRec.tight_layout()
    plt.show()

if __name__ == '__main__':
    
    Theta = 3.0*np.pi/8.0
    disp = np.array([6.0,6.0])

    #Check if Rectangles overlap before creating files
    # Create the square relative to (0, 0)
    #Rectangle 1
    points1 = np.array([
        [-3.0, -4.0],
        [-3.0, 4.0],
        [3.0, 4.0],
        [3.0, -4.0],
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
    
    Rec1 = Rectangle(center1[0],center1[1],0.0,'Rec1')
    Rec2 = Rectangle(center2[0],center2[1],Theta,'Rec2')
    Rec1.Diagonals(points1)
    Rec2.Diagonals(points2)
    Rec1.Projection()
    Rec2.Projection()
    isOverlap = Rec1.Intersects(Rec2)
    

    '''#Projection of spherobot 1 onto spherobot 2
    #Find diagonal vectors from center of rectangle
    diagVec1 = np.zeros((4,2))
    diagVec1[0] = points1[0] - center1
    diagVec1[1] = points1[1] - center1
    diagVec1[2] = points1[2] - center1
    diagVec1[3] = points1[3] - center1
    diagVec2 = np.zeros((4,2))
    diagVec2[0] = points2[0] - center2
    diagVec2[1] = points2[1] - center2
    diagVec2[2] = points2[2] - center2
    diagVec2[3] = points2[3] - center2
    #Normal vector of spherobot 2
    normVec2x = np.array([np.cos(Theta),np.sin(Theta)])
    normVec2y = np.array([np.cos(Theta-np.pi/2.0),np.sin(Theta-np.pi/2.0)])
    print('normVec2 = ',normVec2)
    dotv1_n2 = [np.dot(diagVec1[idx],normVec2) for idx in range(4)]
    dotv2_n2 = [np.dot(diagVec2[idx],normVec2) for idx in range(4)]
    projVec1 = np.zeros((4,2))
    projVec2 = np.zeros((4,2))
    for idx in range(4):
        projVec1[idx] = dotv1_n2[idx]*normVec2
        projVec2[idx] = dotv2_n2[idx]*normVec2
    print(projVec1)
    print(projVec2)
    
    
    
    Rec1 = Rectangle(center1[0],center1[1],2.0*RADIUSSMALL,8.0*RADIUSSMALL,0.0,'Rec1')
    Rec2 = Rectangle(center2[0],center2[1],2.0*RADIUSSMALL,8.0*RADIUSSMALL,Theta,'Rec2')
    print(Rec2.intersects(Rec1))'''