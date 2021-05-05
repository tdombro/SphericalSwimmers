#!/usr/bin/env python2.7.12
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 13:30:12 2019

@author: thomas
"""

#This script will generate images of the spherobots and fluid flow up until
#the last time step that was recorded
#Specify 1) Dist bw spherobots (R) 2) Angle 3) Anti or Para 4) SSL or LSL
#Import databases (fluid, mesh) and Open them
#Add Omega pseudocolor
#Add Nbot many botlow, botup, and skeleton
#Window size based off of box boundaries
#Save Window Attibutes
#Save in specific directory '../Images/'+RLength+'/PI'+str(Theta)+'/'+str(ANTIoPARA)+'/'+str(SSLoLSL)+'/'

import os
import sys
import pathlib

RLength = "5"
Theta = sys.argv[1]
ANTIoPARA = sys.argv[2]
SSLoLSL = sys.argv[3]

#Theta = "0"
#ANTIoPARA = "Anti"
#SSLoLSL = "SSL"

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'

#Parameters to be specified
Nbots = 2
plotNumber = 0

#Database paths
collectivePath = dir_path+"../Structures/"+RLength+"/PI"+Theta+"/"+ANTIoPARA+"/"+SSLoLSL+"/"
mesh = collectivePath+"/viz2D/lag_data.visit"
fluid = collectivePath+"/viz2D/dumps.visit"

db = [mesh,fluid]

def SetMeshAttributes(meshName,plotNumber):
    SetActivePlots(plotNumber)
    #SetMeshAttributes
    MeshAtts = MeshAttributes()
    MeshAtts.legendFlag = 0
    MeshAtts.lineStyle = MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    if("skeleton" in meshName):
        MeshAtts.lineWidth = 1
        MeshAtts.meshColor = (0,255,255,255)
        MeshAtts.pointSizePixels = 4
    else:
        MeshAtts.lineWidth = 1
        MeshAtts.meshColor = (255, 255, 255, 255)
        MeshAtts.pointSizePixels = 2
    MeshAtts.meshColorSource = MeshAtts.MeshCustom  # Foreground, MeshCustom
    MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom
    MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
    MeshAtts.pointSize = 0.05
    MeshAtts.opaqueColor = (255, 255, 255, 255)
    MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
    MeshAtts.pointSizeVarEnabled = 0
    MeshAtts.pointSizeVar = "default"
    MeshAtts.pointType = MeshAtts.Sphere  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    MeshAtts.showInternal = 0
    #MeshAtts.pointSizePixels = 2
    MeshAtts.opacity = 1
    SetPlotOptions(MeshAtts)
    return

if __name__ == '__main__':

    #Open databases
    OpenDatabase("localhost:" + db[0])
    OpenDatabase("localhost:" + db[1])
    
    #Delete sublevel plot
    SetActivePlots(0)
    DeleteActivePlots()
    
    #Correlation between fluid and mesh
    CreateDatabaseCorrelation("Correlation1",('localhost:' + db[0], 'localhost:' + db[1]), 0)
    
    #Save Window Attributes (For images)
    s = SaveWindowAttributes()
    s.fileName = "Image"
    s.family = 1
    s.outputToCurrentDirectory = 0
    #pathlib.Path("../Images/"+RLength+"/PI"+Theta+"/"+ANTIoPARA+"/"+SSLoLSL).mkdir(parents=True, exist_ok=True)
    s.outputDirectory = dir_path+"../Images/"+RLength+"/PI"+Theta+"/"+ANTIoPARA+"/"+SSLoLSL+"/"
    s.format = s.PNG
    s.width, s.height = 1080, 1080
    s.saveTiled = 0
    SetSaveWindowAttributes(s)
    
    #SetActiveWindow(1)
    #Activate Fluid database
    ActivateDatabase("localhost:"+db[1])
    
    #Add Pseudocolor of Vorticity
    AddPlot("Pseudocolor", "Omega", 1, 0)
    plotNumber += 1
    SetActivePlots(0)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = -10
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 10
    PseudocolorAtts.centering = PseudocolorAtts.Nodal  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "difference"
    PseudocolorAtts.invertColorTable = 0
    PseudocolorAtts.legendFlag = 0
    SetPlotOptions(PseudocolorAtts)
    print('Vorticity Plot Created!')
    
    #Activate Mesh database
    ActivateDatabase("localhost:"+db[0])
    
    
    '''Sphere Meshes'''
    for idx in range(1,Nbots+1):
        meshNameSS = "botlow"+str(idx)+"_vertices"
        meshNameLS = "botup"+str(idx)+"_vertices"
        meshNameSk = "skeleton"+str(idx)+"_mesh"
        AddPlot("Mesh",meshNameLS, 1, 0)
        SetMeshAttributes(meshNameLS,plotNumber)
        plotNumber += 1
        AddPlot("Mesh",meshNameSS, 1, 0)
        SetMeshAttributes(meshNameSS,plotNumber)
        plotNumber += 1
        AddPlot("Mesh",meshNameSk, 1, 0)
        SetMeshAttributes(meshNameSk,plotNumber)
        #HideActivePlots()
        plotNumber += 1
    print('Sphere Meshes have been created!')
    DrawPlots()
    
    #Window Size/Position
    # Begin spontaneous state
    View2DAtts = View2DAttributes()
    View2DAtts.windowCoords = (-0.05, 0.05, -0.05, 0.05)
    View2DAtts.viewportCoords = (0.1, 0.9, 0.1, 0.9)
    View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
    View2DAtts.fullFrameAutoThreshold = 100
    View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.windowValid = 1
    SetView2D(View2DAtts)
    # End spontaneous state
    
    #Create Annotation Objects!
    #Title
    title = CreateAnnotationObject("Text2D")
    title.text = RLength+"R: PI"+Theta+": "+ANTIoPARA+": "+SSLoLSL
    title.position = (0.25, 0.95)
    title.fontBold = 1
    # Add a time slider in the lower left corner
    slider = CreateAnnotationObject("TimeSlider")
    slider.height = 0.07
    
    SaveSession(collectivePath+"Images.session")
    
    SetActiveTimeSlider("Correlation1")
    
    for state in range(TimeSliderGetNStates()):
        SetTimeSliderState(state)
        AnnotationAtts = AnnotationAttributes()
        AnnotationAtts.axes2D.visible = 0
        SetAnnotationAttributes(AnnotationAtts)
        SaveWindow()
        
    SaveSession(collectivePath+"Images.session")
    
    exit()
    
