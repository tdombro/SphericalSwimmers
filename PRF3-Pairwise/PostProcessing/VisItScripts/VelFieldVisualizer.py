#!/usr/bin/env python2.7.12
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 13:30:12 2019

@author: thomas
"""

#This script will generate images of the spherobots and fluid flow up until
#the last time step that was recorded
#It will also record the position data for each sphere at each time step
#Specify 1) Theta 2) Anti or Parallel 3) SSL or LSL
#Import databases (fluid, mesh) and Open them
#Add Velocity field vectors (Uniform, 50000 of them, no magnitude scale, hot)
#Add operators (threshold, clip "cylinder")
#Add Omega pseudocolor (Use RYB color scheme and invert)
#Add botlow, botup, and skeleton for each swimmer 
#(Use sphere geometry for point type) pointsize = 0.0003
#Make spheres White in color
#Window size based off of box boundaries (Not sure yet)
#Save Window Attibutes
#For each time step, find centroids of spheres. Will correspond to centers of clips
#Save in specific directory '../Images/Re+str(Re)+'/'

import os
import sys
#import pathlib
#import numpy as np

#ReValue = sys.argv[1]

RLength="5"
Theta = sys.argv[1]
ANTIoPARA = sys.argv[2]
SSLoLSL = sys.argv[3]

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'

#Parameters to be specified
Nbots = 2
plotNumber = 0

#Open File
#f = open(dir_path+'../PosData/'+RLength+'/Antiphase/PI'+Theta+'/'+ANTIoPARA+'/'+SSLoLSL+'/pd.txt','w')
f = open(dir_path+'../PosData/'+RLength+'/PI'+Theta+'/'+ANTIoPARA+'/'+SSLoLSL+'/pd.txt','w')
f.write('aXU aYU aXL aYL bXU bYU bXL bYL time\n')

#Database paths
#collectivePath = dir_path+"../VisitFiles/Antiphase/PI"+Theta+"/"+ANTIoPARA+"/"+SSLoLSL+"/"
collectivePath = dir_path+"../VisitFiles/5/PI"+Theta+"/"+ANTIoPARA+"/"+SSLoLSL+"/"
mesh = collectivePath+"/viz2D/lag_data.visit"
fluid = collectivePath+"/viz2D/dumps.visit"

db = [mesh,fluid]

def SetMeshAttributes(meshName):
    #SetMeshAttributes
    MeshAtts = MeshAttributes()
    MeshAtts.legendFlag = 0
    MeshAtts.lineStyle = MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    MeshAtts.meshColorSource = MeshAtts.MeshCustom  # Foreground, MeshCustom
    MeshAtts.meshColor = (200, 200, 200, 255) #255 #200
    MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom
    MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
    MeshAtts.lineWidth = 1
    MeshAtts.opaqueColor = (255, 255, 255, 255)
    MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
    MeshAtts.pointSizeVarEnabled = 0
    MeshAtts.pointSizeVar = "default"
    MeshAtts.pointType = MeshAtts.SphereGeometry  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    MeshAtts.pointSize = 0.0003
    MeshAtts.showInternal = 0
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
    s.outputDirectory = dir_path+"../Images/Test5PI3AntiLSL_Grey/" #"Test_Grey" #"Test_Hot"
    s.format = s.PNG
    s.width, s.height = 1080, 1080
    s.saveTiled = 0
    SetSaveWindowAttributes(s)
    
    #SetActiveWindow(1)
    #Activate Fluid database
    ActivateDatabase("localhost:"+db[1])
    
    #Add Velocity Field Vectors
    AddPlot("Vector", "U", 1, 0)
    SetActivePlots(0)
    VectorAtts = VectorAttributes()
    VectorAtts.glyphLocation = VectorAtts.UniformInSpace  # AdaptsToMeshResolution, UniformInSpace
    VectorAtts.useStride = 0
    VectorAtts.nVectors = 50000
    VectorAtts.lineStyle = VectorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    VectorAtts.lineWidth = 0
    VectorAtts.scale = 0.03
    VectorAtts.headSize = 0.25
    VectorAtts.headOn = 1
    VectorAtts.scaleByMagnitude = 0
    VectorAtts.autoScale = 1
    VectorAtts.colorByMag = 1
    VectorAtts.colorTableName = "Greys" #"Greys" #"hot" 
    VectorAtts.useLegend = 0
    VectorAtts.invertColorTable = 0
    VectorAtts.vectorOrigin = VectorAtts.Tail  # Head, Middle, Tail
    VectorAtts.minFlag = 1
    VectorAtts.maxFlag = 1
    VectorAtts.limitsMode = VectorAtts.OriginalData  # OriginalData, CurrentPlot
    VectorAtts.min = 0
    VectorAtts.max = 0.003
    VectorAtts.geometryQuality = VectorAtts.High
    SetPlotOptions(VectorAtts)

    #Add Cylindrical Clips for each sphere
    SetActivePlots(0)
    #Swimmer 1
    #R
    AddOperator("Cylinder")
    CylinderAtts = CylinderAttributes()
    CylinderAtts.point1 = (-0.006, 0.0, -1)
    CylinderAtts.point2 = (-0.006, 0.0, 1)
    CylinderAtts.radius = 0.002
    CylinderAtts.inverse = 1
    SetOperatorOptions(CylinderAtts, 0)
    #r
    AddOperator("Cylinder")
    CylinderAtts = CylinderAttributes()
    CylinderAtts.point1 = (-0.006, -0.005, -1)
    CylinderAtts.point2 = (-0.006, -0.005, 1)
    CylinderAtts.radius = 0.001
    CylinderAtts.inverse = 1
    SetOperatorOptions(CylinderAtts, 1)
    #Swimmer 2
    #R
    AddOperator("Cylinder")
    CylinderAtts = CylinderAttributes()
    CylinderAtts.point1 = (0.006, -0.005, -1)
    CylinderAtts.point2 = (0.006, -0.005, 1)
    CylinderAtts.radius = 0.002
    CylinderAtts.inverse = 1
    SetOperatorOptions(CylinderAtts, 2)
    #r
    AddOperator("Cylinder")
    CylinderAtts = CylinderAttributes()
    CylinderAtts.point1 = (0.006, 0.0, -1)
    CylinderAtts.point2 = (0.006, 0.0, 1)
    CylinderAtts.radius = 0.001
    CylinderAtts.inverse = 1
    SetOperatorOptions(CylinderAtts, 3)
    
    #Add a Threshold operator to limit vectors shown
    '''AddOperator("Threshold", 0)
    ThresholdAtts = ThresholdAttributes()
    ThresholdAtts.outputMeshType = 0
    ThresholdAtts.listedVarNames = ("U_magnitude")
    ThresholdAtts.zonePortions = (1)
    ThresholdAtts.lowerBounds = (0.0005) #20x smaller than swimmer
    ThresholdAtts.upperBounds = (1e+37)
    ThresholdAtts.defaultVarName = "U"
    ThresholdAtts.defaultVarIsScalar = 0
    SetOperatorOptions(ThresholdAtts, 4)'''
    print('Velocity Field Vectors Created!')
    
    #Add Pseudocolor of Vorticity
    AddPlot("Pseudocolor", "Omega", 1, 0)
    SetActivePlots(1)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = -10
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 10
    PseudocolorAtts.centering = PseudocolorAtts.Nodal  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "difference" #"difference" #"RdYlBu"
    PseudocolorAtts.invertColorTable = 0 #0 #1
    PseudocolorAtts.legendFlag = 0
    SetPlotOptions(PseudocolorAtts)
    print('Vorticity Plot Created!')
    
    #Activate Mesh database
    ActivateDatabase("localhost:"+db[0])
    
    plotNumber = 2
    '''Sphere Meshes'''
    for idx in range(1,Nbots+1):
        meshNameSS = "botlow"+str(idx)+"_vertices"
        meshNameLS = "botup"+str(idx)+"_vertices"
        AddPlot("Mesh",meshNameLS, 1, 0)
        SetActivePlots(plotNumber)
        SetMeshAttributes(meshNameLS)
        plotNumber += 1
        AddPlot("Mesh",meshNameSS, 1, 0)
        SetActivePlots(plotNumber)
        SetMeshAttributes(meshNameSS)
        plotNumber += 1
    print('Sphere Meshes have been created!')
    DrawPlots()
    
    #Window Size/Position
    # Begin spontaneous state
    View2DAtts = View2DAttributes()
    maxWin = 0.025
    minWin = -1.0*maxWin
    View2DAtts.windowCoords = (minWin, maxWin, minWin, maxWin)
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
    title.text = "Antiphase: "+RLength+"R: PI"+Theta+": "+ANTIoPARA+": "+SSLoLSL
    title.position = (0.25, 0.95)
    title.fontBold = 1
    # Add a time slider in the lower left corner
    slider = CreateAnnotationObject("TimeSlider")
    slider.height = 0.07
    
    SaveSession(collectivePath+"Test5PI3AntiLSL_Grey.session")
    
    SetActiveTimeSlider("Correlation1")
    
    #centerL = np.zeros((Nbots,2))
    #centerS = np.zeros((Nbots,2))
    centerL = [[0.0,0.0],[0.0,0.0]]
    centerS = [[0.0,0.0],[0.0,0.0]]
    for state in range(TimeSliderGetNStates()):
        #Keep Axes Invisible (Test)
        SetTimeSliderState(state)
        AnnotationAtts = AnnotationAttributes()
        AnnotationAtts.axes2D.visible = 0
        SetAnnotationAttributes(AnnotationAtts)
        #Find Centroids of Large and Small spheres
        #Find Centers of Spheres
        plotNumber = 2
        for idx in range(Nbots):
            SetQueryFloatFormat("%g")
            SetActivePlots(plotNumber)
            queryL = Query("Centroid")
            print('query: ',queryL)
            centerL[idx] = GetQueryOutputValue()
            f.write('%.5e %.5e '%(centerL[idx][0],centerL[idx][1]))
            plotNumber += 1
            SetActivePlots(plotNumber)
            queryS = Query("Centroid")
            print('query: ',queryS)
            centerS[idx] = GetQueryOutputValue()
            f.write('%.5e %.5e '%(centerS[idx][0],centerS[idx][1]))
            plotNumber += 1
        time = 0.01*state
        f.write('%.5e\n'%time)
        #Use Centroids to change cylinder operators
        SetActivePlots(0)
        #Swimmer 1
        #R
        CylinderAtts = CylinderAttributes()
        CylinderAtts.point1 = (centerL[0][0], centerL[0][1], -1)
        CylinderAtts.point2 = (centerL[0][0], centerL[0][1], 1)
        CylinderAtts.radius = 0.002
        CylinderAtts.inverse = 1
        SetOperatorOptions(CylinderAtts, 0)
        #r
        CylinderAtts = CylinderAttributes()
        CylinderAtts.point1 = (centerS[0][0], centerS[0][1], -1)
        CylinderAtts.point2 = (centerS[0][0], centerS[0][1], 1)
        CylinderAtts.radius = 0.001
        CylinderAtts.inverse = 1
        SetOperatorOptions(CylinderAtts, 1)
        #Swimmer 2
        #R
        CylinderAtts = CylinderAttributes()
        CylinderAtts.point1 = (centerL[1][0], centerL[1][1], -1)
        CylinderAtts.point2 = (centerL[1][0], centerL[1][1], 1)
        CylinderAtts.radius = 0.002
        CylinderAtts.inverse = 1
        SetOperatorOptions(CylinderAtts, 2)
        #r
        CylinderAtts = CylinderAttributes()
        CylinderAtts.point1 = (centerS[1][0], centerS[1][1], -1)
        CylinderAtts.point2 = (centerS[1][0], centerS[1][1], 1)
        CylinderAtts.radius = 0.001
        CylinderAtts.inverse = 1
        SetOperatorOptions(CylinderAtts, 3)
        #Also Use centroids to change window zoom if need be
        maxCenterX = max(centerL[0][0],centerL[1][0],centerS[0][0],centerS[1][0])
        minCenterX = min(centerL[0][0],centerL[1][0],centerS[0][0],centerS[1][0])
        maxCenterY = max(centerL[0][1],centerL[1][1],centerS[0][1],centerS[1][1])
        minCenterY = min(centerL[0][1],centerL[1][1],centerS[0][1],centerS[1][1])
        maxSwim = max(abs(maxCenterX),abs(minCenterX),abs(maxCenterY),abs(minCenterY))
        #Check if swimmer is outside of window
        if(maxSwim + 0.01 > maxWin): #Swimmer is leaving window area. Make it larger
            maxWin += 0.002
            minWin = -1.0*maxWin
            print("\n\nmaxWin = %.3f\n\n"%maxWin)
        #Change Window Attributes
        View2DAtts = View2DAttributes()
        View2DAtts.windowCoords = (minWin, maxWin, minWin, maxWin)
        View2DAtts.viewportCoords = (0.1, 0.9, 0.1, 0.9)
        View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
        View2DAtts.fullFrameAutoThreshold = 100
        View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.windowValid = 1
        SetView2D(View2DAtts)   
        
        SaveWindow()
        
    SaveSession(collectivePath+"Test5PI3AntiLSL_Grey.session")

    f.close()

    exit()
    
