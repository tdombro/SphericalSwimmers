#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:41:56 2019

@author: thomas
"""

import os
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'

meshSSL = "/Users/thomas/sfw/Paper/Kinematics/HydroForces/VisitFiles/Re2.5/viz2D/lag_data.visit"
fluidSSL = "/Users/thomas/sfw/Paper/Kinematics/HydroForces/VisitFiles/Re2.5/viz2D/dumps.visit"
meshLSL = "/Users/thomas/sfw/Paper/Kinematics/HydroForces/VisitFiles/Re70/viz2D/lag_data.visit"
fluidLSL = "/Users/thomas/sfw/Paper/Kinematics/HydroForces/VisitFiles/Re70/viz2D/dumps.visit"

db = [meshSSL,fluidSSL,meshLSL,fluidLSL]

#Start Times
timeText = [None]*5
timeSSL = 43
timeLSL = 71

OpenDatabase("localhost:" + db[0])
OpenDatabase("localhost:" + db[1])
OpenDatabase("localhost:" + db[2])
OpenDatabase("localhost:" + db[3])

SetActivePlots((0,1))
DeleteActivePlots()

CreateDatabaseCorrelation("CorrelationSSL",('localhost:' + db[0], 'localhost:' + db[1]), 0)
CreateDatabaseCorrelation("CorrelationLSL",('localhost:' + db[2], 'localhost:' + db[3]), 0)


SetWindowLayout(4)

s = SaveWindowAttributes()
s.fileName = "VFC"
s.family = 1
s.outputToCurrentDirectory = 0
s.outputDirectory = "/Users/thomas/sfw/Paper/Kinematics/HydroForces/VisitFiles/VFC"
s.format = s.PNG
s.saveTiled = 1
SetSaveWindowAttributes(s)

for idxWin in range(1,5):
    SetActiveWindow(idxWin)
    print('Window = ',idxWin)
    if(idxWin < 3): #Re = 2.5 SSL
        ActivateDatabase("localhost:" + db[0])
    else:
        ActivateDatabase("localhost:" + db[2])
    '''Sphere Meshes'''
    #LargeSphere
    AddPlot("Mesh", "botup1_vertices", 1, 0)
    #SmallSphere
    AddPlot("Mesh", "botlow1_vertices", 1, 0)
    SetActivePlots((0,1))
    #Apply Mesh Attributes
    #SetMeshAttributes
    MeshAtts = MeshAttributes()
    MeshAtts.legendFlag = 0
    MeshAtts.lineStyle = MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    MeshAtts.lineWidth = 0
    MeshAtts.meshColor = (255, 255, 255, 255)
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
    MeshAtts.pointSizePixels = 10
    MeshAtts.opacity = 1
    SetPlotOptions(MeshAtts)
    print('Sphere Meshes have been created!')
    DrawPlots()
    
    #Find Centers of Spheres
    SetQueryFloatFormat("%g")
    SetActivePlots(0)
    queryL = Query("Centroid")
    print('query: ',queryL)
    centerL = GetQueryOutputValue()
    SetActivePlots(1)
    queryS = Query("Centroid")
    print('query: ',queryS)
    centerS = GetQueryOutputValue()
    print('LS: ',centerL[1])
    print('SS: ',centerS[1])

    '''Velocity Field Vectors'''
    if(idxWin < 3): #Re = 2.5 SSL
        ActivateDatabase("localhost:" + db[1])
    else:
        ActivateDatabase("localhost:" + db[3])
    #Add Velocity Field Vectors
    AddPlot("Vector", "U", 1, 0)
    SetActivePlots(2)
    VectorAtts = VectorAttributes()
    VectorAtts.glyphLocation = VectorAtts.AdaptsToMeshResolution  # AdaptsToMeshResolution, UniformInSpace
    VectorAtts.useStride = 1
    VectorAtts.stride = 5
    VectorAtts.lineStyle = VectorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    VectorAtts.lineWidth = 0
    VectorAtts.scale = 0.1
    VectorAtts.headSize = 0.25
    VectorAtts.headOn = 1
    if(idxWin == 1 or idxWin == 3): #Re = 2.5 SSL
        VectorAtts.scaleByMagnitude = 1
        VectorAtts.autoScale = 1
        VectorAtts.colorByMag = 1
        VectorAtts.colorTableName = "inferno"
    else:
        VectorAtts.scaleByMagnitude = 0
        VectorAtts.autoScale = 0
        VectorAtts.colorByMag = 0
        VectorAtts.vectorColor = (0, 0, 0, 255)   
    VectorAtts.useLegend = 0
    VectorAtts.invertColorTable = 0
    VectorAtts.vectorOrigin = VectorAtts.Tail  # Head, Middle, Tail
    VectorAtts.minFlag = 1
    VectorAtts.maxFlag = 1
    VectorAtts.limitsMode = VectorAtts.OriginalData  # OriginalData, CurrentPlot
    VectorAtts.min = 0
    VectorAtts.max = 1.5
    SetPlotOptions(VectorAtts)

    #Add Clip of Radius R or r in velocity field
    SetActivePlots(2)
    #Radius R
    AddOperator("Clip", 0)
    ClipAtts = ClipAttributes()
    ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
    ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
    ClipAtts.center = (0, centerL[1], 0)
    ClipAtts.radius = 0.3
    ClipAtts.sphereInverse = 0
    SetOperatorOptions(ClipAtts, 0)
    #Radius r
    AddOperator("Clip", 0)
    ClipAtts = ClipAttributes()
    ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
    ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
    ClipAtts.center = (0, centerS[1], 0)
    ClipAtts.radius = 0.15
    ClipAtts.sphereInverse = 0
    SetOperatorOptions(ClipAtts, 1)
    #Resample Velocity Field
    AddOperator("Resample",0)
    ResampleAtts = ResampleAttributes()
    ResampleAtts.useExtents = 1
    ResampleAtts.startX = -4
    ResampleAtts.endX = 4
    ResampleAtts.samplesX = 512
    ResampleAtts.startY = -4
    ResampleAtts.endY = 4
    ResampleAtts.samplesY = 512
    ResampleAtts.is3D = 0
    ResampleAtts.startZ = 0
    ResampleAtts.endZ = 1
    ResampleAtts.samplesZ = 10
    ResampleAtts.tieResolver = ResampleAtts.random  # random, largest, smallest
    ResampleAtts.tieResolverVariable = "default"
    ResampleAtts.defaultValue = 0
    ResampleAtts.distributedResample = 1
    ResampleAtts.cellCenteredOutput = 0
    SetOperatorOptions(ResampleAtts,2)
    print('Velocity Field Created!')
    
    #Add Pseudocolor of Vorticity
    AddPlot("Pseudocolor", "Omega", 1, 0)
    SetActivePlots(3)
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
    
    '''WINDOW SIZE/ZOOM'''
    if(idxWin < 3):
        # Begin spontaneous state
        View2DAtts = View2DAttributes()
        View2DAtts.windowCoords = (-0.523793, 0.511746, -1.48907, -0.5063)
        View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
        View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
        View2DAtts.fullFrameAutoThreshold = 100
        View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.windowValid = 1
        SetView2D(View2DAtts)
        # End spontaneous state
    else:
        # Begin spontaneous state
        View2DAtts = View2DAttributes()
        View2DAtts.windowCoords = (-0.52383, 0.51171, 0.0922078, 1.07498)
        View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
        View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
        View2DAtts.fullFrameAutoThreshold = 100
        View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
        View2DAtts.windowValid = 1
        SetView2D(View2DAtts)
        # End spontaneous state
    print('Window has been fit/zoomed correctly')
    
    DrawPlots()
    
    #Create Annotation Objects!
    title = CreateAnnotationObject("Text2D")
    if(idxWin == 1):
        title.text = "Re = 2.5: Scaled Velocity Vectors"
    elif(idxWin == 2):
        title.text = "Re = 2.5: Unscaled Velocity Vectors"
    elif(idxWin == 3):
        title.text = "Re = 70: Scaled Velocity Vectors"
    else:
        title.text = "Re = 70: Unscaled Velocity Vectors"
    title.position = (0.25, 0.95)
    title.fontBold = 1
    timeText[idxWin] = CreateAnnotationObject("Text2D")
    timeText[idxWin].position = (0.01,0.01)
    timeText[idxWin].fontBold = 1

    SaveSession(dir_path+"VF_Figures.session")
    
for idxTime in range(0,20):
    for idxWin in range(1,5):
        SetActiveWindow(idxWin)
        print('Window = ',idxWin)
        if(idxWin < 3): #Re = 2.5 SSL
            ActivateDatabase("localhost:" + db[0])
            #Set correct time correlation
            SetActiveTimeSlider("CorrelationSSL")
            timestep = timeSSL + idxTime
            SetTimeSliderState(timestep)
            sliderText = 'Tau = ' + str(float((timestep - timeSSL)/20.0))
        else:
            ActivateDatabase("localhost:" + db[2])
            #Set correct time correlation
            SetActiveTimeSlider("CorrelationLSL")
            timestep = timeLSL + idxTime
            SetTimeSliderState(timestep)
            sliderText = 'Tau = ' + str(float((timestep - timeLSL)/20.0))
        timeText[idxWin].text = sliderText
        
        #Find Centers of Spheres
        SetQueryFloatFormat("%g")
        SetActivePlots(0)
        queryL = Query("Centroid")
        print('query: ',queryL)
        centerL = GetQueryOutputValue()
        SetActivePlots(1)
        queryS = Query("Centroid")
        print('query: ',queryS)
        centerS = GetQueryOutputValue()
        print('LS: ',centerL[1])
        print('SS: ',centerS[1])
        
        if(idxWin < 3): #Re = 2.5 SSL
            ActivateDatabase("localhost:" + db[1])
        else:
            ActivateDatabase("localhost:" + db[3])
        SetActivePlots(2)
        ClipAtts = ClipAttributes() 
        ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
        ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
        ClipAtts.radius = 0.30
        ClipAtts.sphereInverse = 0
        ClipAtts.center = (0, centerL[1], 0)
        SetOperatorOptions(ClipAtts, 0)
        ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
        ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
        ClipAtts.center = (0, centerS[1], 0)
        ClipAtts.radius = 0.15
        ClipAtts.sphereInverse = 0
        SetOperatorOptions(ClipAtts, 1)

    #2) Save Figure!
    SaveWindow()
    
SaveSession(dir_path+"VF_Figures.session")

exit()

