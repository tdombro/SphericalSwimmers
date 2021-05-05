#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 13:46:06 2018

@author: thomas
"""

import os
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
cwd = os.getcwd()

St = os.path.basename(cwd)
Re = os.path.split(os.path.dirname(cwd))[-1]

OpenDatabase("localhost:"+dir_path+"VTK/AVG.vtk")

DefineScalarExpression("AverageUx", "AverageVelocity[0]")

DefineScalarExpression("AverageUy", "AverageVelocity[1]")

#Begin Average Omega and Velocity Operations
AddPlot("Pseudocolor", "AverageOmega")
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = -10
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = 10
PseudocolorAtts.centering = PseudocolorAtts.Nodal  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "difference"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 1
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.tubeResolution = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.tubeRadiusVarEnabled = 0
PseudocolorAtts.tubeRadiusVar = ""
PseudocolorAtts.tubeRadiusVarRatio = 10
PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 0.125
PseudocolorAtts.endPointRadiusBBox = 0.05
PseudocolorAtts.endPointResolution = 10
PseudocolorAtts.endPointRatio = 5
PseudocolorAtts.endPointRadiusVarEnabled = 0
PseudocolorAtts.endPointRadiusVar = ""
PseudocolorAtts.endPointRadiusVarRatio = 10
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 0
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)

AddPlot("Vector", "AverageVelocity")
VectorAtts = VectorAttributes()
VectorAtts.glyphLocation = VectorAtts.UniformInSpace  # AdaptsToMeshResolution, UniformInSpace
VectorAtts.useStride = 0
VectorAtts.stride = 1
VectorAtts.nVectors = 10000
VectorAtts.lineStyle = VectorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
VectorAtts.lineWidth = 0
VectorAtts.scale = 0.05
VectorAtts.scaleByMagnitude = 0
VectorAtts.autoScale = 1
VectorAtts.headSize = 0.25
VectorAtts.headOn = 1
VectorAtts.colorByMag = 0
VectorAtts.useLegend = 0
VectorAtts.vectorColor = (0, 0, 0, 255)
VectorAtts.colorTableName = "Default"
VectorAtts.invertColorTable = 0
VectorAtts.vectorOrigin = VectorAtts.Tail  # Head, Middle, Tail
VectorAtts.minFlag = 0
VectorAtts.maxFlag = 0
VectorAtts.limitsMode = VectorAtts.OriginalData  # OriginalData, CurrentPlot
VectorAtts.min = 0
VectorAtts.max = 1
VectorAtts.lineStem = VectorAtts.Line  # Cylinder, Line
VectorAtts.geometryQuality = VectorAtts.Fast  # Fast, High
VectorAtts.stemWidth = 0.08
VectorAtts.origOnly = 1
VectorAtts.glyphType = VectorAtts.Arrow  # Arrow, Ellipsoid
SetPlotOptions(VectorAtts)

AddOperator("Box")
SetActivePlots(1)
BoxAtts = BoxAttributes()
BoxAtts.amount = BoxAtts.Some  # Some, All
BoxAtts.minx = -2
BoxAtts.maxx = 2
BoxAtts.miny = -2
BoxAtts.maxy = 2
BoxAtts.minz = 0
BoxAtts.maxz = 1
BoxAtts.inverse = 0
SetOperatorOptions(BoxAtts, 0)
SetActivePlots(0)
BoxAtts = BoxAttributes()
BoxAtts.amount = BoxAtts.Some  # Some, All
BoxAtts.minx = -2
BoxAtts.maxx = 2
BoxAtts.miny = -2
BoxAtts.maxy = 2
BoxAtts.minz = 0
BoxAtts.maxz = 1
BoxAtts.inverse = 0
SetOperatorOptions(BoxAtts, 0)

#Remove Database and User Info
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
AnnotationAtts.userInfoFont.scale = 1
AnnotationAtts.userInfoFont.useForegroundColor = 1
AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
AnnotationAtts.userInfoFont.bold = 0
AnnotationAtts.userInfoFont.italic = 0
AnnotationAtts.databaseInfoFlag = 0
SetAnnotationAttributes(AnnotationAtts)

DrawPlots()

# Begin spontaneous state
View2DAtts = View2DAttributes()
View2DAtts.windowCoords = (-2, 2, -2, 2)
View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
View2DAtts.fullFrameAutoThreshold = 100
View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.windowValid = 1
SetView2D(View2DAtts)
# End spontaneous state

SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = dir_path
SaveWindowAtts.fileName = "AvgFluidFlow"+Re+St
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1920
SaveWindowAtts.height = 1080
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.advancedMultiWindowSave = 0
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

#Begin Lineout Operations
DeleteAllPlots()
AddPlot("Pseudocolor", "AverageUx", 1, 0)
DrawPlots()
SetQueryFloatFormat("%g")

#Create a Lineout operation (Horizontal)
Query("Lineout", end_point=(4, 0.01, 0), num_samples=500, start_point=(0.0, 0.01, 0), use_sampling=0, vars=("AverageUx"))

SetActiveWindow(2)
CurveAtts = CurveAttributes()
CurveAtts.showLines = 0
CurveAtts.lineStyle = CurveAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
CurveAtts.lineWidth = 0
CurveAtts.showPoints = 1
CurveAtts.symbol = CurveAtts.Point  # Point, TriangleUp, TriangleDown, Square, Circle, Plus, X
CurveAtts.pointSize = 5
CurveAtts.pointFillMode = CurveAtts.Static  # Static, Dynamic
CurveAtts.pointStride = 1
CurveAtts.symbolDensity = 50
CurveAtts.curveColorSource = CurveAtts.Cycle  # Cycle, Custom
CurveAtts.curveColor = (255, 0, 0, 255)
CurveAtts.showLegend = 1
CurveAtts.showLabels = 1
CurveAtts.designator = "A"
CurveAtts.doBallTimeCue = 0
CurveAtts.ballTimeCueColor = (0, 0, 0, 255)
CurveAtts.timeCueBallSize = 0.01
CurveAtts.doLineTimeCue = 0
CurveAtts.lineTimeCueColor = (0, 0, 0, 255)
CurveAtts.lineTimeCueWidth = 0
CurveAtts.doCropTimeCue = 0
CurveAtts.timeForTimeCue = 0
CurveAtts.fillMode = CurveAtts.NoFill  # NoFill, Solid, HorizontalGradient, VerticalGradient
CurveAtts.fillColor1 = (255, 0, 0, 255)
CurveAtts.fillColor2 = (255, 100, 100, 255)
CurveAtts.polarToCartesian = 0
CurveAtts.polarCoordinateOrder = CurveAtts.R_Theta  # R_Theta, Theta_R
CurveAtts.angleUnits = CurveAtts.Radians  # Radians, Degrees
SetPlotOptions(CurveAtts)

SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = dir_path
SaveWindowAtts.fileName = "Horizontal."+Re+"."+St
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 1024
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.advancedMultiWindowSave = 0
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

DeleteAllPlots()
SetActiveWindow(1)
DeleteAllPlots()
AddPlot("Pseudocolor", "AverageUy", 1, 0)
DrawPlots()

#Create a Lineout operation (Vertical)
Query("Lineout", end_point=(0.01, 4.0, 0), num_samples=500, start_point=(0.01, 0.0, 0), use_sampling=0, vars=("AverageUy"))

SetActiveWindow(2)
SetPlotOptions(CurveAtts)

SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = dir_path
SaveWindowAtts.fileName = "Vertical."+Re+"."+St
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 1024
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.advancedMultiWindowSave = 0
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

DeleteAllPlots()
SetActiveWindow(1)
DeleteAllPlots()

SaveSession(dir_path+"AverageFlow.session")

DeleteAllPlots()
exit()
