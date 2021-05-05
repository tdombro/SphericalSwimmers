#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 16:40:50 2018

@author: thomas
"""

import os
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
cwd = os.getcwd()

St = os.path.basename(cwd)
Re = os.path.split(os.path.dirname(cwd))[-1]
print('Re = ',Re)

start=0
end=19

OpenDatabase("localhost:"+dir_path+"viz2D/dumps.visit")

DeleteActivePlots()

SetTimeSliderState(0)


AddPlot("Vector", "U")
AddPlot("Pseudocolor","Omega")
'''AddOperator("Slice",1)
SliceAtts = SliceAttributes()
SliceAtts.originType = SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
SliceAtts.originPoint = (0, 0, 0)
SliceAtts.originIntercept = 0
SliceAtts.originPercent = 0
SliceAtts.originZone = 0
SliceAtts.originNode = 0
SliceAtts.normal = (0, 0, 1)
SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
SliceAtts.upAxis = (0, 1, 0)
SliceAtts.project2d = 1
SliceAtts.interactive = 1
SliceAtts.flip = 0
SliceAtts.originZoneDomain = 0
SliceAtts.originNodeDomain = 0
SliceAtts.meshName = "amr_mesh"
SliceAtts.theta = 0
SliceAtts.phi = 90
SetOperatorOptions(SliceAtts)'''

AddOperator("Resample",1)
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
SetOperatorOptions(ResampleAtts)

DrawPlots()

#for i in range(20):
# Change time states
#SetTimeSliderState(start+i)
#DrawPlots()
# Export the database (ALL TIMES)
dbAtts = ExportDBAttributes()
dbAtts.allTimes = 1 #ALL TIMES
dbAtts.db_type = "VTK"
dbAtts.dirname = dir_path+"VTK/"
#dbAtts.filename = "Re5.Time%04d" % (start+i)
dbAtts.filename = "RAW"
# Tuple of variables to export. "default" means the plotted variable.
# Note that we're also exporting the vector expression "vec"
dbAtts.variables = ("default", "amr_mesh","levels", "patches","U")
ExportDatabase(dbAtts)

SaveSession(dir_path+"Export2VTK2D.session")
	
exit()
