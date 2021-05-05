#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 16:40:50 2018

@author: thomas
"""

import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
cwd = os.getcwd()

Re = sys.argv[1]
Config = sys.argv[2]

Nbots = 2

missingString = sys.argv[3]
print 'missing string = ',missingString
missingList = missingString.split()
print 'missing List = ',missingList
missingValues = [int(a) for a in missingList]
print 'missing values = ',missingValues

collectivePath=dir_path+"../"+Config+"/SweepRe/Re"+Re+"/"
mesh = collectivePath+"/viz2D/lag_data.visit"
fluid = collectivePath+"/viz2D/dumps.visit"

db = [mesh,fluid]

OpenDatabase("localhost:"+db[0])
OpenDatabase("localhost:"+db[1])

#Correlation between fluid and mesh
CreateDatabaseCorrelation("Correlation1",('localhost:' + db[0], 'localhost:' + db[1]), 0)

#Delete Sublevel Plot
SetActivePlots(0)
DeleteActivePlots()

AddPlot("Vector", "U")
AddPlot("Pseudocolor","Omega")

AddOperator("Resample",1)
ResampleAtts = ResampleAttributes()
ResampleAtts.useExtents = 1
ResampleAtts.startX = -0.05
ResampleAtts.endX = 0.05
ResampleAtts.samplesX = 1024
ResampleAtts.startY = -0.05
ResampleAtts.endY = 0.05
ResampleAtts.samplesY = 1024
ResampleAtts.is3D = 0
ResampleAtts.startZ = -1.0
ResampleAtts.endZ = 1.0
ResampleAtts.samplesZ = 128
ResampleAtts.tieResolver = ResampleAtts.random  # random, largest, smallest
ResampleAtts.tieResolverVariable = "default"
ResampleAtts.defaultValue = 0
ResampleAtts.distributedResample = 1
ResampleAtts.cellCenteredOutput = 0
SetOperatorOptions(ResampleAtts)

DrawPlots()

SetActiveTimeSlider("Correlation1")

for state in missingValues:
    # Change time states
    SetTimeSliderState(state)
    DrawPlots()
    # Export the database (ALL TIMES)
    dbAtts = ExportDBAttributes()
    dbAtts.allTimes = 0 #ALL TIMES
    dbAtts.db_type = "VTK"
    dbAtts.dirname = collectivePath+"VTK/"
    #dbAtts.filename = "Re5.Time%04d" % (start+i)
    dbAtts.filename = "DATA%05d"%(state)
    # Tuple of variables to export. "default" means the plotted variable.
    # Note that we're also exporting the vector expression "vec"
    dbAtts.variables = ("default", "amr_mesh","levels", "patches","U","P")
    ExportDatabase(dbAtts)

SaveSession(collectivePath+"Export2VTK2D_specific.session")
	
exit()
'''
OpenDatabase("localhost:"+db[0])
#Sphere Meshes
for idx in range(1,Nbots+1):
    meshNameSS = "botlow"+str(idx)+"_vertices"
    meshNameLS = "botup"+str(idx)+"_vertices"
    AddPlot("Mesh",meshNameLS, 1, 0)
    AddPlot("Mesh",meshNameSS, 1, 0)
print('Sphere Meshes have been created!')
DrawPlots()

SetActivePlots((2,3,4,5))

for state in range(TimeSliderGetNStates()):
    # Change time states
    SetTimeSliderState(state)
    DrawPlots()
    # Export the database (ALL TIMES)
    dbAtts = ExportDBAttributes()
    dbAtts.allTimes = 0 #ALL TIMES
    dbAtts.db_type = "VTK"
    dbAtts.dirname = collectivePath+"VTK/"
    #dbAtts.filename = "Re5.Time%04d" % (start+i)
    dbAtts.filename = "MESH%05d"%(state)
    # Tuple of variables to export. "default" means the plotted variable.
    # Note that we're also exporting the vector expression "vec"
    dbAtts.variables = ("default","botlow1_vertices","botlow2_vertices","botup1_vertices","botup2_vertices")
    ExportDatabase(dbAtts)
    '''
