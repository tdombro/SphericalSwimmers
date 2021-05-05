#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 16:40:50 2018

@author: thomas
"""

import os
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
cwd = os.getcwd()

Re = os.path.basename(cwd)
print('Re = ',Re)

#OpenDatabase("localhost:"+dir_path+"viz2D/dumpsExp.visit")
#OpenDatabase("loacalhost:"+dir_path+"viz2D/dumpsComp.visit")
OpenDatabase("localhost:"+dir_path+"viz2D/dumpsWhole.visit")

DeleteActivePlots()

SetTimeSliderState(0)

AddPlot("Vector", "U")
AddPlot("Pseudocolor","Omega")

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

# Export the database (ALL TIMES)
dbAtts = ExportDBAttributes()
dbAtts.allTimes = 1 #ALL TIMES
dbAtts.db_type = "VTK"
dbAtts.dirname = dir_path+"VTK/"
#dbAtts.filename = "Re5.Time%04d" % (start+i)
#dbAtts.filename = "Exp"
#dbAtts.filename = "Comp"
dbAtts.filename = "ExpComp"
# Tuple of variables to export. "default" means the plotted variable.
# Note that we're also exporting the vector expression "vec"
dbAtts.variables = ("default", "amr_mesh","levels", "patches","U","Omega")
ExportDatabase(dbAtts)

#SaveSession(dir_path+"Export2VTKExp.session")
#SaveSession(dir_path+"Export2VTKComp.session")
SaveSession(dir_path+"Export2VTKWhole.session")

exit()
