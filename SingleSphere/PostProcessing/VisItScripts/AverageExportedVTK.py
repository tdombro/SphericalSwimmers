#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 16:33:50 2018

@author: thomas
"""

import os
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
cwd = os.getcwd()

St = os.path.basename(cwd)
Re = os.path.split(os.path.dirname(cwd))[-1]

OpenDatabase("localhost:"+dir_path+"VTK/RAW"+"*.vtk database")

DeleteActivePlots()
SetTimeSliderState(0)

DefineVectorExpression("AverageVelocity", "average_over_time(U, "+str(0)+", "+str(19)+", 1)")

DefineScalarExpression("AverageUMAG","sqrt(AverageVelocity[0]*AverageVelocity[0] + AverageVelocity[1]*AverageVelocity[1])")

DefineScalarExpression("AverageOmega", "average_over_time(Omega, "+str(0)+", "+str(19)+", 1)")

AddPlot("Vector", "AverageVelocity", 1, 1)
AddPlot("Pseudocolor","AverageOmega", 1, 1)
DrawPlots()

ExportDBAtts = ExportDBAttributes()
ExportDBAtts.allTimes = 0
ExportDBAtts.db_type = "VTK"
ExportDBAtts.db_type_fullname = "VTK_1.0"
ExportDBAtts.filename = "AVG"
ExportDBAtts.dirname = dir_path+"VTK/"
ExportDBAtts.variables = ("default","amr_mesh", "levels", "patches", "AverageVelocity", "AverageUMAG", "AverageOmega")
ExportDBAtts.writeUsingGroups = 0
ExportDBAtts.groupSize = 48
ExportDBAtts.opts.types = (0, 0)
ExportDBAtts.opts.help = ""
ExportDatabase(ExportDBAtts)

SaveSession(dir_path+"Average2D.session")

exit()
