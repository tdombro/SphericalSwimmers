#!/bin/sh

Re=$1
String=$2

/nas/longleaf/home/tdombro/sfw/visit/visit2.12.3/src/bin/visit -nowin -cli -s $PWD/ExportResampledDataVTKScript_specific.py $Re "$String"
