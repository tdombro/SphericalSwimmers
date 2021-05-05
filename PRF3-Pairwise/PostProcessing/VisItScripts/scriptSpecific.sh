#!/bin/sh

Re=$1
Config=$2
String=$3

/nas/longleaf/home/tdombro/sfw/visit/visit2.12.3/src/bin/visit -nowin -cli -s $PWD/ExportResampledDataVTKScript_specific.py $Re $Config "$String"
