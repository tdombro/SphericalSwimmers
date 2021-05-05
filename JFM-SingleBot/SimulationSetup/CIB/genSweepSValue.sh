#!/bin/sh

eps="1.0"
SList="2.0 3.0 4.0 5.0 10.0 20.0 30.0"

for sval in $SList
do
    mkdir -p S_Value/s$sval
    cd Template
    cp main.cpp PETScOptions.dat small.vertex large.vertex Makefile $(PWD)'/../S_Value/s'$sval'/'
    sed 's/EPS_VALUE/'$eps'/g;s/SVALUE/'$sval'/g' input3D > $(PWD)'/../S_Value/s'$sval'/input3D'
    sed 's/EPS_VALUE/'$eps'/g;s/SVALUE/'$sval'/g' script.sh > $(PWD)'/../S_Value/s'$sval'/script.sh'
    cd ../
done
