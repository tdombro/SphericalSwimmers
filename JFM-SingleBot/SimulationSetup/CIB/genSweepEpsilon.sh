#!/bin/sh

epsilonList="0.2 0.4 0.6 0.8 1.0"
S_VALUE="1.0"

for eps in $epsilonList
do
    mkdir -p epsilon/eps$eps
    cd Template
    cp main.cpp PETScOptions.dat small.vertex large.vertex Makefile $(PWD)'/../epsilon/eps'$eps'/'
    sed 's/EPS_VALUE/'$eps'/g;s/SVALUE/'$S_VALUE'/g' input3D > $(PWD)'/../epsilon/eps'$eps'/input3D'
    sed 's/EPS_VALUE/'$eps'/g;s/SVALUE/'$S_VALUE'/g' script.sh > $(PWD)'/../epsilon/eps'$eps'/script.sh'
    cd ../
done
