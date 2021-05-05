#!/bin/sh

ReList="2.5 70"

path=$(pwd)

cd ..

for Re in $ReList
do
    dir=Re$Re
    cd $dir
    #bsub < export.bsub
    /Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s $path/../$dir/AverageExportedVTK.py
    cd ..
    
    # Commenting this out because sbatch is cluster specific
    #sbatch sbatch_temp.sh infile.py
done

