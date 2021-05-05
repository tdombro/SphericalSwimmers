#!/bin/sh

path=$(pwd)

for dir in $( ls -d */ )
do
  cd $dir
  for dir2 in $(ls -d */)
  do
    cd $dir2
    #bsub < export.bsub
    /Applications/VisIt.app/Contents/Resources/bin/visit -nowin -cli -s $path/$dir/$dir2/ExportResampledDataVTKScript.py
    cd ..
    
    # Commenting this out because sbatch is cluster specific
    #sbatch sbatch_temp.sh infile.py
    
  done
  cd ..
done

