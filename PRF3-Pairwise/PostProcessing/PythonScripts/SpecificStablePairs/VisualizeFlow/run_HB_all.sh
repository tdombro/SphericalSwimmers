#!/bin/sh        

ReList="0.5 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 25.0 30.0 35.0 40.0 50.0 60.0"

for Re in $ReList
do
    #mkdir -p '../PosData/HB/Re'$Re
    echo 'HB: Re'$Re
    #sbatch -J 'PD_Extract_HB_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' scriptVisit.sh $Re 'HB'
    #mkdir -p '../HB/SweepRe/Re'$Re'/VTK/AVG'
    #sbatch -J 'VTK_HB_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptExportVTK.sh $Re 'HB'
    #sbatch -J 'ExportAVG_HB_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptExportAVG.sh 'HB' $Re
    sbatch -J 'Fig1_HB_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptGenFig1.sh 'HB' $Re
done
