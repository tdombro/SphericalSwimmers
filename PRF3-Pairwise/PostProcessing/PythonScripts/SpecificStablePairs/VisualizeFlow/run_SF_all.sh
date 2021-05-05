#!/bin/sh        

ReList="0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 5.5 6.0 6.5 7.0 7.5 10.0 12.5 15.0"

for Re in $ReList
do
    #mkdir -p '../PosData/SF/Re'$Re
    echo 'SF: Re'$Re
    #sbatch -J 'PD_Extract_SF_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' scriptVisit.sh $Re 'SF'
    #mkdir -p '../SF/SweepRe/Re'$Re'/VTK/AVG'
    #sbatch -J 'VTK_SF_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptExportVTK.sh $Re 'SF'
    #sbatch -J 'ExportAVG_SF_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptExportAVG.sh 'SF' $Re
    sbatch -J 'Fig1_SF_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptGenFig1.sh 'SF' $Re
done
