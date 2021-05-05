#!/bin/sh        

ReList="0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 5.5 6.0 6.5 7.0 7.5 10.0 12.5 15.0 17.5 20.0 25.0 30.0 35.0 40.0 50.0 60.0"

for Re in $ReList
do
    #mkdir -p '../PosData/V/Re'$Re
    echo 'V: Re'$Re
    #sbatch -J 'PD_Extract_V_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' scriptVisit.sh $Re 'V'
    #mkdir -p '../V/SweepRe/Re'$Re'/VTK/AVG'
    #sbatch -J 'VTK_V_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptExportVTK.sh $Re 'V'
    #sbatch -J 'ExportAVG_V_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptExportAVG.sh 'V' $Re
    sbatch -J 'Fig1_V_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptGenFig1.sh 'V' $Re
done
