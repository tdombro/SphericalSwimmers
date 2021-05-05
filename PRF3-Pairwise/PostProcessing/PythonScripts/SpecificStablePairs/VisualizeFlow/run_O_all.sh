#!/bin/sh        

ReList="0.5 1.0 2.0 3.0 4.0 5.0 5.5 6.0 6.5 7.0 7.5 10.0 12.5 15.0 17.5 20.0 25.0 30.0 35.0 40.0 50.0 60.0"

for Re in $ReList
do
    #mkdir -p '../PosData/O/Re'$Re
    echo 'O: Re'$Re
    #sbatch -J 'PD_Extract_O_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' scriptVisit.sh $Re 'O'
    #mkdir -p '../O/SweepRe/Re'$Re'/VTK/AVG'
    #sbatch -J 'VTK_O_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptExportVTK.sh $Re 'O'
    #sbatch -J 'ExportAVG_O_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptExportAVG.sh 'O' $Re
    sbatch -J 'Fig1_O_Re'$Re -t1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptGenFig1.sh 'O' $Re
done

