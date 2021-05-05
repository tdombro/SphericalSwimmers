#!/bin/sh        

ReList="0.5 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 25.0 30.0 35.0 40.0 50.0 60.0"

for Re in $ReList
do
    mkdir -p '../PosData/O/Re'$Re
    echo $Theta
    sbatch -J 'PD_Extractor_O_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' scriptVisit.sh $Re 'O'
done
