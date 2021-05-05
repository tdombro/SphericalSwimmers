#!/bin/sh        

ReList="0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 5.5 6.0 6.5 7.0 7.5 10.0 12.5 15.0 17.5 20.0 25.0 30.0 35.0 40.0 50.0 60.0"

for Re in $ReList
do
    sbatch -J 'AVG_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' --mem-per-cpu=25000 scriptAVG.sh $Re
    #sbatch -J 'Plot_Re'$Re -t 0-2 -n 1 -p general -o '%x.out' scriptPlot.sh $Re
done
