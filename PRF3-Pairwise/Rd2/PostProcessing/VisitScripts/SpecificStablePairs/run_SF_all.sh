#!/bin/sh        

ReList="20.0 35.0"

for Re in $ReList
do
    if [ $Re == "20.0"]; then
        timeVTK="15.0"
    elif [ $Re == "35.0"]; then
        timeVTK="10.0"
    fi
    sbatch -J 'ExtractVTK_SF_Re'$Re -t 1-0 -n 1 -p general -o '%x.out' scriptVisit.sh $Re 'SF' $timeVTK
done
