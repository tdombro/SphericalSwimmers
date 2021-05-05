#!/bin/sh

RList="3 5 6 7"
#RList="3 5 6 7"
CONFIGList="Anti Parallel PerpL PerpS"
#CONFIGList="PerpL PerpS"
RE_VALUE="SSL LSL Stat"

for R in $RList
do
    cd $R
    for Angle in $(ls -d */)
    do
	echo $Angle
	Angle2=${Angle::${#Angle}-1}
	echo $Angle2
	#cd ../
        echo $(PWD)
	for Config in $CONFIGList
	do
            #Change job name and copy over to directory
	    sed 's/RVALUE/'$R'/g;s/ANGLE/'$Angle2'/g;s/CONFIG/'$Config'/g;s/REVALUE/SSL/g' ../scriptSSL.sh > ../$R'/'$Angle2'/'$Config'/SSL/scriptSSL.sh'
	    sed 's/RVALUE/'$R'/g;s/ANGLE/'$Angle2'/g;s/CONFIG/'$Config'/g;s/REVALUE/LSL/g' ../scriptLSL.sh > ../$R'/'$Angle2'/'$Config'/LSL/scriptLSL.sh'
	    sed 's/RVALUE/'$R'/g;s/ANGLE/'$Angle2'/g;s/CONFIG/'$Config'/g;s/REVALUE/Stat/g' ../scriptStat.sh > ../$R'/'$Angle2'/'$Config'/Stat/scriptStat.sh'
	#cd $R
	done
    done
    cd ../
done
