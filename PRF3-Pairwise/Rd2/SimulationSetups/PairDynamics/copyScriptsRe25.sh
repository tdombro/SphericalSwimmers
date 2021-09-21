#!/bin/sh

ThetaList="0.0 22.5 45.0 67.5 90.0 112.5 135.0 157.5 180.0 202.5 225.0 247.5 270.0 292.5 315.0 337.5"
#HxList="-0.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5"
#HyList="-13 -11 -9 -7 -5 -3 -1 1 3 5 7 9"
ReList="25"

for Re in $ReList
do
    cd 'Re'$Re
    for Theta in $ThetaList
    do
        cd 'Theta'$Theta'/Hx8.5/'
        for Hy in $(ls -d */)
        do
            Hy2=${Hy::${#Hy}-1}
            cd $Hy2
            cd $(PWD)'/../../../..'
            #echo $(PWD)
            echo 'Re'$Re': T'$Theta': '$Hy2
            sed 's/THETA/Theta'$Theta'/g;s/HX/8.5/g;s/HY/'$Hy2'/g;s/RE/Re'$Re'/g' script.sh > 'Re'$Re'/Theta'$Theta'/Hx8.5/'$Hy2'/script.sh'
            sed 's/RE_VALUE/'$Re'/g' input2D > 'Re'$Re'/Theta'$Theta'/Hx8.5/'$Hy2'/input2D'
            cd 'Re'$Re'/Theta'$Theta'/Hx8.5/'

        done
        cd ../../
    done
    cd ../
done
