#!/bin/sh

ThetaList="0.0 22.5 45.0 67.5 90.0 112.5 135.0 157.5 180.0 202.5 225.0 247.5 270.0 292.5 315.0 337.5"
HxList="-0.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5"
HyList="-13 -11 -9 -7 -5 -3 -1 1 3 5 7 9"

for Theta in $ThetaList
do
    cd 'Theta'$Theta
    for Hx in $HxList
    do
        cd 'Hx'$Hx
        for Hy in $(ls -d */)
        do
            echo $Hy
            Hy2=${Hy::${#Hy}-1}
            echo $(PWD)
            cd $Hy2
            echo $(PWD)
            cd $(PWD)'/../../..'
            echo $(PWD)
            sed 's/THETA/Theta'$Theta'/g;s/HX/Hx'$Hx'/g;s/HY/'$Hy2'/g' scriptSSL.sh > 'Theta'$Theta'/Hx'$Hx'/'$Hy2'/script.sh'
            cd 'Theta'$Theta'/Hx'$Hx
        done
        cd ../
    done
    cd ../
done
