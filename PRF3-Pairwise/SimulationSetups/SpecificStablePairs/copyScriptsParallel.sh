#!/bin/sh

RE_VALUE="0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 5.5 6.0 6.5 7.0 7.5 10.0 12.5 15.0"

for Re in $RE_VALUE
do
    mkdir -p 'SweepRe/Re'$Re
    cp Template/botlow* 'SweepRe/Re'$Re'/'
    cp Template/botup* 'SweepRe/Re'$Re'/'
    cp Template/skeleton* 'SweepRe/Re'$Re'/'
    cp Template/main.C 'SweepRe/Re'$Re'/'
    cp Template/Makefile 'SweepRe/Re'$Re'/'
    cp Template/update_springs* 'SweepRe/Re'$Re'/'
    #Change input2D Re_Value
    #Change script.sh Re_Value
    sed 's/RE_VALUE/'$Re'/g' Template/input2dSSL > 'SweepRe/Re'$Re'/input2D'
    sed 's/RE_VALUE/'$Re'/g' Template/scriptSSL.sh > 'SweepRe/Re'$Re'/script.sh'
done
