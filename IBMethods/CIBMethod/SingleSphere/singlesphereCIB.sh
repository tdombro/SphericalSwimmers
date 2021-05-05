#!/bin/sh

#SBATCH --job-name=SweepSS_StST_VALUEReRE_VALUE
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=general
#SBATCH --output=%x.out

./main2d input2D > 1sphere.log
