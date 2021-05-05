#!/bin/sh
#SBATCH --job-name=Pair_O_ReRE_VALUE
#SBATCH --ntasks=1
#SBATCH --time=4-0
#SBATCH --partition=general
#SBATCH --output=IB2D.out

./main2d.exe input2D
