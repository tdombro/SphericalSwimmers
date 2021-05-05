#!/bin/sh

#SBATCH --job-name=SweepCIBrsl1.0v2.02amp0.8
#SBATCH --ntasks=1
#SBATCH --time=7-0
#SBATCH --partition=general
#SBATCH --output=%x.out

./main2d input2D > 1bot.log
