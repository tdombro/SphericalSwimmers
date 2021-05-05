#!/bin/sh
#SBATCH --job-name=NangiaBhalla_Standard
#SBATCH --ntasks=1
#SBATCH --time=7-0
#SBATCH --partition=general
#SBATCH --output=%x.out

./main2d input2D
