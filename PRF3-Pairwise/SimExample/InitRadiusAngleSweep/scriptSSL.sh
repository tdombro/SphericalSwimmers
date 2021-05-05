#!/bin/sh
#SBATCH --job-name=5R_PI2_Anti_SSL
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --partition=general
#SBATCH --output=IB2D.out

./main2d.exe input2dSSL
