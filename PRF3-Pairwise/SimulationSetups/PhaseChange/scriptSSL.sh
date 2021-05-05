#!/bin/sh
#SBATCH --job-name=PI_RVALUER_ANGLE_CONFIG_REVALUE
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --partition=general
#SBATCH --output=IB2D.out

./main2d.exe input2dSSL
