#!/bin/sh
#SBATCH --job-name=PD_Theta45.0_Hx5.5_Hy3_Re2
#SBATCH --ntasks=1
#SBATCH --time=7-0
#SBATCH --partition=general
#SBATCH --output=IB2D.out

./main2d.exe input2D
