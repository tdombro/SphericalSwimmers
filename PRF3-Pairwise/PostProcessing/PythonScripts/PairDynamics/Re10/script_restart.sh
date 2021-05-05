#!/bin/sh
#SBATCH --job-name=PD_Theta337.5_Hx12.5_Hy9_Re10
#SBATCH --ntasks=1
#SBATCH --time=7-0
#SBATCH --partition=general

./main2d.exe input2D_restart restart_IB2d 2000000