#!/bin/bash
#SBATCH --job-name=100botSSL_11235
#SBATCH --ntasks=128
#SBATCH --time=3-0
#SBATCH --partition=528_queue
#SBATCH --output=IB2D.out

mpirun ./main2d input2dSSL
