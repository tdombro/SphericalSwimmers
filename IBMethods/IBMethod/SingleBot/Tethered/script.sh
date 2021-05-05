#!/bin/bash
#SBATCH --job-name=TetheredRe2.5
#SBATCH --ntasks=64
#SBATCH --time=3-0
#SBATCH --partition=528_queue
#SBATCH --output=IB2D.out

mpirun ./main2d input2D
