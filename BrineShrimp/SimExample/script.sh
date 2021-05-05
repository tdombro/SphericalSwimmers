#!/bin/bash
#SBATCH --job-name=flow_tank
#SBATCH --ntasks=64
#SBATCH --time=1-0
#SBATCH --mem=100000
#SBATCH --partition=528_queue
#SBATCH --output IBFE3D.out

mpirun ./main3d input3d
