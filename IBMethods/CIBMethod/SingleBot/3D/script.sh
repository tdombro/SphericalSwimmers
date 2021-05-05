#!/bin/bash
#SBATCH --job-name=3D_3dX_eps1.0_s1.0
#SBATCH --ntasks=24
#SBATCH --time=7-0
#SBATCH --partition=bigmem
#SBATCH --output=%x.out

/nas/longleaf/home/tdombro/longleaf/sfw/linux/openmpi/2.1.1/bin/mpirun ./main3d input3D
