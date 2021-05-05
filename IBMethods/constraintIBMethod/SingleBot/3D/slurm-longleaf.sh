#!/bin/sh
#SBATCH --job-name=ConstraintIB_3D_Test
#SBATCH --ntasks=2
#SBATCH -t 1-00:00:00
#SBATCH -p general
#SBATCH -o %x.out
$HOME/longleaf/sfw2/linux/openmpi/4.0.2/bin/mpirun ./main2d input2d
