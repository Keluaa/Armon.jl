#!/bin/bash
#SBATCH -J Armon_SPIN
#SBATCH -e ./Armon_SPIN_long_check_%j.err
#SBATCH -o ./Armon_SPIN_long_check_%j.out
#SBATCH -c 128
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 24:00:00
#SBATCH --signal=INT

cd ./run_long_check
srun -n 1 -c 128 -- ./pan_bitstate -m10000000 -c10 -k1 -w34 -q -I

