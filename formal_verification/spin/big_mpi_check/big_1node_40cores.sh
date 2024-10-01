#!/bin/bash
#SBATCH -J Armon_SPIN_1node
#SBATCH -e ./Armon_SPIN_1node_%j.err
#SBATCH -o ./Armon_SPIN_1node_%j.out
#SBATCH -c 4
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 64:00:00
#SBATCH --signal=INT
#SBATCH --distribution=cyclic:cyclic

cd ./run_1node
rm -f ./swarm_done_s* ./swarm_times_up 
srun -n 32 -c 4 -- bash -c "sh ./spin_script_\$SLURM_PROCID > spin_script_\$SLURM_PROCID.out || true"
