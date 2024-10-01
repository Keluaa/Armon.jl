#!/bin/bash
#SBATCH -J Armon_SPIN
#SBATCH -e ./Armon_SPIN_4nodes_%j.err
#SBATCH -o ./Armon_SPIN_4nodes_%j.out
#SBATCH -c 2
#SBATCH -n 256
#SBATCH -N 4
#SBATCH --exclusive
#SBATCH -t 64:00:00
#SBATCH --signal=INT
#SBATCH --distribution=cyclic:cyclic

cd ./run
rm -f ./swarm_done_s* ./swarm_times_up
srun -n 200 -c 2 -- bash -c "sh ./spin_script_\$SLURM_PROCID > spin_script_\$SLURM_PROCID.out || true"
