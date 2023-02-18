#!/bin/bash
#SBATCH --requeue
#SBATCH --partition=scavenge
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=sp_0202_P0
#SBATCH --mem-per-cpu=4096
#SBATCH --time=1-
#SBATCH -o vT.out

module load MATLAB/2020b
/ysm-gpfs/apps/software/dSQ/1.05/dSQBatch.py --suppress-stats-file --job-file /gpfs/ysm/project/ohern/yz974/OverDamped/com_in/0202_P0/JobList_sp.txt
