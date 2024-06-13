#!/bin/bash 
# 
#SBATCH --job-name=AR_grasp
# 
#SBATCH --mail-type=end
#SBATCH --mail-user=seppemic@students.zhaw.ch
#
#SBATCH --nodes=1 
#SBATCH --ntasks=24 
#SBATCH --time=4-00:00:00 
#SBATCH --partition=earth-1 
#SBATCH --constraint=rhel7
#SBATCH --mem=32G

module load USS/2020
module load gcc/7.3.0

python /cfs/earth/scratch/seppemic/TM/scripts/exec/run_grasp.py
