#!/bin/bash 
# 
#SBATCH --job-name=AR_ArPIP
# 
#SBATCH --mail-type=end
#SBATCH --mail-user=seppemic@students.zhaw.ch
#SBATCH --output=/cfs/earth/scratch/seppemic/TM/ArPIP_AR.txt
#
#SBATCH --nodes=1 
#SBATCH --ntasks=24 
#SBATCH --time=2-00:00:00 
#SBATCH --partition=earth-1 
#SBATCH --constraint=rhel7
#SBATCH --mem=64G

module load USS/2020
module load gcc/7.3.0

eval "$(conda shell.bash hook)"
conda activate MT_env

python /cfs/earth/scratch/seppemic/TM/scripts/exec/run_ArPIP.py
