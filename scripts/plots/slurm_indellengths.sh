#!/bin/bash 
# 
#SBATCH --job-name=Indellengths_plot
# 
#SBATCH --mail-type=end
#SBATCH --mail-user=sepe@zhaw.ch
#SBATCH --output=/cfs/earth/scratch/sepe/TM/Indellengths_plot.txt
#
#SBATCH --nodes=1 
#SBATCH --ntasks=24 
#SBATCH --time=2:00:00 
#SBATCH --partition=earth-1
#SBATCH --constraint=rhel8
#SBATCH --mem=32G

module load USS/2020
module load gcc/7.3.0
module load miniconda3/4.8.2

eval "$(conda shell.bash hook)"
conda activate MT_env

python /cfs/earth/scratch/sepe/TM/scripts/plots/plot_indel_length.py

conda deactivate