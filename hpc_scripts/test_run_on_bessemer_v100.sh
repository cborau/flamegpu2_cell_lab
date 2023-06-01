#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=carlos.dcsv100.sh

# Request a DCS v100
#SBATCH --partition=dcs-gpu
#SBATCH --account=dcs-res
#SBATCH --gres=gpu:1

# 10 CPU cores (1/4th of the node) and 1 GPUs worth of memory < 1/4th of the node)
#SBATCH --cpus-per-task=4
#SBATCH --mem=34G
## Name of file to which standard output will be redirected
#SBATCH --output="carlos_run_output.out"
## Name of file to which the standard error stream will be redirected
#SBATCH --error="carlos_run_error.err"
## Send me emails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robert.chisholm@sheffield.ac.uk

## Default pre-job stuff
## srun /bin/hostname
cd $SLURM_SUBMIT_DIR

## Start job

# Load modules
module load Anaconda3/5.3.0
module load GCC/8.3.0
module load CUDAcore/11.0.2

# Activate conda environment
source activate pyflamegpu

nvidia-smi

python moving_boundaries_grid3D_diffusion_ensemble.py