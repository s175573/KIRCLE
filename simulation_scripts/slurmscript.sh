#!/bin/bash
#
# Batch script to run commands on BioHPC cluster.
# Submit to cluster using `sbatch` command.
# Lines starting with # are comments, and will not be run.
# Lines starting with #SBATCH specify options for the scheduler.
# Lines that do not start with # or #SBATCH are commands that will run.

# Name for the job (visible in job queue & accounting tools)
#SBATCH --job-name kir_wxs_sim

# SLURM partition to run this job
#SBATCH -p 32GB       # partition (queue)
# Number of nodes
#SBATCH -N 1

# Memory (RAM) requirement/limit in MB.
#SBATCH --mem 28672      # Memory Requirement (MB)

# Wall time limit for the job (format is Days-H:M:S)
#SBATCH -t 2-0:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o kir_wxs_sim_%j.out
#SBATCH -e kir_wxs_sim_%j.err

# Send an email when the job status changes, to the specified address.
#SBATCH --mail-type ALL
#SBATCH --mail-user galen.gao@utsouthwestern.edu

# Load appropriate modules

module load python/3.6.4-anaconda
#module load samtools/gcc/1.8
module load parallel

PATH=$PATH:$HOME/bin:$HOME/bin/ncbi-blast-2.8.1+/bin/
export BLASTDB=$HOME/blastdb

# Run script
cd /project/bioinformatics/Li_lab/s430614/EM_sims/
bash genotyper_sims.sh KIR2DS1
