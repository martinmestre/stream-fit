#!/bin/bash

#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:05:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your email address>

julia test_distributed.jl