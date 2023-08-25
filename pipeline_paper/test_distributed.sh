#!/bin/bash

#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=52G

julia test_distributed.jl