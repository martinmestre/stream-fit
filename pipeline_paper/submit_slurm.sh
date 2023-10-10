#!/bin/bash
#SBATCH --partition=short
#SBATCH -o short_test.out  #Nombre del archivo de salida
#SBATCH -J m56  #Nombre del trabajo
#SBATCH --nodes=2  #Numero de nodos para correr el trabajo
#SBATCH --ntasks-per-node=64  #Numero de procesos


#Prepara el ambiente de trabajo
ulimit -l unlimited

#Ejecuta el programa paralelo con MPI
# . /etc/profile
julia test_slurm_manager.jl
