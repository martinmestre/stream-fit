#!/bin/bash
#SBATCH --partition=multi
#SBATCH -o optim_pot_m360.out  #Nombre del archivo de salida
#SBATCH -J optim_pot_m360  #Nombre del trabajo
#SBATCH --nodes=3  #Numero de nodos para correr el trabajo
#SBATCH --ntasks=3  #Numero de procesos
#SBATCH --cpus-per-task=64  #Numero de cpus por proceso

#Prepara el ambiente de trabajo
ulimit -l unlimited

#Ejecuta el programa paralelo con MPI
. /etc/profile
julia -p 192 optim_pot.jl
