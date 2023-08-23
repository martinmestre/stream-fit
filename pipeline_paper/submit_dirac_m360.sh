#!/bin/bash
#SBATCH --partition=debug
#SBATCH -o optim_pot_m360.out  #Nombre del archivo de salida
#SBATCH -J optim_pot_m360  #Nombre del trabajo
#SBATCH --nodes=1  #Numero de nodos para correr el trabajo
#SBATCH --ntasks=4  #Numero de procesos
#SBATCH --tasks-per-node=4  #Numero de tareas por nodo

#Prepara el ambiente de trabajo
ulimit -l unlimited

#Ejecuta el programa paralelo con MPI
. /etc/profile
prun julia -p 4 optim_pot.jl
