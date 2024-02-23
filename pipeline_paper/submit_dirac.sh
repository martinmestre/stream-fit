#!/bin/bash
#SBATCH --job-name=chi2full
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmestre@fcaglp.unlp.edu.ar
#SBATCH --partition=batch

### Cantidad de nodos a usar
### mono/gpu: 1
### multi:    2-8
#SBATCH --nodes=1

### Cores a utilizar por nodo = procesos por nodo * cores por proceso
### mono/gpu: <= 8
### multi:    16-20
### Cantidad de procesos a lanzar por nodo
#SBATCH --ntasks-per-node=48
### Cores por proceso (para MPI+OpenMP)
###SBATCH --cpus-per-task=1

### GPUs por nodo
### mono/multi:        0
### gpu (Tesla K20X):  1-2
### gpu (Tesla M2090): 1-3
###SBATCH --gres=gpu:0

### Tiempo de ejecucion. Formato dias-horas:minutos.
### mono/gpu: <= 7 días
### multi:    <= 4 días
#SBATCH --time 2-0:00

. /etc/profile

python fit_pot_from_IbataPolysGaiaDR2-data_chi2full_2.py
