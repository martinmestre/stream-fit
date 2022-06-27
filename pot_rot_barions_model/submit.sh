#!/bin/bash
###SBATCH --name=fit_pot-slice-rot-barions
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmestre@fcaglp.unlp.edu.ar
#SBATCH --partition=multi

### Cantidad de nodos a usar
### mono/gpu: 1
### multi:    2-8
#SBATCH --nodes=1

### Cores a utilizar por nodo = procesos por nodo * cores por proceso
### mono/gpu: <= 8
### multi:    16-20
### Cantidad de procesos a lanzar por nodo
###SBATCH --ntasks-per-node=
### Cores por proceso (para MPI+OpenMP)
###SBATCH --cpus-per-task=

### GPUs por nodo
### mono/multi:        0
### gpu (Tesla K20X):  1-2
### gpu (Tesla M2090): 1-3
###SBATCH --gres=gpu:0

### Tiempo de ejecucion. Formato dias-horas:minutos.
### mono/gpu: <= 7 días
### multi:    <= 4 días
#SBATCH --time 2-00:00

. /etc/profile


python fit_pot-slice-rot-barions_from_IbataPolysGaiaDR2-data.py
