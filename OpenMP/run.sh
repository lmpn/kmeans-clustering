#!/bin/sh

#Dar nome ao processo
#PBS -N openmp 

#Tempo maximo do processo
#PBS -l walltime=5:00

#PBS -l nodes=6:r662

#Fila de espera para ir
#PBS -q mei

#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
#PBS -M a77763@alunos.uminho.pt
module load papi/5.5.0 && module load gcc/4.8.2
cd /home/a77763/OpenMP
./bin/kmeans par 1 32 /home/a77763/OpenMP/datasets/input32.data
./bin/kmeans par 1 32 /home/a77763/OpenMP/datasets/input32.data
./bin/kmeans par 1 32 /home/a77763/OpenMP/datasets/input32.data