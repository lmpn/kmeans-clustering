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


touch /home/a77763/log/output.log 
for i in {1..4}; do
    export OMP_NUM_THREADS=$i
    /home/a77763/bin/run > /home/a77763/log/output.log
done
