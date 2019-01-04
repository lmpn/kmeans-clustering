#!/./bin/sh

#Dar nome ao processo
#PBS -N Final

#Tempo maximo do processo
#PBS -l walltime=02:00:00

#PBS -l nodes=1:r641:ppn=32

#Fila de espera para ir
#PBS -q mei

#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
#PBS -M a77763@alunos.uminho.pt
module load papi/5.5.0 && module load gcc/7.2.0
cd /home/a77763/PCP/OpenMP
export VERSION=final
make clean
make
./bin/kmeans_sf seq 5 10 1966080 datasets/input1966080.data
./bin/kmeans_sf seq 5 10 62914560 datasets/input62914560.data
