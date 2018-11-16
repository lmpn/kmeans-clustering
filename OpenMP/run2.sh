#!/./bin/sh

#Dar nome ao processo
#PBS -N openmp2

#Tempo maximo do processo
#PBS -l walltime=02:00:00

#PBS -l nodes=24:r662

#Fila de espera para ir
#PBS -q mei

#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
#PBS -M a77763@alunos.uminho.pt
module load papi/5.5.0 && module load gcc/4.8.2
cd /home/a77763/PCP/OpenMP
echo Tempo versão sequencial final
./bin/kmeans_sf seq 5 5 4194304 datasets/input4194304.data nopapi

echo Papi final
./bin/kmeans_spf seq 5 5 4194304 datasets/input4194304.data l1mr
./bin/kmeans_spf seq 5 5 4194304 datasets/input4194304.data l2mr
./bin/kmeans_spf seq 5 5 4194304 datasets/input4194304.data l3mr

echo Percentagem de tempo por etapa
./bin/kmeans_sfs seq 5 5 4194304 datasets/input4194304.data nopapi

