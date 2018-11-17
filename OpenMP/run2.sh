#!/./bin/sh

#Dar nome ao processo
#PBS -N openmp2

#Tempo maximo do processo
#PBS -l walltime=02:00:00
#PBS -l nodes=4:r662:ppn=2
#Fila de espera para ir
#PBS -q mei

#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
#PBS -M a77763@alunos.uminho.pt
module load papi/5.5.0 && module load gcc/4.8.2
cd /home/a77763/PCP/OpenMP



echo Tempo versão sequencial final 1048576
./bin/kmeans_sf seq 5 10 1048576 datasets/input1048576.data
echo Papi final 1048576
./bin/kmeans_spf seq 5 10 1048576 datasets/input1048576.data 

echo Tempo versão sequencial final 2097152
./bin/kmeans_sf seq 5 10 2097152 datasets/input2097152.data
echo Papi final 2097152
./bin/kmeans_spf seq 5 10 2097152 datasets/input2097152.data 

echo Tempo versão sequencial final 4194304
./bin/kmeans_sf seq 5 10 4194304 datasets/input4194304.data
echo Papi final 4194304
./bin/kmeans_spf seq 5 10 4194304 datasets/input4194304.data


