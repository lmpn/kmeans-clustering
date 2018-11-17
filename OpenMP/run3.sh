#!/./bin/sh

#Dar nome ao processo
#PBS -N openmp3 

#Tempo maximo do processo
#PBS -l walltime=02:00:00
#Fila de espera para ir
#PBS -q mei
#PBS -l nodes=4:r662:ppn=2
#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
#PBS -M a77763@alunos.uminho.pt
module load papi/5.5.0 && module load gcc/4.8.2
cd /home/a77763/PCP/OpenMP


echo Percentagem de tempo por etap
./bin/kmeans_sfs seq 5 10 1048576 datasets/input1048576.data
echo Percentagem de tempo por etapa
./bin/kmeans_sfs seq 5 10 2097152 datasets/input2097152.data
echo Percentagem de tempo por etapa
./bin/kmeans_sfs seq 5 10 4194304 datasets/input4194304.data
