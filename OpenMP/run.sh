#!/./bin/sh

#Dar nome ao processo
#PBS -N openmp 

#Tempo maximo do processo
#PBS -l walltime=02:00:00

#PBS -l nodes=4:r662

#Fila de espera para ir
#PBS -q mei

#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
#PBS -M a77763@alunos.uminho.pt
module load papi/5.5.0 && module load gcc/4.8.2
cd /home/a77763/PCP/OpenMP

echo Tempo versão sequencial inicial
./bin/kmeans_si seq 20 5 4194304 datasets/input4194304.data nopapi

echo Papi inicial
./bin/kmeans_spi seq 20 5 4194304 datasets/input4194304.data flops
./bin/kmeans_spi seq 20 5 4194304 datasets/input4194304.data l2mr
./bin/kmeans_spi seq 20 5 4194304 datasets/input4194304.data l3mr

echo Tempo versão sequencial final
./bin/kmeans_sf seq 20 5 4194304 datasets/input4194304.data nopapi

echo Papi final
./bin/kmeans_spf seq 20 5 4194304 datasets/input4194304.data flops
./bin/kmeans_spf seq 20 5 4194304 datasets/input4194304.data l2mr
./bin/kmeans_spf seq 20 5 4194304 datasets/input4194304.data l3mr

echo Percentagem de tempo por etapa
./bin/kmeans_sfs seq 20 5 4194304 datasets/input4194304.data nopapi

array=(1048576 2097152 4194304)
echo Tempo versão final com varios datasets e particoes
for p in 4 5 6 7 8 9 10
do
	for item in ${array[*]}
    do
        echo Tempo versão final com $item dataset e $p particoes 
       ./bin/kmeans_sf seq 20 $p $item datasets/input$item.data nopapi
    done
done
