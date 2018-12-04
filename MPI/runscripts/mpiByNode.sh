#!/./bin/sh

#Dar nome ao processo
#PBS -N Inicial

#Tempo maximo do processo
#PBS -l walltime=02:00:00

#PBS -l nodes=2:r641:ppn=32

#Fila de espera para ir
#PBS -q mei

#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
module load gcc/5.3.0
module load gnu/openmpi_eth/1.8.4
cd /home/a77211/trabalho/PCP/MPI


for datasets in input1966080.data input62914560.data
	do
	for processos in 2 4 6 8 10 12 16 20 24 30 32
		do
		mpirun -bynode -np $processos --mca btl self,sm,tcp bin/kmeans_mpi par 8 10 $datasets
	done
done
