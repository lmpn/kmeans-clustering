#!/./bin/sh

#Dar nome ao processo
#PBS -N arnode

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
cd /home/a77763/PCP/MPI
export ALLRED=yes
export NODE=yes
make clean
make

for datasets in 1966080 62914560
do
	for processos in 2 4 8 12 16 24 32
	do
		mpirun --map-by node -np $processos -report-bindings --mca btl self,sm,tcp bin/arnode par 8 10 $datasets "datasets/input"$datasets".data" > "out/arnode"$processos"_"$datasets
	done
done
