#Dar nome ao processo
#PBS -N parl3 

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
export OMP_PROC_BIND=close
export VERSION=par
make clean
make
for t in 2 4 8 12 16 24 32
do
	export OMP_NUM_THREADS=$t
	./bin/kmeans_par par 5 10 1966080 datasets/input1966080.data
done
