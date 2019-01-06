#!/./bin/sh

#Dar nome ao processo
#PBS -N ParL3

#Tempo maximo do processo
#PBS -l walltime=02:00:00

#PBS -l nodes=1:r662:ppn=48

#Fila de espera para ir
#PBS -q mei

#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
#PBS -M a77211@alunos.uminho.pt
module load papi/5.5.0 && module load gcc/7.2.0
cd /home/a77211/trabalho/pcp/OpenMP
export OMP_PLACES=1
export OMP_PROC_BIND=close
for t in 2 4 8 12 16 24
do
	export OMP_NUM_THREADS=$t
	./bin/kmeans_par par 5 10 786432 datasets/input786432.data
done
