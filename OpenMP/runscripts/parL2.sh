#!/./bin/sh

#Dar nome ao processo
#PBS -N XEON

#Tempo maximo do processo
#PBS -l walltime=02:00:00

#PBS -l nodes=1:ppn=32

#Fila de espera para ir
#PBS -q mei

#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
cd /home/a77763/lab3
for t in 2 4 8 16 32
do
	export OMP_THREAD_NUM=$t
    ./helloflops3o_xeon
done
