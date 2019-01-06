#!/./bin/sh

#Dar nome ao processo
#PBS -N arcore

#Tempo maximo do processo
#PBS -l walltime=02:00:00

#PBS -l nodes=2:r641:ppn=32

#Fila de espera para ir
#PBS -q mei

#mandar mail no principio(b) no final(e) e em caso de aborto(a)
#PBS -m bea

#para onde mandar mails
#PBS -M pachedo@gmail.com


module load gcc/7.2.0
module load gnu/openmpi_eth/1.8.4
cd /home/a77763/PCP/Hybrid

make clean
make

mpirun --map-by ppr:1:core \
  -np 32 \
  -report-bindings \
  --mac btl self,sm,tcp \
  bin/rbcore_omp par 8 10 62914560 \
  datasets/input62914560.data
/
