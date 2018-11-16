#!/./bin/sh

#Dar nome ao processo
#PBS -N openmp3 

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


array=(1048576 2097152 4194304)
echo Tempo versão final com varios datasets e particoes
for item in ${array[*]}
do
        echo Tempo versão final com $item dataset e 7 particoes 
       ./bin/kmeans_sf seq 5 7 $item datasets/input$item.data nopapi
done
