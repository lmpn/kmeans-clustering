#include <omp.h>
#include <stdio.h>

int main() {
	int i=0;
	int j=0;
	int thread[4];
	#pragma omp parallel
	{
		do{
			#pragma omp for
				for(int i=0;i<4;i++)
					thread[i]=omp_get_thread_num()*(i+1);
			#pragma omp master
			{
					printf("Single\n");
					j++;

			}

			#pragma omp barrier

		} while(j!=2);



	}

	for(int a =0;a<4;a++) {
		printf("Thread num: %d\n",thread[a] );
	}
	return 0;
}