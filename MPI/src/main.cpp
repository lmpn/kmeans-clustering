#include <utils.h>
#include <mpi.h>
#include <kmeansCluster.h>

using namespace std;

// kmeans PAR/SEQ REPS SIZE DATASET PAPI_OPT

int main(int argc, char *argv[])
{
    int repetitions;
    int size;
    int clusters;
    char const *filename = NULL;
    char const *mode = NULL;
    char const *papiOpt = NULL;
    double *xcomp = NULL;
    double *ycomp = NULL;
    int *sets = NULL;
    if (argc < 6)
    {
        fprintf(stderr, "Usage: ./bin/kmeans PAR|SEQ #REPETITIONS #CLUSTERS #SIZE DATASET_PATH \n");
        return -1;
    }

    mode = argv[1];
    repetitions = atoi(argv[2]);
    clusters = atoi(argv[3]);
    size = atoi(argv[4]);
    filename = argv[5];
    xcomp = (double *)_mm_malloc(sizeof(double) * size,64);
    ycomp = (double *)_mm_malloc(sizeof(double) * size,64);
    int utils_error = utils_read_dataset(filename, xcomp, ycomp);

    if (utils_error == -1)
    {
        return -1;
    }

    if (!strcmp(mode, PAR))
    {
        int myrank, nprocesses;
        double times[repetitions];
        double reptimes[repetitions*8];
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        for (int i = 0; i < repetitions; i++)
        {
            utils_clear_cache();
            double x;
            if (myrank == 0)
            {
		        utils_start_section_timer();
	        }
            kmc_mpi(clusters, size, xcomp, ycomp, myrank, nprocesses, &sets, &(reptimes[i*8]));
            if (myrank == 0)
            {
                printf("%llu\n",utils_stop_section_timer());
            }
        }
        MPI_Finalize();
        if (myrank != 0)
        {
            return 0;
        }
    }
    
    return 0;
}
