#include <kmeansCluster.h>
using namespace std;
#include <mpi.h>
#include <mm_malloc.h>

int *kmc_seq_final(int clusters, int size, double *xcomp, double *ycomp)
{

    std::mt19937 rng;
    uint32_t seed_val;
    rng.seed(seed_val);
    int *sets = (int *)calloc(size, sizeof(int));
    int point_set_idx = 0;
    double *sets_counter = (double *)calloc(clusters, sizeof(double));
    double max = -DBL_MAX;
    double error = DBL_MAX;
    double c_error = DBL_MAX;

    double centroid_x[clusters];
    double centroid_y[clusters];

    for (size_t i = 0; i < size; i++)
    {
        if (max < xcomp[i])
            max = xcomp[i];
        if (max < ycomp[i])
            max = ycomp[i];
    }
    uniform_real_distribution<double> urd_g(0, max);
    for (int i = 0; i < clusters; i++)
    {
        centroid_x[i] = urd_g(rng);
        centroid_y[i] = urd_g(rng);
    }

    do
    {
        c_error = error;
        error = 0.0;

        for (int point_idx = 0; point_idx < size; point_idx++)
        {
            int current_point_cluster_idx = -1;
            double minimun_distance = DBL_MAX;

            for (int cluster_idx = 0; cluster_idx < clusters; cluster_idx++)
            {
                double dx = centroid_x[cluster_idx] - xcomp[point_idx];
                double dy = centroid_y[cluster_idx] - ycomp[point_idx];
                double centroid_point_distance = dx * dx + dy * dy;
                if (minimun_distance > centroid_point_distance)
                {
                    minimun_distance = centroid_point_distance;
                    current_point_cluster_idx = cluster_idx;
                }
            }

            sets_counter[current_point_cluster_idx] += 1.0;
            sets[point_idx] = current_point_cluster_idx;
        }

        for (int cluster_idx = 0; cluster_idx < clusters; cluster_idx++)
        {
            error = error - centroid_y[cluster_idx] - centroid_x[cluster_idx];
            centroid_x[cluster_idx] = 0.0;
            centroid_y[cluster_idx] = 0.0;
            sets_counter[cluster_idx] = 1 / sets_counter[cluster_idx];
        }

        for (int i = 0; i < size; i++)
        {
            int point_set_idx = sets[i];
            double set_size = sets_counter[point_set_idx];
            centroid_x[point_set_idx] += xcomp[i] * set_size;
            centroid_y[point_set_idx] += ycomp[i] * set_size;
        }

        utils_start_section_timer();
        for (int k = 0; k < clusters; k++)
        {
            sets_counter[k] = 0;
            error = error + centroid_x[k] + centroid_y[k];
        }

    } while (error != c_error);

    return sets;
}


void kmc_mpi(int clusters, int size, double *xcomp, double *ycomp, int myrank, int nprocesses, int **result, double * times)
{



//utils_stop_section_timer
    /*
    * Declaring global variables
    */
    long long unsigned t[8] = {0,0,0,0,0,0,0,0};
    double centroid_x_global[clusters];
    double centroid_y_global[clusters];
    int *sets_global = (int *)calloc(size, sizeof(int));
    double *sets_counter_global = (double *)calloc(clusters, sizeof(double));
    int chunk_size = size / nprocesses;
    double error = DBL_MAX;
    double c_error = DBL_MAX;
    double start, end;
    int *sets_local = (int *)calloc(chunk_size, sizeof(int));
    double *xcomp_local = (double *)_mm_malloc(chunk_size * sizeof(double), 64);
    double *ycomp_local = (double *)_mm_malloc(chunk_size * sizeof(double), 64);
    if (myrank == 0)
    {
        utils_start_section_timer();
        std::mt19937 rng;
        uint32_t seed_val;
        rng.seed(seed_val);

        double max = -DBL_MAX;

        for (size_t i = 0; i < size; i++)
        {
            if (max < xcomp[i])
                max = xcomp[i];
            if (max < ycomp[i])
                max = ycomp[i];
        }

        uniform_real_distribution<double> urd_g(0, max);

        for (int i = 0; i < clusters; i++)
        {
            centroid_x_global[i] = urd_g(rng);
            centroid_y_global[i] = urd_g(rng);
        }
        printf("%llu;",utils_stop_section_timer());
    }

    /*
     * Declaring local variables
    */


    utils_start_section_timer();
    MPI_Scatter(xcomp, chunk_size, MPI_DOUBLE, xcomp_local, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(ycomp, chunk_size, MPI_DOUBLE, ycomp_local, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(centroid_x_global, clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(centroid_y_global, clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    printf("%llu;",utils_stop_section_timer());
    
    do
    {
        if (myrank == 0)
        {
            //phase1
            utils_start_section_timer();
            c_error = error;
            error = 0.0;
        }

            
        /**
         * Phase 2 Master + Workers
        */
        for (int point_idx = 0; point_idx < chunk_size; point_idx++)
        {
            int current_point_cluster_idx = -1;
            double minimun_distance = DBL_MAX;
            for (int cluster_idx = 0; cluster_idx < clusters; cluster_idx++)
            {
                double dx = centroid_x_global[cluster_idx] - xcomp_local[point_idx];
                double dy = centroid_y_global[cluster_idx] - ycomp_local[point_idx];
                double centroid_point_distance = dx * dx + dy * dy;
                if (minimun_distance > centroid_point_distance)
                {
                    minimun_distance = centroid_point_distance;
                    current_point_cluster_idx = cluster_idx;
                }
            }

            sets_counter_global[current_point_cluster_idx] += 1.0;
            sets_local[point_idx] = current_point_cluster_idx;
        }
        /**
         * Phase 3
         */

        if (myrank == 0)
        {
            
            times[2] += utils_stop_section_timer();

        }
#ifdef REDUCEBCAST
        if (myrank == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, sets_counter_global, clusters, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Reduce(sets_counter_global, sets_counter_global, clusters, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        MPI_Bcast(sets_counter_global, clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#ifdef ALLRED
        MPI_Allreduce(MPI_IN_PLACE, sets_counter_global, clusters, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        if (myrank == 0)
        {
            
            
            //COMM1
            times[3] += utils_stop_section_timer();;
            utils_start_section_timer();
            for (int cluster_idx = 0; cluster_idx < clusters; cluster_idx++)
            {
                error = error - centroid_y_global[cluster_idx] - centroid_x_global[cluster_idx];
                centroid_x_global[cluster_idx] = 0.0;
                centroid_y_global[cluster_idx] = 0.0;
                sets_counter_global[cluster_idx] = 1 / sets_counter_global[cluster_idx];
            }
        }
        else
        {
            for (int cluster_idx = 0; cluster_idx < clusters; cluster_idx++)
            {
                centroid_x_global[cluster_idx] = 0.0;
                centroid_y_global[cluster_idx] = 0.0;
                sets_counter_global[cluster_idx] = 1 / sets_counter_global[cluster_idx];
            }
        }
        for (int i = 0; i < chunk_size; i++)
        {
            int point_set_idx = sets_local[i];
            double set_size = sets_counter_global[point_set_idx];
            centroid_x_global[point_set_idx] += xcomp_local[i] * set_size;
            centroid_y_global[point_set_idx] += ycomp_local[i] * set_size;
        }
        if (myrank == 0)
        {
            //PHASE32
            
            times[4] += utils_stop_section_timer();;
            start=end;
        }

#ifdef REDUCEBCAST
        if (myrank == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, centroid_x_global, clusters, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Reduce(centroid_x_global, centroid_x_global, clusters, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        MPI_Bcast(centroid_x_global, clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (myrank == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, centroid_y_global, clusters, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Reduce(centroid_y_global, centroid_y_global, clusters, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        MPI_Bcast(centroid_y_global, clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#ifdef ALLRED
        MPI_Allreduce(MPI_IN_PLACE, centroid_x_global, clusters, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, centroid_y_global, clusters, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        if (myrank == 0)
        {
            
            //COMM2
            times[5] += utils_stop_section_timer();;
            start=end;
            for (int k = 0; k < clusters; k++)
            {
                sets_counter_global[k] = 0;
                error = error + sets_counter_global[k] + sets_counter_global[k];
            }
    	    
            times[6] += utils_stop_section_timer();;
            start=end;
            int msg;
            if (error == c_error)
            {
                msg = 1;
                MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(sets_local, chunk_size, MPI_INT, sets_global, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
                
                //PHASE33
                times[7] += utils_stop_section_timer();;
                start=end;
                break;
            }
            else
            {
                msg = 0;
                MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
                
                //PHASE33
                times[7] += utils_stop_section_timer();;
                start=end;
            }
        }
        else
        {
            for (int k = 0; k < clusters; k++)
            {
                sets_counter_global[k] = 0;
            }
            int msg;
            MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (msg == 1)
            {
                MPI_Gather(sets_local, chunk_size, MPI_INT, sets_global, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
                break;
            }
        }
    } while (true);

    if (myrank == 0)
    {
        *result = sets_global;
    }

    return;
}
