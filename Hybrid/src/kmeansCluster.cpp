#include <kmeansCluster.h>
using namespace std;
#include <mpi.h>
#include <mm_malloc.h>
#include <omp.h>



void kmc_mpi(int clusters, int size, double *xcomp, double *ycomp, int myrank, int nprocesses, int **result)
{
  double centroid_x_global[clusters];
  double centroid_y_global[clusters];
  int *sets_global = (int *)_mm_malloc(size, sizeof(int));
  double *sets_counter_global = (double *)_mm_malloc(clusters, sizeof(double));
  int chunk_size = size / nprocesses;
  double error = DBL_MAX;
  double c_error = DBL_MAX;
  long long unsigned start;
  int *sets_local = (int *)_mm_malloc(chunk_size, sizeof(int));
  double *xcomp_local = (double *)_mm_malloc(chunk_size * sizeof(double), 64);
  double *ycomp_local = (double *)_mm_malloc(chunk_size * sizeof(double), 64);
  if (myrank == 0) {
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
    	sets_global[i]=0;
    	sets_counter_global[i]=0;
    	sets_local[i]=0;
    }

    uniform_real_distribution<double> urd_g(0, max);

    for (int i = 0; i < clusters; i++)
    {
      centroid_x_global[i] = urd_g(rng);
      centroid_y_global[i] = urd_g(rng);
    }
  }

  MPI_Scatter(xcomp, chunk_size, MPI_DOUBLE, xcomp_local, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(ycomp, chunk_size, MPI_DOUBLE, ycomp_local, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(centroid_x_global, clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(centroid_y_global, clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  do
  {
    #pragma omp parallel for reduction(+: sets_counter_global[: clusters]) num_threads(2)
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

#ifdef REDUCEBCAST
    if (myrank                               == 0)
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
    for (int cluster_idx                  = 0; cluster_idx < clusters; cluster_idx++)
    {
      error                             = error - centroid_y_global[cluster_idx] - centroid_x_global[cluster_idx];
      centroid_x_global[cluster_idx]    = 0.0;
      centroid_y_global[cluster_idx]    = 0.0;
      sets_counter_global[cluster_idx]  = 1 / sets_counter_global[cluster_idx];
    }
    #pragma omp parallel for reduction(+: centroid_x_global[: clusters]) reduction(+: centroid_y_global[: clusters])
    for (int i                                = 0; i < chunk_size; i++)
    {
      int point_set_idx                     = sets_local[i];
      double set_size                       = sets_counter_global[point_set_idx];
      centroid_x_global[point_set_idx]     += xcomp_local[i] * set_size;
      centroid_y_global[point_set_idx]     += ycomp_local[i] * set_size;
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

    for (int k = 0; k < clusters; k++)
    {
      sets_counter_global[k] = 0;
      error = error + sets_counter_global[k] + sets_counter_global[k];
    }

    if (myrank == 0)
    {
      int msg;
      if (error == c_error)
      {
        msg = 1;
        MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(sets_local, chunk_size, MPI_INT, sets_global, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
        break;
      }
      else
      {
        msg = 0;
        MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
    }
    else
    {
      int msg;
      MPI_Bcast(&msg, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (msg == 1)
      {
        MPI_Gather(sets_local, chunk_size, MPI_INT, sets_global, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
        break;
      }
    }
  } while (true);


  return;
}
