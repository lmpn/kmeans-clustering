#include <kmeansCluster.h>
using namespace std;
#include <mpi.h>

int * kmc_seq_final(int clusters, int size, double *xcomp, double *ycomp)
{
    
    std::mt19937 rng;
    uint32_t seed_val;
    rng.seed(seed_val);
    int*    sets                      = (int *) calloc(size,sizeof(int));
    int     point_set_idx             = 0;
    double* sets_counter              = (double *) calloc(clusters,sizeof(double));
    double  max                       = -DBL_MAX;
    double  error                     = DBL_MAX;
    double  c_error = DBL_MAX;

    double centroid_x[clusters];
    double centroid_y[clusters];



    for(size_t i = 0; i < size; i++)
    {
        if(max < xcomp[i])
            max = xcomp[i];
        if(max < ycomp[i])
            max = ycomp[i];  
    }
    uniform_real_distribution<double> urd_g(0,max); 
    for(int i = 0 ; i < clusters; i++){
        centroid_x[i] = urd_g(rng);
        centroid_y[i] = urd_g(rng);
    }


    

    do
    {
        c_error = error;
        error = 0.0;
        
        for (int point_idx = 0; point_idx < size ;  point_idx++)
        {
            int current_point_cluster_idx = -1;
            double minimun_distance = DBL_MAX;

            for (int cluster_idx = 0; cluster_idx < clusters; cluster_idx++)
            {
                double dx = centroid_x[cluster_idx] - xcomp[point_idx];
                double dy = centroid_y[cluster_idx] - ycomp[point_idx];
                double centroid_point_distance = dx*dx + dy*dy;
                if(minimun_distance > centroid_point_distance)
                {
                    minimun_distance = centroid_point_distance; 
                    current_point_cluster_idx = cluster_idx;
                }
            }

            sets_counter[current_point_cluster_idx] += 1.0;
            sets[point_idx] = current_point_cluster_idx;
        }
        




        
        for(int cluster_idx = 0; cluster_idx < clusters ; cluster_idx++)
        {
            error = error - centroid_y[cluster_idx] - centroid_x[cluster_idx];
            centroid_x[cluster_idx] = 0.0;
            centroid_y[cluster_idx] = 0.0;
            sets_counter[cluster_idx] = 1/sets_counter[cluster_idx];
        }
        
        

        for(int i = 0; i < size; i++)
        {
            int point_set_idx = sets[i];
            double set_size = sets_counter[point_set_idx];
            centroid_x[point_set_idx] += xcomp[i]*set_size;
            centroid_y[point_set_idx] += ycomp[i]*set_size;
            
        }
        
        utils_start_section_timer(); 
        for(int k = 0; k < clusters; k++)
        {
            sets_counter[k] = 0;
            error = error + centroid_x[k] + centroid_y[k];
        }
        
    }while(error != c_error);
    
    return sets;
}


int * kmc_seq_initial(int clusters, int size, double *xcomp, double *ycomp)
{
    std::mt19937 rng;
    uint32_t seed_val;
    rng.seed(seed_val);

    int*    sets                      = (int *) calloc(size,sizeof(int));
    int     point_set_idx             = 0;
    int*    sets_counter              = (int *) calloc(clusters,sizeof(int));
    int     current_point_cluster_idx = -1;
    double  max                       = -DBL_MAX;
    double  norm                      = 0.0;
    double  error                     = DBL_MAX;
    double  minimun_distance          = DBL_MAX;
    double  random_real               = 0.0;


    double* centroid                  = (double*) malloc(sizeof(double)*clusters*2); 
    double* centroid_old              = (double*) malloc(sizeof(double)*clusters*2);
    double  centroid_point_distance   = 0.0;
    for(size_t i = 0; i < size; i++)
    {
        if(max < xcomp[i])
            max = xcomp[i];
        if(max < ycomp[i])
            max = ycomp[i];  
    }
    uniform_real_distribution<double> urd_g(0,max); 
    for(int i = 0 ; i < clusters; i++){
        centroid[i*2] = urd_g(rng);
        centroid[1+i*2] = urd_g(rng);
    }

    while( error != norm)
    {
        error = norm;
		norm = 0.0;

        for (int point_idx = 0; point_idx < size ;  point_idx++)
        {
            for (int cluster_idx = 0; cluster_idx < clusters; cluster_idx++)
            {
                double dx = centroid[cluster_idx*2]   - xcomp[point_idx];
                double dy = centroid[cluster_idx*2+1] - ycomp[point_idx];
                double centroid_point_distance = sqrt(pow(dx,2.0) + pow(dy,2.0));
                if(minimun_distance > centroid_point_distance)
                {
                    minimun_distance = centroid_point_distance; 
                    current_point_cluster_idx = cluster_idx;
                }
            }

            sets_counter[current_point_cluster_idx]++;
            sets[point_idx] = current_point_cluster_idx;
            current_point_cluster_idx = -1;
            minimun_distance = DBL_MAX;
        }




        for(int k_idx = 0; k_idx < clusters*2 ; k_idx++)
        {
            centroid_old[k_idx] = centroid[k_idx];
            centroid[k_idx] = 0.0;
        }


        for(int i = 0; i < size; i++)
        {
            int point_set_idx = sets[i]*2;
            centroid[point_set_idx] += xcomp[i];
            centroid[point_set_idx+1] += ycomp[i];
            
        }
        for(int k = 0; k < clusters; k++)
        {
            int set_size = sets_counter[k]; 
            sets_counter[k] = 0;
            centroid[k*2] = centroid[k*2]/ set_size;
            centroid[k*2+1] = centroid[k*2+1]/ set_size;
            norm += centroid[k*2]-centroid_old[k*2];
            norm += centroid[k*2+1]-centroid_old[k*2+1];
        }
        current_point_cluster_idx      = -1;
        point_set_idx          = 0;
        minimun_distance       = DBL_MAX;
    }
    
    return sets;
}




int * kmc_mpi(int clusters, int size, double * xcomp, double * ycomp) {

  
  
  
  std::mt19937 rng;
  uint32_t seed_val;
  rng.seed(seed_val);
  int * sets = (int * ) calloc(size, sizeof(int));
  int point_set_idx = 0;
  double * sets_counter = (double * ) calloc(clusters, sizeof(double));
  double max = -DBL_MAX;
  double error = DBL_MAX;
  double c_error = DBL_MAX;

  double centroid_x[clusters];
  double centroid_y[clusters];
  

  
  for (size_t i = 0; i < size; i++) {
    if (max < xcomp[i])
      max = xcomp[i];
    if (max < ycomp[i])
      max = ycomp[i];
  }
  
  uniform_real_distribution < double > urd_g(0, max);
  
  for (int i = 0; i < clusters; i++) {
    centroid_x[i] = urd_g(rng);
    centroid_y[i] = urd_g(rng);
  }
  

  do

  {
    c_error = error;
    error = 0.0;
    
    for (int point_idx = 0; point_idx < size; point_idx++) {
      int current_point_cluster_idx = -1;
      double minimun_distance = DBL_MAX;

      for (int cluster_idx = 0; cluster_idx < clusters; cluster_idx++) {
        double dx = centroid_x[cluster_idx] - xcomp[point_idx];
        double dy = centroid_y[cluster_idx] - ycomp[point_idx];
        double centroid_point_distance = dx * dx + dy * dy;
        if (minimun_distance > centroid_point_distance) {
          minimun_distance = centroid_point_distance;
          current_point_cluster_idx = cluster_idx;
        }
      }

      sets_counter[current_point_cluster_idx] += 1.0;
      sets[point_idx] = current_point_cluster_idx;
    }

    
    

    for (int cluster_idx = 0; cluster_idx < clusters; cluster_idx++) {
      error = error - centroid_y[cluster_idx] - centroid_x[cluster_idx];
      centroid_x[cluster_idx] = 0.0;
      centroid_y[cluster_idx] = 0.0;
      sets_counter[cluster_idx] = 1 / sets_counter[cluster_idx];
    }
    
    

    
    for (int i = 0; i < size; i++) {
      int point_set_idx = sets[i];
      double set_size = sets_counter[point_set_idx];
      centroid_x[point_set_idx] += xcomp[i] * set_size;
      centroid_y[point_set_idx] += ycomp[i] * set_size;

    }

    
    

    for (int k = 0; k < clusters; k++) {
      sets_counter[k] = 0;
      error = error + centroid_x[k] + centroid_y[k];
    }
    

  } while (error != c_error);

return sets;

}


