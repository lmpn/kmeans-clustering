#include <kmeansCluster.h>
using namespace std;
int * kmc_seq(int clusters, int size, double *xcomp, double *ycomp)
{
    //Random
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
    //double reduction + sets = 0
    for(size_t i = 0; i < size; i++)
    {
        if(max < xcomp[i])
            max = xcomp[i];
        if(max < ycomp[i])
            max = ycomp[i];  
    }
    //sets_counter = 0
    uniform_real_distribution<double> urd_g(0,max); 
    for(int i = 0 ; i < clusters; i++){
        centroid[i*2] = urd_g(rng);
        centroid[1+i*2] = urd_g(rng);
        //sets_counter[i] = 0;
    }

    while( error != 0.0)
    {


        //"ASSIGNMENT STEP"
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




        //DEEP COPY 
        for(int k_idx = 0; k_idx < clusters*2 ; k_idx++)
        {
            centroid_old[k_idx] = centroid[k_idx];
            centroid[k_idx] = 0.0;
        }


        //"UPDATE STEP"
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
        error = norm;
        current_point_cluster_idx      = -1;
        point_set_idx          = 0;
        norm           = 0.0;
        minimun_distance       = DBL_MAX;
    }
    
    return sets;
}
void kmc_par(){}