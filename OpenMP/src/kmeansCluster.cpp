#include <kmeansCluster.h>
using namespace std;
int * kmc_seq(int clusters, int size, double *xcomp, double *ycomp)
{
    srand(time(NULL));
    int p_set_idx         = -1;
    int p_set             = 0;
    int iter              = 0;
    bool convergence      = false;
    double error          = DBL_MAX;
    double norm           = 0.0;
    double min_dist       = DBL_MAX;
    double xd             = 0.0;
    double yd             = 0.0;
    double k_p_dist       = 0.0;
    double size_factor    = ((double) 1)/((double) size);
    int *sets             = (int *) calloc(size,sizeof(int));
    double *centroid      = (double*) malloc(sizeof(double)*clusters*2); 
    double *centroid_old  = (double*) malloc(sizeof(double)*clusters*2);
    double max            = -DBL_MAX;

    
    for(size_t i = 0; i < size; i++)
    {
        if(max < xcomp[i])
            max = xcomp[i];
         if(max < ycomp[i])
            max = ycomp[i];  
    }
    

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uni(0,max); 
    for(int i = 0 ; i < clusters; i++){
        double random_integer = uni(generator);
        centroid[i*2] = random_integer ;
        random_integer = uni(generator);
        centroid[i*2+1] = random_integer;
        cout <<centroid[i*2] << " "<<centroid[i*2+1]<<endl;
    }

    while( convergence != true  && error > 0.0 && iter < 1000000)
    {
        iter ++;
        convergence = true;
        for (int p_idx = 0; p_idx < size ;  p_idx++)
        {
            for (int p_k_idx = 0; p_k_idx < clusters; p_k_idx++)
            {
                xd = xcomp[p_idx] - centroid[p_k_idx*2];
                yd = ycomp[p_idx] - centroid[p_k_idx*2+1];
                k_p_dist = sqrt(pow(xd,2) + pow(yd,2));
                if(min_dist > k_p_dist)
                {
                    min_dist = k_p_dist; 
                    p_set_idx = p_k_idx;
                }
            }
            convergence = convergence && (sets[p_idx] == p_set_idx); 
            sets[p_idx] = p_set_idx;
            p_set_idx = -1;
            min_dist = DBL_MAX;
        }
        for(int k_idx = 0; k_idx < clusters*2 ; k_idx++)
        {
            centroid_old[k_idx] = centroid[k_idx];
            centroid[k_idx] = 0.0;
        }
        
        for(int i = 0; i < size; i++)
        {
            p_set = sets[i]*2;
            centroid[p_set] =centroid[p_set] + xcomp[i]*size_factor;
            centroid[p_set+1] =centroid[p_set+1] + ycomp[i]*size_factor;
        }
        
        for(int k = 0; k < clusters*2; k++)
        {
            norm += pow((abs(centroid[k] - centroid_old[k])),2);
        }
        if ( error > sqrt(norm))
        {
            error = sqrt(norm);
            printf("error change: \n%lf\n\n", error);
        }
        if(iter %1000 == 0){
            printf("iter: %d\n", iter);
            printf("error: %lf\n", error);
        }
        p_set_idx      = -1;
        p_set          = 0;
        norm           = 0.0;
        min_dist       = DBL_MAX;
    }
    
    return sets;
}
void kmc_par(){}