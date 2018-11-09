#include <kmeansCluster.h>

unsigned char * kmc_seq(int clusters, int size, double *xcomp, double *ycomp)
{
    bool convergence = true;
    unsigned char *sets = NULL;
    unsigned char *temp = (unsigned char *) malloc(sizeof(unsigned char) * size);
    double dist = -DBL_MAX;
    double centroid[clusters*2];
    int p_set = -1;
    do
    {
        convergence = true;
/***************************************************
        Assignment Step
***************************************************/

        for (int p = 0; p < size ;  p++)
        {
            for (int k = 0; k < clusters; k++)
            {
                double xd = xcomp[p] - centroid[k];
                double yd = ycomp[p] - centroid[k+1];
                double k_p_dist = xd*xd + yd*yd;
                dist > k_p_dist ? dist = k_p_dist, p_set = k : 0;
            }
            //fazer verificação 
            temp[p] = p_set;
            p_set = -1;
            dist = -DBL_MAX;
        }
        if(sets == NULL) {
            convergence = false;
            sets = temp;
            temp =(unsigned char *) malloc(sizeof(unsigned char) * size);
        }
        else
        {
            for(int p = 0; p < size && convergence; p++)
            {
                convergence = convergence && (sets[p] != temp[p]);
            }
            if(!convergence)
            {
                for(int k = 0; k < clusters ; k++)
                {
                    centroid[k] = 0.0;
                    centroid[k+1] = 0.0;
                }
                
                for(int i = 0; i < size; i++)
                {
                    int p_set = temp[i];
                    centroid[p_set] += xcomp[i];
                    centroid[p_set+1] += ycomp[i]; 
                }
                for(int k = 0; k < clusters ; k++)
                {
                    centroid[k] /= size;
                    centroid[k+1] /= size;
                } 
                free(sets);
                sets = temp;
                temp = (unsigned char *) malloc(sizeof(unsigned char) * size); 
            }
            else{
                return temp;
            }
            
        }
    }
    while(!convergence);
    return NULL;
}
void kmc_par(){}