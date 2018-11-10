#include <kmeansCluster.h>
using namespace std;
int * kmc_seq(int clusters, int size, double *xcomp, double *ycomp)
{
    bool convergence = true;
    int *sets = (int *) malloc(sizeof(int) * size);
    int p_set = -1;
    double dist = DBL_MAX;
    double centroid[3*2] = {7.02904,9.13869,5.72611,5.63526,5.32975,6.69299 };
    double centroid_old[clusters*2];
    double size_factor = ((double) 1)/((double) size);
    /*for(int i = 0 ; i < clusters; i++){
        int sc = rand() % size;
        centroid[i*2] = xcomp[sc];
        centroid[i*2+1] = ycomp[sc];
    }*/

    do
    {
        convergence = true;
/***************************************************
        Assignment Step
***************************************************/
        cout << "Assignment Step" << endl;
        for (int p = 0; p < size ;  p++)
        {
            for (int k = 0; k < clusters; k++)
            {
                double xd = xcomp[p] - centroid[k*2];
                double yd = ycomp[p] - centroid[k*2+1];
                double k_p_dist = sqrt(pow(xd,2) + pow(yd,2));
                printf("%lf\n", k_p_dist);
                min(dist, k_p_dist) == k_p_dist ?
                    dist = k_p_dist, p_set = k : 
                    0;
            }
            printf("%d-%lf\n",p_set, dist);
            sets[p] = p_set;
            p_set = -1;
            dist = DBL_MAX;
        }
        printf("\n");
        
        cout << "Update Step" << endl;
        for(int k = 0; k < clusters ; k++)
        {
            centroid_old[k*2] = centroid[k*2];
            centroid[k*2] = 0.0;
            centroid_old[k*2+1] = centroid[k*2+1];
            centroid[k*2+1] = 0.0;
            //printf("old%d: %.10f %.10f\n", k, centroid_old[k*2], centroid_old[k*2+1]);
        }
        
        for(int i = 0; i < size; i++)
        {
            int p_set = sets[i]*2;
            centroid[p_set] += xcomp[i]*size_factor;
            centroid[p_set+1] += ycomp[i]*size_factor; 
        }
        //cout << "Centroids"<< endl;
        double norm = 0.0;
        for(int k = 0; k < clusters; k++)
        {
            norm += pow((centroid[k*2] - centroid_old[k*2]),2);
            norm += pow((centroid[k*2+1] - centroid_old[k*2+1]),2);
        }
        convergence = sqrt(norm) == 0;
    }
    while(!convergence );
    return sets;
}
void kmc_par(){}