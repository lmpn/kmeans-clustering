#include <kmeansCluster.h>
using namespace std;
int * kmc_seq(int clusters, int size, double *xcomp, double *ycomp)
{
    srand(time(NULL));
    bool convergence = true;
    int *sets = (int *) malloc(sizeof(int) * size);
    int p_set_idx = -1;
    double min_dist = DBL_MAX;
    double centroid[clusters*2]; 
    double centroid_old[clusters*2];
    double size_factor = ((double) 1)/((double) size);
    cout<< size<<endl;
    for(int i = 0 ; i < clusters; i++){
        
        int sc = rand() % size;
        centroid[i*2] = xcomp[sc];
        centroid[i*2+1] = ycomp[sc];
    }

    do
    {
        convergence = true;
        cout << "Assignment Step" << endl;
        for (int p_idx = 0; p_idx < size ;  p_idx++)
        {
            cout << p_idx<<endl;
            for (int p_k_idx = 0; p_k_idx < clusters; p_k_idx++)
            {
                cout << "k:"<< p_k_idx<<endl;
                double xd = xcomp[p_k_idx] - centroid[p_k_idx*2];
                double yd = ycomp[p_k_idx] - centroid[p_k_idx*2+1];
                double k_p_dist = sqrt(pow(xd,2) + pow(yd,2));
                min(min_dist, k_p_dist) == k_p_dist ?
                    min_dist = k_p_dist, p_set_idx = p_k_idx : 
                    0;
            }
            sets[p_idx] = p_set_idx;
            p_set_idx = -1;
            min_dist = DBL_MAX;
        }
        
        cout << "Update Step" << endl;
        for(int k_idx = 0; k_idx < clusters*2 ; k_idx++)
        {
            centroid_old[k_idx] = centroid[k_idx];
            centroid[k_idx] = 0.0;
        }
        
        for(int i = 0; i < size; i++)
        {
            int p_set = sets[i]*clusters;
            centroid[p_set] += xcomp[i]*size_factor;
            centroid[p_set+1] += ycomp[i]*size_factor; 
        }
        double norm = 0.0;
        for(int k = 0; k < clusters*2; k++)
        {
            norm += pow((abs(centroid[k]) - abs(centroid_old[k])),2);
        }
        convergence = sqrt(norm) == 0;
        cout << sqrt(norm)<< endl;
    }
    while(!convergence );
    return sets;
}
void kmc_par(){}


int * kmc_seq2(int clusters, int size, double *xcomp, double *ycomp)
{
    bool convergence = true;
    int *sets = NULL;
    int *temp = (int *) malloc(sizeof(int) * size);
    double dist = DBL_MAX;
    double centroid[clusters*2];
    int p_set = -1;
    do
    {
        convergence = true;
/***************************************************
        Assignment Step
***************************************************/
        for (int p_idx = 0; p_idx < size ;  p_idx++)
        {
            for (int p_k_idx = 0; p_k_idx < clusters; p_k_idx++)
            {
                double xd = xcomp[p_idx] - centroid[p_k_idx];
                double yd = ycomp[p_idx] - centroid[p_k_idx+1];
                double k_p_dist = xd*xd + yd*yd;
                dist > k_p_dist ? dist = k_p_dist, p_set = p_k_idx : 0;
            }
            //fazer verificação 
            temp[p_idx] = p_set;
            p_set = -1;
            dist = -DBL_MAX;
        }
        if(sets == NULL) {
            convergence = false;
            sets = temp;
            temp =(int *) malloc(sizeof(int) * size);
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
                temp = (int *) malloc(sizeof(int) * size); 
            }
            else{
                return temp;
            }
            
        }
    }
    while(!convergence);
    return NULL;
}