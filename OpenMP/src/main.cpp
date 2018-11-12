#include <utils.h>
#include <kmeansCluster.h>

using namespace std;


// kmeans PAR/SEQ REPS SIZE DATASET PAPI_OPT


int main(int argc, char const *argv[])
{
    int repetitions;
    int size;
    int clusters;
    char const * filename = NULL;
    char const * mode = NULL;
    char const * papiOpt = NULL;
    double* xcomp = NULL;
    double* ycomp = NULL;
    int* sets = NULL; 
    if(argc < 6){
        fprintf(stderr, "Usage: ./bin/kmeans PAR|SEQ #REPETITIONS #CLUSTERS #SIZE DATASET_PATH PAPI_OPT(optional)\n");
        return -1;
    }
    #ifdef PAPI
    cout << "PAPI"<<endl;
    #endif
    
    mode = argv[1];
    repetitions = atoi(argv[2]);
    clusters = atoi(argv[3]); 
    size = atoi(argv[4]);
    filename = argv[5];
    papiOpt = argv[6];
    xcomp = (double *) malloc(sizeof(double)* size);
    ycomp = (double *) malloc(sizeof(double)* size);
    int utils_error = utils_read_dataset(filename,xcomp,ycomp);

    if(utils_error == -1)
    {
        return -1;
    }
    
    if( !strcmp(mode, PAR) && argc == 6)
    {
        for(int i = 0; i < repetitions; i++ )
        {
            cout << "it" << endl;
            utils_clear_cache();
            utils_start_timer();
            kmc_par();
            utils_stop_timer();
        }
    }
    else if( !strcmp(mode, SEQ) && argc == 7)
    {
        utils_setup_papi(repetitions, papiOpt);
        for(int i = 0; i < repetitions; i++ )
        {
            utils_start_papi();
            utils_start_timer();
            utils_clear_cache();
            sets = kmc_seq(clusters, size, xcomp, ycomp);
            utils_stop_timer();
            utils_stop_papi(i);
        }
    }

    utils_results(papiOpt);
    utils_save_results("bin/kmc_out.csv", xcomp, ycomp, sets, size);
    utils_clean_memory(xcomp, ycomp); 
    return 0;
}



