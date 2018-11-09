#include <utils.h>
#include <kmeansCluster.h>
#define PAR "par"
#define SEQ "seq"
using namespace std;


// kmeans PAR/SEQ REPS SIZE DATASET PAPI_OPT


int main(int argc, char const *argv[])
{
    int repetitions;
    int size;
    char const * filename = NULL;
    char const * mode = NULL;
    char const * papiOpt = NULL;
    double* xcomp = NULL;
    double* ycomp = NULL;
    unsigned char* sets = NULL; 
    if(argc < 5){
        fprintf(stderr, "Usage: ./bin/kmeans PAR|SEQ #REPETITIONS #SIZE DATASET_PATH PAPI_OPT(optional)\n");
        return -1;
    }
    
    mode = argv[1];
    repetitions = atoi(argv[2]);
    size = atoi(argv[3]);
    filename = argv[4];
    papiOpt = argv[5];
    xcomp = (double *) malloc(sizeof(double)* size);
    ycomp = (double *) malloc(sizeof(double)* size);
    int utils_error = utils_read_dataset(filename,xcomp,ycomp);

    if(utils_error == -1)
    {
        return -1;
    }
    
    if( !strcmp(mode, PAR) && argc == 5)
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
    else if( !strcmp(mode, SEQ) && argc == 6)
    {
        utils_setup_papi(repetitions, papiOpt);
        for(int i = 0; i < repetitions; i++ )
        {
            utils_clear_cache();
            utils_start_papi();
            utils_start_timer();
            sets = kmc_seq(3, size, xcomp, ycomp);
            utils_stop_timer();
            utils_stop_papi(i);
        }
    }

    utils_results(papiOpt);
    utils_save_results("bin/kmc_out", xcomp, ycomp, sets);
    utils_clean_memory(xcomp, ycomp); 
    return 0;
}



