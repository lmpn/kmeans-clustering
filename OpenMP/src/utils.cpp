#include <utils.h>
using namespace std;


struct timeval t;
long long unsigned initial_time;
double clearcache [30000000];
vector<long long unsigned> *time_measurement = new vector<long long unsigned>(); 
int numEvents;
long long ** values;
int eventSet = PAPI_NULL;
int *events;








void utils_start_timer (void) 
{
	gettimeofday(&t, NULL);
	initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}


void utils_stop_timer (void) 
{
	gettimeofday(&t, NULL);
	long long unsigned final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
	time_measurement->push_back(final_time - initial_time);
}
void utils_clear_cache (void) 
{
	for (unsigned i = 0; i < 30000000; ++i)
		clearcache[i] = i;
}

int utils_read_dataset(char const * filename, double* xcomp, double* ycomp)
{
	int size = 0;

    ifstream fin(filename, std::ios::in);
    if(!fin) {cout<<"Error in read: "<< filename <<endl; return -1; }
	double x,y;
    while(fin >> x >> y){
		*(xcomp++) = x;
		*(ycomp++) = y;
    }
	return 0;
}




void utils_setup_papi(int repetitions, char const * type)
{
	if (!strcmp(type,FLOPS))
	{
		numEvents = 2;
		events = (int*) malloc(numEvents * sizeof(int));
    	events[0] = PAPI_FP_OPS; 
 		events[1] = PAPI_TOT_CYC; 
	}
	else if (!strcmp(type, L2MR))
	{
		numEvents = 2;
		events = (int*) malloc(numEvents * sizeof(int ));
 		events[0] = PAPI_L2_TCM; 
    	events[1] = PAPI_L1_DCM;
  	}
  	else if (!strcmp(type, L2MR))
  	{
  	  	numEvents = 2;
		events = (int*) malloc(numEvents * sizeof( int));
 		events[0] = PAPI_L3_TCM; 
    	events[1] = PAPI_L2_DCM;
  	}
	values = (long long **) malloc(sizeof(long long) * repetitions);	
	for( int i = 0; i < repetitions; i++)
	{
		values[i] = (long long *) malloc(sizeof(long long) * numEvents);
	}
	PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_create_eventset(&eventSet);
 	PAPI_add_events(eventSet, events, 2); /* Start the counters */
}

void utils_results(char const * type)
{
	int repetitions = time_measurement->size();
	long long avg1 = 0, avg2 = 0, avg3 = 0;
	for(size_t i = 0; i < repetitions; i++)
	{
		avg1 += (*values)[0];
		avg2 += (*values)[1];
		avg3 += time_measurement->at(i);
	}

	if(repetitions != 0)
	{
		cout << "Execution Time:"<<  avg3 / (double) 1000 / (double) repetitions << "ms"<< endl;
	}
	if(type != NULL && !strcmp(type,L3MR) && avg2 != 0)
		{cout << "Level 3 Miss Rate:"<< avg1/avg2 << endl;}
	if(type != NULL && !strcmp(type,L2MR) && avg2 != 0)
		{cout << "Level 2 Miss Rate:"<< avg1/avg2 << endl;}
	else if(type != NULL && !strcmp(type,FLOPS) && avg2 != 0)
		{cout << "FLOPS:"<< avg1/avg2 << endl;}
}

void utils_start_papi()
{
	PAPI_start(eventSet);
}


void utils_stop_papi(int rep)
{
	PAPI_stop(eventSet, values[rep]);
}


void  utils_clean_memory(void * xc, void * yc)
{
	if(xc != NULL)
	{
		free(xc);
	}
	if(yc != NULL)
	{
		free(yc);
	}
	int repetitions = time_measurement->size();
	for(int i = 0; i < repetitions; i++)
	{
		free(values[i]);
	}
	free(events);
	free(values);
	time_measurement->clear();
	time_measurement->~vector();
}