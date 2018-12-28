#include <utils.h>
using namespace std;


struct timeval t;
long long unsigned initial_time, partial_time;
double clearcache [30000000];
vector<long long unsigned> *time_measurement = new vector<long long unsigned>(); 
int numEvents;
long long ** values;
#ifdef PAPI
int eventSet = PAPI_NULL;
#endif
int *events;


void utils_start_timer (void) 
{
	gettimeofday(&t, NULL);
	initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}
void utils_start_section_timer (void) 
{
	gettimeofday(&t, NULL);
	partial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}

long long unsigned utils_stop_section_timer (void) 
{
	gettimeofday(&t, NULL);
	long long unsigned final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
	return (final_time - partial_time);
}


long long unsigned start_time (void) 
{
	gettimeofday(&t, NULL);
	long long unsigned partial = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
	return partial;
}

long long unsigned stop_time (long long unsigned partial) 
{
	gettimeofday(&t, NULL);
	long long unsigned final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
	return (final_time - partial);
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
 	FILE *fin;
	fin = fopen(filename, "r");
    if(!fin) {cout<<"Error in read: "<< filename <<endl; return -1; }
	double x = 0.0, y = 0.0;
    while(fscanf(fin,"%lf %lf", &x, &y) != EOF){
		*(xcomp++) = x;
		*(ycomp++) = y;
    }
	return 0;
}




void utils_setup_papi(int repetitions)
{
	#ifdef PAPI
	numEvents = 4;
	events = (int*) malloc(numEvents * sizeof(int));
	events[0] = PAPI_L1_TCM; 
	events[1] = PAPI_L2_TCM;
	events[2] = PAPI_L3_TCM;
	events[3] = PAPI_TOT_INS;
	values = (long long **) malloc(sizeof(long long) * repetitions);	
	for( int i = 0; i < repetitions; i++)
	{
		values[i] = (long long *) malloc(sizeof(long long) * numEvents);
	}
	PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_create_eventset(&eventSet);
 	PAPI_add_events(eventSet, events, numEvents); /* Start the counters */
	#endif
}

void utils_results()
{
	int repetitions = time_measurement->size();
	long long avg1 = 0, avg2 = 0, avg3 = 0;
	for(size_t i = 0; i < repetitions; i++)
	{
		#ifdef PAPI
		cout << values[i][0] << endl;
		cout << values[i][1] << endl;
		cout << values[i][2] << endl;
		cout << values[i][3] << endl;
		#endif
		double tm = time_measurement->at(i) / (double ) 1000;
		cout << "Execution Time #"<< i << ": "<< tm  << "ms"<<  endl;
	}
}


void utils_start_papi()
{
	#ifdef PAPI
		PAPI_start(eventSet);
	#endif	
}


void utils_stop_papi(int rep)
{
	#ifdef PAPI
		PAPI_stop(eventSet, values[rep]);
	#endif
}



void utils_save_results(char const * filename, double *xcomp, double *ycomp, int * sets, int size)
{
    ofstream fout(filename, std::ios::out);
	if( !fout ) 
	{
		cout << "error on write" << endl;
	 	return;
	}
	
	for(size_t i = 0; i < size; i++)
	{
		fout << xcomp[i] << "," << ycomp[i] << "," << sets[i] << endl;
	}
	fout.close();
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
	if(values != NULL)
	{
		for(int i = 0; i < repetitions; i++)
		{
			free(values[i]);
		}
		free(events);
		free(values);
	}
	time_measurement->clear();
	time_measurement->~vector();
}
