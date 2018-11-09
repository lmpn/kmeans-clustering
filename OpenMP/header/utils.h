#ifndef UTILS_H 
#define UTILS_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <vector>
#include <string.h>
#include <papi.h>
#define FLOPS "flops"
#define L2MR "l2mr"
#define L3MR "l3mr"
#define TIME_RESOLUTION 1000000 // us
#define MAX_THREADS 48 // maximum number of threads to run

void utils_stop_timer (void); 
void utils_clean_memory(void * xc, void * yc);
void utils_stop_papi(int rep);
void utils_start_papi(void);
void utils_results(char const*);
void utils_setup_papi(int repetitions, char const * type);
int utils_read_dataset(char const * filename, double* xcomp, double* ycomp);
void utils_clear_cache (void); 
void utils_start_timer (void); 
#endif 