#ifndef UTILS_H 
#define UTILS_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <vector>
#include <string.h>
#ifdef PAPI
    #include <papi.h>
#endif
#define PAR "par"
#define SEQ "seq"
#define TIME_RESOLUTION 1000000 // us
#define MAX_THREADS 48 // maximum number of threads to run



void utils_start_section_timer (void);
long long unsigned utils_stop_section_timer (void);
void utils_stop_timer (void);
void utils_start_timer (void); 
void utils_clean_memory(void * xc, void * yc);
void utils_stop_papi(int rep);
void utils_start_papi(void);
void utils_setup_papi(int repetitions);
void utils_results();
void utils_save_results(char const * , double *, double *, int * , int);
void utils_clear_cache (void); 
int utils_read_dataset(char const * filename, double* xcomp, double* ycomp);
void printMedian(double*, double*,int);
#endif 
