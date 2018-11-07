#ifndef UTILS_H 
#define UTILS_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <vector>
#include <string>

#define TIME_RESOLUTION 1000000 // us
#define MAX_THREADS 48 // maximum number of threads to run

using namespace std;
struct timeval t;
long long unsigned initial_time;
double clearcache [30000000];
vector<long long unsigned> sequential_measurements;
vector<long long unsigned> measurements[MAX_THREADS/4 + 1];	// position 0 is for the sequential version
#endif 