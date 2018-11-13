#ifndef KMCLUSTER_H
#define KMCLUSTER_H
#include <utils.h>
#include <random>
#include <float.h>
#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#define ZERO 1e-12
int * kmc_seq_final(int clusters, int size, double *xcomp, double *ycomp);
int * kmc_seq_initial(int clusters, int size, double *xcomp, double *ycomp);
void kmc_par();
#endif