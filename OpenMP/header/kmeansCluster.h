#ifndef KMCLUSTER_H
#define KMCLUSTER_H
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
unsigned char * kmc_seq(int clusters, int size, double *xcomp, double *ycomp);
void kmc_par();
#endif