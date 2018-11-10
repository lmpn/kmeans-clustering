#ifndef KMCLUSTER_H
#define KMCLUSTER_H
#include <random>
#include <float.h>
#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
int * kmc_seq(int clusters, int size, double *xcomp, double *ycomp);
void kmc_par();
#endif