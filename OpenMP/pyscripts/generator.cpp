#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <random>

using namespace std;





int main(void)
{
    int array[3] = {1048576,2097152,4194304};

    std::mt19937 rng;
    uint32_t seed_val;
    rng.seed(seed_val);
    long long max = 1000000;
    for(int i = 0; i < 3; i++)
    {
        uniform_real_distribution<double> urd_g(0,max); 
        std::ofstream fout("input"+to_string(array[i])+".data");
        //#pragma omp parallel for ordered
        for(int j = 0; j < array[i]; j++)
        {
            double x = urd_g(rng);
            double y = urd_g(rng);
            //#pragma omp ordered
            fout <<x << " " << y<< endl; 
        }
        cout << "done"<<endl;
        max *= 2;
    }
    


}
