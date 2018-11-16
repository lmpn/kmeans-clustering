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
<<<<<<< HEAD
    int array[3] = {1048576,2097152,4194304};
=======
    int array[4] = {1048576,2097152,4194304, 33554432};
>>>>>>> 3e035836501dcf345e3b56d2b15115dd5921d5ee

    std::mt19937 rng;
    uint32_t seed_val;
    rng.seed(seed_val);
<<<<<<< HEAD
    long long max = 1000000;
    for(int i = 0; i < 3; i++)
=======
    long long max = 30000000;
    for(int i = 3; i < 4; i++)
>>>>>>> 3e035836501dcf345e3b56d2b15115dd5921d5ee
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
