#include <utils.h>

void start (void) {
	gettimeofday(&t, NULL);
	initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}


long long unsigned stop (int thread) {
	gettimeofday(&t, NULL);
	long long unsigned final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;

	if (thread == -1)
		sequential_measurements.push_back(final_time - initial_time);
	else
		measurements[thread].push_back(final_time - initial_time);

	return final_time - initial_time;
}

void clearCache (void) {
	for (unsigned i = 0; i < 30000000; ++i)
		clearcache[i] = i;
}