/*
Function to perform bootstrap

add flags -lgls -lgslcblas
*/

#include<stdio.h>
#include<gsl/gsl_rng.h>
#include<R.h>
#include "random.h"

//double bootstrap(double data[], size_t numElements, size_t B) {
//
//}

void sample_rng(int *n, int res[]) {

	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	for (size_t i = 0; i < *n; i++) {
		res[i] = gsl_rng_uniform_int(r, *n);
		Rprintf("%i\n", res[i]);
	}

	gsl_rng_free(r);
}


