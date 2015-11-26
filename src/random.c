#include <stdio.h>
#include <gsl/gsl_rng.h>

int my_random(int *n, int a[]) {

	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	for (int i = 0; i < *n; i++) {
		a[i] = gsl_rng_uniform_int(r, *n);
	}

	gsl_rng_free(r);
	return 0;
}
