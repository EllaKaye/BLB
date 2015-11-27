#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <time.h>

void samp_k_from_n (int *k, int *n, int *a)
{
  int b[*n];

  static gsl_rng *restrict r = NULL;

  if(r == NULL) { // First call to this function, setup RNG
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_mt19937);
  }

  /* sample k objects from an array of size n with replacement */
  for (int i = 0; i < *n; i++)
  {
    b[i] = i;
  }

  gsl_ran_sample (r, a, *k, b, *n, sizeof(int));
}

void bootstrap(double x[], double *result, int *B, int *n) {

  double T_boot[*B];

  int a[*n];

  for (int i = 0; i < *B; i++) {

    samp_k_from_n(n, n, a);

    double x_star[*n];

    for (int j = 0; j < *n; j++)
    {
      x_star[j] = x[a[j]];
    }

    T_boot[i] = gsl_stats_mean(x_star, 1, *n);
  }

  *result = sqrt(gsl_stats_variance(T_boot, 1, *B));
}
