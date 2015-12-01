#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

void bootstrap_b(double x[], double *result, int *b, int *B, int *n, gsl_rng *restrict r);

void bootstrap_b_omp(double x[], double *result, int *b, int *B, int *n, gsl_rng *restrict r)
{
  //printf("\nI am now here.\n\n");
  double T_boot[*B]; // to store result of doing statistic (mean) on each replicate
  
  for (int i = 0; i < *B; i++) {
    
    double x_star[*n];
    gsl_ran_sample(r, x_star, *n, x, *b, sizeof(double));
    
    T_boot[i] = gsl_stats_mean(x_star, 1, *n);
    //printf("\n\n T_[%i] = %f\n\n", i, T_boot[i]);
  }
  
  *result = sqrt(gsl_stats_variance(T_boot, 1, *B));
}

  

void BLB_omp(double x[], double *result, float *gamma, int *s, int *R, int *n)
  // n is length of the data 
{
  static gsl_rng *restrict r = NULL;
  
  if(r == NULL) { // First call to this function, setup RNG
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_mt19937);
  }
  
  int b;
  b = (int) floor(pow(*n, *gamma));
  
  double xi;
  
  #pragma omp parallel for reduction(+:xi)
  for (int i = 0; i < *s; i++)
  {
    double xi_temp;
    gsl_rng_set(r, time(NULL) + i * 1024);
    
    // select subsample of size b of the data
    double subsamp[b];
    gsl_ran_choose(r, subsamp, b, x, *n, sizeof(double));
    
    // run bootstrap on that subsample (resamples n for each replicate)
    bootstrap_b(subsamp, &xi_temp, &b, R, n, r);
    xi += xi_temp;
  }
  
  // take average of results from the s subsamples
  xi = xi / (double) *n;
  *result = xi;	
}
