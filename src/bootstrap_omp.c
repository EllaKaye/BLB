#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

void bootstrap_b(double x[], double *result, int *b, int *B, int *n, gsl_rng *restrict r)
{
  //printf("\nI am now here.\n\n");
  double T_boot[*B]; // to store result of doing statistic (mean) on each replicate
  
  //int a[n];
  
  for (int i = 0; i < *B; i++) {
    
    //printf("\n i = %d\n", i);
    //samp_k_from_n(n, n, a, r);
    //printf("B\n\n");
    
    double x_star[*n];
    gsl_ran_sample(r, x_star, *n, x, *b, sizeof(double));
    
    T_boot[i] = gsl_stats_mean(x_star, 1, *n);
    //printf("\n\n T_[%i] = %f\n\n", i, T_boot[i]);
  }
  
  
  //printf("\n\n I've finished the for loop\n\n");
  //double temp = sqrt(gsl_stats_variance(T_boot, 1, B));
  //printf("\n\ntemp is %f\n\n", temp);
  
  *result = sqrt(gsl_stats_variance(T_boot, 1, *B));
  //printf("\n\nresult is %f\n\n", *result);
  //printf("\n\nresult is %f\n\n", sqrt(gsl_stats_variance(T_boot, 1, B)));
  
}

void bootstrap_b_omp(double x[], double *result, int *b, int *B, int *n)
{
  //printf("\nI am now here.\n\n");
  // double T_boot[*B]; // to store result of doing statistic (mean) on each replicate
  
  double theta_star = 0.0;
  double theta_star_sq = 0.0;
  
  #pragma omp parallel reduction(+:theta_star, theta_star_sq)
  {
    gsl_rng *restrict r = NULL;
    if(r == NULL) { // First call to this function, setup RNG
      gsl_rng_env_setup();
      r = gsl_rng_alloc(gsl_rng_mt19937);
    }
    gsl_rng_set(r, time(NULL) + omp_get_thread_num() * 1024);
    
    #pragma omp for 
    for (int i = 0; i < *B; i++) {
    
      double theta_star_temp;
      double theta_star_sq_temp;
    
      double x_star[*n];
      gsl_ran_sample(r, x_star, *n, x, *b, sizeof(double));
    
      theta_star_temp = gsl_stats_mean(x_star, 1, *n);
      theta_star_sq_temp = pow(theta_star_temp, 2);
      theta_star += theta_star_temp;
      theta_star_sq += theta_star_sq_temp;
    //printf("\n\n T_[%i] = %f\n\n", i, T_boot[i]);
    }
  }
  double theta_star_bar = theta_star / (double) *B;
  double var;
  var = (theta_star_sq - *B * pow(theta_star_bar, 2)) / (*B - 1);
  *result = sqrt(var);
}

void BLB_serial(double x[], double *result, double *gamma, int *s, int *R, int *n)
  // n is length of the data
{
  static gsl_rng *restrict r = NULL;
  
  if(r == NULL) { // First call to this function, setup RNG
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, time(NULL));
  }
  
  int b;
  b = (int) floor(pow(*n, *gamma));
  
  double xis[*s];
  
  for (int i = 0; i < *s; i++)
  {
    // select subsample of size b of the data
    double subsamp[b];
    gsl_ran_choose(r, subsamp, b, x, *n, sizeof(double));
    
    // run bootstrap on that subsample (resamples n for each replicate)
    bootstrap_b(subsamp, &xis[i], &b, R, n, r);
  }
  
  // take average of results from the s subsamples
  *result = gsl_stats_mean(xis, 1, *s);
}

void BLB_omp_on_s(double x[], double *result, double *gamma, int *s, int *R, int *n)
  // n is length of the data 
{
  int b;
  b = (int) floor(pow(*n, *gamma));
  
  double xi = 0;
  
  #pragma omp parallel reduction(+:xi)
  {
    gsl_rng *restrict r = NULL;
    if(r == NULL) { // First call to this function, setup RNG
      gsl_rng_env_setup();
      r = gsl_rng_alloc(gsl_rng_mt19937);
    }
    gsl_rng_set(r, time(NULL) + omp_get_thread_num() * 1024);
    
    #pragma omp for
    for (int i = 0; i < *s; i++)
    {
      double xi_temp;
      
      // select subsample of size b of the data
      double subsamp[b];
      gsl_ran_choose(r, subsamp, b, x, *n, sizeof(double));
      
      // run bootstrap on that subsample (resamples n for each replicate)
      bootstrap_b(subsamp, &xi_temp, &b, R, n, r);
      xi += xi_temp;
    }
  }
  // take average of results from the s subsamples
  xi = xi / (double) *s;
  *result = xi;	
}

void BLB_omp_on_B(double x[], double *result, double *gamma, int *s, int *R, int *n)
  // n is length of the data
{
  gsl_rng *restrict r = NULL;
  
  if(r == NULL) { // First call to this function, setup RNG
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, time(NULL));
  }
  
  int b;
  b = (int) floor(pow(*n, *gamma));
  
  double xis[*s];
  
  for (int i = 0; i < *s; i++)
  {
    // select subsample of size b of the data
    double subsamp[b];
    gsl_ran_choose(r, subsamp, b, x, *n, sizeof(double));
    
    // run bootstrap on that subsample (resamples n for each replicate)
    bootstrap_b_omp(subsamp, &xis[i], &b, R, n);
  }
  
  // take average of results from the s subsamples
  *result = gsl_stats_mean(xis, 1, *s);
}
