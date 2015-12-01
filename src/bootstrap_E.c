#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <time.h>

void samp_k_from_n_E (int *k, int *n, int *a)
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

    samp_k_from_n_E(n, n, a);

    double x_star[*n];

    for (int j = 0; j < *n; j++)
    {
      x_star[j] = x[a[j]];
    }

    T_boot[i] = gsl_stats_mean(x_star, 1, *n);
  }

  *result = sqrt(gsl_stats_variance(T_boot, 1, *B));
}

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

/*
void bootstrap_b_multi(double subsamp[], double *result, int *b, int *B, int *n, gsl_rng *restrict r)
{
  //printf("\nI am now here.\n\n");
  double theta_star[*B]; // to store result of doing statistic (mean) on each replicate

  for (int i = 0; i < *B; i++) {

    // indices for resample
    unsigned int w[*b]; // where the weights (# of each index) will be stored
    double p[*b]; //probs for multinomial resampling, then set equal
    for (int j = 0; j < *b; j++) {
      p[j] = 1;
    }
    gsl_ran_multinomial(r, *b, n, *p, w);

    double v[*b];
    for (int j = 0; j < *b; j++) {
      v[j] = (double) w[j];

      //printf("W[%i] = %i\n", j, w[j]);
      //printf("V[%i] = %f\n", j, v[j]);
      //printf("subsamp[%i] = %f\n", j, subsamp[j]);
    }
    //printf("\n);

    // calculate theta_star
    theta_star[i] = gsl_stats_wmean(v, 1, subsamp, 1, *b);
    //printf("%f\n", theta_star[i]);
    //printf("\n\n T_[%i] = %f\n\n", i, T_boot[i]);
  }


  //printf("\n\n I've finished the for loop\n\n");
  //double temp = sqrt(gsl_stats_variance(T_boot, 1, B));
  //printf("\n\ntemp is %f\n\n", temp);

  *result = sqrt(gsl_stats_variance(theta_star, 1, *B));
  //printf("\n\nresult is %f\n\n", *result);
  //printf("\n\nresult is %f\n\n", sqrt(gsl_stats_variance(T_boot, 1, B)));

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
    bootstrap_b(subsamp, &xis[i], b, *R, *n, r);
  }

  // take average of results from the s subsamples
  *result = gsl_stats_mean(xis, 1, *s);
}

void BLB_serial_multi(double x[], double *result, double *gamma, int *s, int *R, int *n)
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
    double subsamp[*b];
    gsl_ran_choose(r, subsamp, b, x, *n, sizeof(double));

    // run bootstrap on that subsample (uses weights)
    bootstrap_b_multi(subsamp, &xis[i], b, *R, *n, r);
  }

  // take average of results from the s subsamples
  *result = gsl_stats_mean(xis, 1, *s);
}

*/
