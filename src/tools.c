#include "tools.h"
//samples k out of n elements with replacement and stores them in the vector a
void samp_k_from_n (int* k, int* n, int* a, gsl_rng *restrict r)
{

  //b stores the numbers 0 to -1
  int b[ *n ];

  for (int i = 0; i < (*n) ; i++)
  {
    b[ i ] = i;
  }

  gsl_ran_sample (r, a, *k, b, *n, sizeof( int ));

}
