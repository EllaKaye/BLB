#ifndef _TOOLS_H_
#define _TOOLS_H_


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


void samp_k_from_n (int* k, int* n, int* a, gsl_rng *restrict r);

#endif
