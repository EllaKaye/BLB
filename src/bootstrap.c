/*
Function for bootstrap
*/

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <math.h>


void samp_k_from_n (int k, int n, int * a)
{
	printf("\nI got further.\n\n");
  int b[n];
  
  static gsl_rng *restrict r = NULL;
  
  if(r == NULL) { // First call to this function, setup RNG
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_mt19937);
  }
	
  /* sample k objects from an array of size n with replacement */
	printf("\n n = %d\n", n);		
  

	printf("\n\nZ\n\n");
	
  for (int i = 0; i < n; i++) 
    {
      b[i] = i;
    }
    
  printf("\n\nA\n\n");
  gsl_ran_sample (r, a, k, b, n, sizeof(int));
  printf("\n\nA\n\n");
}

void bootstrap(double x[], double *result, int B, int n) {
	printf("\nI am now here.\n\n");
	
	double T_boot[B];

	int a[n];
	
	for (int i = 0; i < B; i++) {

		printf("\n i = %d\n", i);		
		samp_k_from_n(n, n, a);
		printf("B\n\n");
		
		double x_star[n];
		
		for (int j = 0; j < n; j++)
		{
		printf("a[%d]=%d\n", j, a[j]);
			x_star[j] = x[a[j]];
		}
				printf("B\n\n");

		T_boot[i] = gsl_stats_mean(x_star, 1, n);
	}
	
	*result = sqrt(gsl_stats_variance(T_boot, 1, B));	
}

int main(void) {
	printf("\nI am here.\n\n");
	double data[6] = {9, 2, 4, 8, 1, 10};
	
	double *result;
	
	bootstrap(data, result, 1000, 6);
	
	printf("se is %f.\n", *result);
	return 0;
}