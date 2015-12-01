//bootstrap for multivariate data that performs regression
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cblas.h>
#include <math.h>
#include <time.h>
#include <omp.h>


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



void bootstrap(double x[], double y[], double* result, int* b, int* B, int *n, int* d) {

	static gsl_rng *restrict r = NULL;

    if(r == NULL) { // First call to this function, setup RNG
        gsl_rng_env_setup();
        r = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(r, time(NULL));
    }


    //a stores the sampled indices
    int a[ *n ];

	//allocate memory for the regression step
	gsl_matrix * pred = gsl_matrix_alloc ( *n, *d );
  	gsl_vector * resp = gsl_vector_alloc( *n );
    gsl_multifit_linear_workspace * work  = gsl_multifit_linear_alloc ( *n, *d );
    gsl_vector* coef = gsl_vector_alloc ( *d );
    gsl_matrix* cov = gsl_matrix_alloc ( *d, *d );
    gsl_matrix * T_boot = gsl_matrix_alloc ( *B, *d );
    double chisq;




    //create bootstrap samples
    for ( int i = 0; i < *B; i++ ) {

        //sample the indices
        samp_k_from_n( n, b, a, r);
        printf("dfdfdfd");

      	//transfer x to a matrix pred and y to a vector resp
  		for ( int i = 0; i < *n; i++ ) {
  			gsl_vector_set (resp, i, y[ a[i] ]);

    		for (int j = 0; j < *d; j++)
      			gsl_matrix_set (pred, i, j, x[ j + ( a[ i ] * (*d) ) ]);
      	}


		//linera regression
      	gsl_multifit_linear ( pred, resp, coef, cov, &chisq,  work );


      	//pass the elements of coef to the ith row of T_boot
      	gsl_matrix_set_row ( T_boot, i, coef );
	}


	//compute the standard deviation of each coefficient accros the bootstrap repetitions
	for ( int j = 0; j < *d; j++){

		result[ j ] = sqrt( gsl_stats_variance( gsl_matrix_ptr ( T_boot, 0, j ), 1, *B ) );
	}

	//free the memory
    gsl_matrix_free (pred);
    gsl_vector_free(resp);
	gsl_multifit_linear_free ( work);
    gsl_vector_free (coef);
    //gsl_vector_free (w);
    gsl_matrix_free (cov);
	printf("\nI AM DONE\n\n");


}
