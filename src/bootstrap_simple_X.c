// select a random sample of K numbers from the set 0 to 99 without replacement,

#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <stdio.h>
void sample_fun( double data[], int K, int N, int indices[],double new_data[]);


//write a function that accepts a vector x and resamples k out of n

void sample_fun( double data[], int K, int N, int indices[], double new_data[]){
	// the indices from 0 to N-1
	int b[ N ];
	//size of int
	size_t size = sizeof(int) ;
	// the vector of size K of the new indices
	//int indices[ K ];
	// the new sampled vector
	//double newvec [ K ];
	//define r for the random generation
	static gsl_rng *restrict r = NULL;

	if(r == NULL) {
		gsl_rng_env_setup();
		r = gsl_rng_alloc( gsl_rng_mt19937);

	}
	//create the indices

	for (int i = 0; i < N; i++){
   	b[i] =  i;
 	}
 	//sample the new indices
	gsl_ran_sample (r, indices, K, b, N, size);

	/*
	for (int i = 0; i < K; i++ ){
		new_vec [i] =  X[indices[i]];
	printf("mpla mpla ampl %f\n", new_vec[i]);
	}
	*/
	for (int i = 0; i < K; i++){
		new_data[i] = data[indices[i]];
	}

}



int main() {

	// set the seed
	//gsl_rng_set(r,20)


	/*
	//print the vector of the new indices
	for (int i = 0; i < K; i++ ){
	printf("The new indices are %d \n", indices[i]);
	}
	*/
	// create the datavector to be multinormal
	int N = 100;
	double X[ N ];
	int K = 50;
	int index [ K ];
	double my_new_data[K];
	int B = 1000;
	double mean[B];

	//define r for the random generation
	gsl_rng* r = gsl_rng_alloc( gsl_rng_default);


	for (int i =0; i < N; i++){
		X[i] = gsl_ran_gaussian(r,3.5);
	}


	for (int j =0; j < B; j++){

	sample_fun(X,K,N, index, my_new_data);
	mean[j] = gsl_stats_mean(my_new_data, 1, K);

	 //store the newdata in a matrix
	}


	//double gsl_stats_mean (const double data[], size_t stride, size_t n)
	printf("the variance is %f\n",gsl_stats_variance(mean, 1, B));


}
	/*
	// now the new_data contains the new _data
	for (int i = 0; i < K; i++ ){
	printf("%f\n",my_new_data[i]);
*/
