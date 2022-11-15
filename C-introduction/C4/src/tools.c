#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"

void constant_multiplication(double* res,
							double* v1,
							double a1,
							unsigned int len)
{
	for(int i=0;i<len;i++){
		res[i] = a1*v1[i];
	}
}
    
void 
elementwise_addition(
		     double *res,
		     double *v1,
		     double *v2,
		     unsigned int len
	            )
{
	for(int i=0;i<len;i++){
		res[i] = v1[i] + v2[i];
	}
}

void
elementwise_multiplication(
			   double *res,
			   double *v1,
			   double *v2,
			   unsigned int len
	                  )
{
	for(int i=0;i<len;i++){
		res[i] = v1[i] * v2[i];
	}
}

double
dot_product(
	    double *v1,
	    double *v2,
	    unsigned int len
	   )
{
	double res = 0;
	for(int i=0;i<len;i++){
		res += v1[i] * v2[i];
	}
    return res;
}

void
create_2D_array(
		double ***array,
		unsigned int nrows,
		unsigned int ncols
	       )
{
	// allocate 1D array of doubles containing the whole matrix
	double* linear_array = (double*)malloc(nrows*ncols*sizeof(double));
	// allocate 1D array of pointers to doubles containing the pointers to each row starting point
	*array = (double**)malloc(nrows*sizeof(double*));
	// let each row pointer point to the correct address
	for(int row=0;row<nrows;row++){
		(*array)[row] = linear_array + row*ncols;
	}
}

void
destroy_2D_array(
		 double **array
		)
{
	// free the linear_array
	free(array[0]);
	// free the pointers array
	free(array);
}

void
matrix_multiplication(
		      double **result,
		      double **m1,
		      double **m2,
		      unsigned int m,
		      unsigned int n
		     )
{
	// https://en.wikipedia.org/wiki/Computational_complexity_of_matrix_multiplication#Schoolbook_algorithm
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			// calculate matrix element in row i, col j
			result[i][j] = 0;
			for(int k=0;k<m;k++){
				result[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
}

double
vector_norm(
	    double *v1,
	    unsigned int len
	   )
{
    return sqrt(dot_product(v1, v1, len));
}


void
normalize_vector(
		 double *v1,
		 unsigned int len
		)
{
	double norm = vector_norm(v1, len);
	constant_multiplication(v1, v1, 1./norm, len);
}

double
average(
	double *v1,
	unsigned int len
       )
{
	double res = 0;
	for(int i=0;i<len;i++){
		res += v1[i];
	}
    return res/len;
}


double
standard_deviation(
	    	   double *v1,
	    	   unsigned int len
	          )
{
	/* https://numpy.org/doc/stable/reference/generated/numpy.std.html
	 * The standard deviation is the square root of the average of the squared deviations from the mean, 
	 * i.e., std = sqrt(mean(x)), where x = abs(a - a.mean())**2. 
	 * std(v1) = sqrt(sum(v1 - v1_mean)^2 / len(v1)) */
	double mean = average(v1, len);
	double res = 0;
	for(int i=0;i<len;i++){
		res += (v1[i] - mean)*(v1[i] - mean);
	}
    return sqrt(res/len);
}

double
distance_between_vectors(
			 double *v1,
			 double *v2,
			 unsigned int len
	                )
{
	// dist(v1,v2) = |v1 - v2|
	double res = 0;
	for(int i=0;i<len;i++){
		res += (v1[i] - v2[i])*(v1[i] - v2[i]);
	}
    return sqrt(res);
}
