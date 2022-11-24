/* **********************************************
 *
 * Multiply vector v1 by scalar a1
 * results are stored in res
 *
 * **********************************************/
void constant_multiplication(double* res,
							double* v1,
							double a,
							unsigned int len);


/* **********************************************
 *
 * Add v1 and v2 elementwise
 * results is stored in res
 *
 * res should be properly initialised to zero
 * for this function to work correctly
 *
 * **********************************************/
void
elementwise_addition(
		     double *res,
		     double *v1,
		     double *v2,
		     unsigned int len
	            );

/* **********************************************
 *
 * Multiply v1 and v2 elementwise
 * results is stored in res
 *
 * res should be properly initialised to zero
 * for this function to work correctly
 *
 * **********************************************/
void
elementwise_multiplication(
			   double *res,
			   double *v1,
			   double *v2,
			   unsigned int len
	                  );

/* **********************************************
 *
 * Calculate the dot product between
 * v1 and v2
 *
 * the result is returned as a double
 *
 * **********************************************/
double
dot_product(
	    double *v1,
	    double *v2,
	    unsigned int len
	   );

/* **********************************************
 *
 * Allocate the memory to a 2D array
 *
 * **********************************************/
double**
create_2D_array(
		unsigned int n_rows,
		unsigned int n_cols
	       );

/* **********************************************
 *
 * Free the memory to a 2D array created by create_2D_array
 *
 * **********************************************/
void
destroy_2D_array(
		 double **array
	        );

/* **********************************************
 *
 * Calculate the matrix product of the matrices
 * m1 (n x m) and m2 (m x n)
 * so that the resulting matrix is (n x n)
 * and store the result in "result"
 *
 * **********************************************/
void
matrix_multiplication(
		      double **result,
		      double **m1,
		      double **m2,
		      unsigned int m,
		      unsigned int n
		     );

/* **********************************************
 *
 * Calculate the L2 norm of vector v1
 *
 * **********************************************/
double
vector_norm(
	    double *v1,
	    unsigned int len
	   );

/* **********************************************
 *
 * Normalize a given vector v1
 * and store the results in v1
 * v1_normed = v1/|v1|
 *
 * **********************************************/
void
normalize_vector(
		 double *v1,
		 unsigned int len
		);

/* **********************************************
 *
 * Calculate the mean value of vector v1
 * mean = sum(v1)/len(v1)
 *
 * **********************************************/
double
average(
	double *v1,
	unsigned int len
       );

/* **********************************************
 *
 * Calculate the standard deviation value of vector v1
 * I use the following definition:
 * https://numpy.org/doc/stable/reference/generated/numpy.std.html
 * The standard deviation is the square root of the average of the squared deviations from the mean, 
 * i.e., std = sqrt(mean(x)), where x = abs(a - a.mean())**2.
 *
 * **********************************************/
double
standard_deviation(
	    	   double *v1,
	    	   unsigned int len
	          );

/* **********************************************
 *
 * Calculate the euclidian distance between vectors v1 and v2
 *
 * **********************************************/
double
distance_between_vectors(
			 double *v1,
			 double *v2,
			 unsigned int len
	                );
