#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#include "linalg.h"

// execution: runc main.c linalg.c -lm -lgsl -lgslcblas && python plot_hist.py

int main(){
    // read in the length
    int length;
    printf("Enter the vector length: ");
    scanf("%d", &length);
    printf("length = %i\n", length);

    // initialize the vectors (with zeros)
    double* vec1 = (double*)calloc(length, sizeof(double));
    double* vec2 = (double*)calloc(length, sizeof(double));

    // initialize the vectors (with numbers)
    for(int i=0;i<length;i++){
        vec1[i] = i;
        vec2[i] = 2*i;
    }

    // print the vectors
    printf("\nvec1 = ");
    print_vector(vec1,length);
    printf("vec2 = ");
    print_vector(vec2,length);
    
    // calculate and print the scalar product
    printf("vec1*vec2 = %4f\n", scalar_product(vec1, vec2, length));

    // dynamically allocate a nxm array by having pointers on pointers
    int n = length;
    int m = 3;
    double** mat1 = (double**)calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++)
        mat1[i] = (double*)calloc(m, sizeof(double));

    // initialize the matrix
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            mat1[i][j] = i + j*0.1;
        }
    }

    printf("\nmat1 = \n");
    print_matrix(mat1, n, m);

    // calculate all distances
    for(int i=0;i<n;i++){
        for(int j=i;j<n;j++){
            printf("dist(mat1[%d], mat1[%d]) = %.4f\n", i, j, distance(mat1[i], mat1[j], m));
        }
    }

    // save matrix to file
    FILE *file = fopen("mat1.csv", "w+");
    fprint_matrix(file, mat1, n, m);
    fclose(file);

    // generate N random numbers (uniform)
    // https://www.gnu.org/software/gsl/doc/html/rng.html#examples
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);

    int N = 100000;
    double* vec3 = (double*)calloc(N, sizeof(double));
    for (int i=0;i<N;i++){
        vec3[i] = gsl_rng_uniform(rng);
    }
    FILE *file2 = fopen("random_numbers.csv", "w+");
    fprint_vector(file2, vec3, N);
    fclose(file2);

    gsl_rng_free(rng);

    return 0;
}