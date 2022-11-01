#include <stdio.h>
#include <math.h>

double scalar_product(double* vec1, double* vec2, int length){
    double res = 0;
    for(int i=0;i<length;i++){
        res += vec1[i]*vec2[i];
    }
    return res;
}

double distance(double* vec1, double* vec2, int length){
    double res = 0;
    for(int i=0;i<length;i++){
        res += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
    }
    res = sqrt(res);
    return res;
}

void print_vector(double* vec, int length){
    printf("[");
    for(int i=0; i < length; i++){
        printf("%.2f, ", vec[i]);
    }
    printf("\b\b]\n");
}

void fprint_vector(FILE* file, double* vec, int length){
    for(int i=0; i < length; i++){
        fprintf(file, "%.6f", vec[i]);
        if(i<length-1) fprintf(file, ", ");
    }
    fprintf(file, "\n");
}

void print_matrix(double** mat, int n, int m){
    printf("[");
    for(int i=0; i < n; i++){
        printf("[");
        for(int j=0; j < m; j++){
            printf("%.2f, ", mat[i][j]);
        }
        printf("\b\b],");
        if(i<n-1) printf("\n");
    }
    printf("\b]\n");
}

void fprint_matrix(FILE* file, double** mat, int n, int m){
    // write matrix to a file
    for(int i=0; i < n; i++){
        for(int j=0; j < m; j++){
            fprintf(file, "%.6f", mat[i][j]);
            if(j<m-1) fprintf(file, ", ");
            else fprintf(file, "\n");
        }
    }
}