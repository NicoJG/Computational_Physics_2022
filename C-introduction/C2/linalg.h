#pragma once // This is the include guard
#include <stdio.h>

double scalar_product(double* vec1, double* vec2, int length);
double distance(double* vec1, double* vec2, int length);

void print_vector(double* vec, int length);
void print_matrix(double** mat, int n, int m);
void fprint_vector(FILE* file, double* vec, int length);
void fprint_matrix(FILE* file, double** mat, int n, int m);
