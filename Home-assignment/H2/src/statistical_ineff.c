#include <stdio.h>
#include <math.h>
#include "tools.h"

#include "statistical_ineff.h"

double calc_corr_func(double* f, double f2_mean, double f_mean2, int k, int n_f) {
    // f2_mean = <f^2> and f_mean2 = <f>^2
    // calculate <f_{i+k}*f_i> (each product is saved in ff)
    int n_ff = n_f - k;
    double *ff = (double*)malloc(n_f*sizeof(double));
    elementwise_multiplication(ff, f+k, f, n_ff);
    double ff_mean = average(ff, n_ff);
    free(ff);

    return (ff_mean - f_mean2)/(f2_mean - f_mean2);
}

double calc_s_corr(double* f, int n_f) {
    // calculate the statistical inefficiency through the correlation function
    // use binary search to find the k so that the correlation function is at exp(-2)

    // first shift f  
    double f_mean = average(f, n_f);
    double* f_shifted = (double*)malloc(n_f*sizeof(double));
    for (int i=0; i<n_f; i++) {
        f_shifted[i] = f[i] - f_mean;
    }
    // calculate <f^2> and <f>^2
    double* f_squared = (double*)malloc(n_f*sizeof(double));
    elementwise_multiplication(f_squared, f_shifted, f_shifted, n_f);
    double f2_mean = average(f_squared, n_f);
    double f_mean2 = 0;
    free(f_squared);

    // do binary search to find k0 
    // for which phi(k0) is just below exp(-2)
    // and phi(k0-1) is just above exp(-2)
    // assumption: phi(k<k0)>= exp(-2) and phi(k>k0)< exp(-2)
    // (hopefully true)
    double Phi_wanted = exp(-2);
    int k_low = 1;
    int k_high = n_f;
    double Phi_k_low = calc_corr_func(f_shifted, f2_mean, f_mean2, k_low, n_f);
    double Phi_k_high = calc_corr_func(f_shifted, f2_mean, f_mean2, k_high, n_f);
    int n_steps_binary = 0;
    while (k_high-k_low > 1) {
        if (Phi_k_low < Phi_wanted || Phi_k_high >= Phi_wanted) {
            // assumption probably does not hold
            // Phi_wanted is not in the range k_low, k_high anymore
            printf("ERROR: Binary search failed. (Phi_wanted = %.5f)\n", Phi_wanted);
            printf("ERROR: k_low = %i, Phi_k_low = %.5f\n", k_low, Phi_k_low);
            printf("ERROR: k_high = %i, Phi_k_high = %.5f\n", k_high, Phi_k_high);
            exit(1);
        }
        int k_mid = (k_high+k_low)/2;
        double Phi_k_mid = calc_corr_func(f_shifted, f2_mean, f_mean2, k_mid, n_f);
        if (Phi_k_mid >= Phi_wanted) {
            k_low = k_mid;
            Phi_k_low = Phi_k_mid;
        } else {
            k_high = k_mid;
            Phi_k_high = Phi_k_mid;
        }
        //printf("k_low = %i; k_high = %i ; Phi_k_low = %f; Phi_k_high = %f\n",k_low,k_high,Phi_k_low,Phi_k_high);
        n_steps_binary++;
    }
    printf("n_steps_binary_search = %i; s_corr = %i\n", n_steps_binary, k_high);
    return k_high;
}