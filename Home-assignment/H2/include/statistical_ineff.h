#pragma once

double calc_corr_func(double* f, double f2_mean, double f_mean2, int k, int n_f);
void calc_all_corr_func(double* f, int n_f, int* k, double* Phi, int n_k);
double calc_s_corr(double* f, int n_f);