#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "tools.h"

// morse potential
double V_morse(double x) {
    return 0.5*(1-exp(-x))*(1-exp(-x));
}

void task1(gsl_rng* rng) {
    printf("Perform Task 1...\n");
    // initialize walkers
    int N0 = 200; // initial number of walkers
    int N_max = 10*N0; // maximum number of walkers
    double x_min = -5;
    double x_max = 5;
    double dx = (x_max - x_min)/(N0-1);
    double* x = (double*)calloc(N_max, sizeof(double));
    int N = N0;
    // initialize N walkers
    for (int i=0; i<N; i++) {
        x[i] = x_min + i*dx;
    }
    // fill the rest with NAN
    for (int i=N; i<N_max; i++) {
        x[i] = NAN;
    }

    // simulation constants
    double E_T = 0.5;
    //E_T = 3./8.;
    double gamma = 0.5;
    double dtau = 0.02;
    double sqrt_dtau = sqrt(dtau);
    int n_eq_steps = 10000;
    int n_prod_steps = 100000;
    int n_steps = n_eq_steps + n_prod_steps; // plus the step 0

    // x histogram constants
    x_min = -4;
    x_max = 10;
    int n_bins = 100;
    double x_bin_width = (x_max-x_min)/n_bins;
    int* bin_counts = (int*)calloc(n_bins, sizeof(int));
    int total_bin_count = 0;

    int* N_arr = (int*)malloc((n_steps+1)*sizeof(int));
    double* E_T_arr = (double*)malloc((n_steps+1)*sizeof(double));

    N_arr[0] = N;
    E_T_arr[0] = E_T;

    double E_T_sum = E_T; // sum all previous E_T for the cummulative average
    // we need a temporary array of walkers in each step
    // because the number of walkers changes
    double* temp_x = (double*)calloc(N_max, sizeof(double));

    // do diffusion monte carlo without importance sampling
    // step through imaginary time 
    print_progress(0,0,n_steps, true);
    for (int i_step=1; i_step<n_steps+1; i_step++) {
        int i_temp = 0;
        for (int i_walker=0; i_walker<N; i_walker++) {
            // diffusive part
            x[i_walker] = x[i_walker] + sqrt_dtau*gsl_ran_gaussian(rng, 1);

            // reactive part
            double W = exp(-(V_morse(x[i_walker]) - E_T)*dtau);
            int m = (int)(W + gsl_rng_uniform(rng));
            // copy the walker m times into temp_x for the next step
            for (int i_copy=0; i_copy<m; i_copy++) {
                temp_x[i_temp] = x[i_walker];
                i_temp++;

                if (i_temp >= N_max) {
                    printf("ERROR: there are too many walkers!\n");
                    exit(1);
                }
            }
        }
        N = i_temp;
        // copy temp_x to x so that we have the new list of walkers
        for (int i_walker=0; i_walker<N; i_walker++) {
            x[i_walker] = temp_x[i_walker];
        }
        // fill the rest with NAN for cleanliness
        for (int i_walker=N; i_walker<N_max; i_walker++) {
            x[i_walker] = NAN;
        }

        // adjustment of E_T
        if (i_step < n_eq_steps+1) {
            // in equilibration
            E_T = E_T_sum/i_step - gamma*log(((double)N)/N0);
        } else {
            if (i_step == n_eq_steps+1) {
                // first step of production
                E_T_sum = E_T;
            }
            // in production
            E_T = E_T_sum/(i_step-n_eq_steps) - gamma*log(((double)N)/N0);
        }
        E_T_sum += E_T;
        
        // save the results
        N_arr[i_step] = N;
        E_T_arr[i_step] = E_T;

        // construct the historgram
        if (i_step >= n_eq_steps+1) {
            for (int i_walker=0; i_walker<N; i_walker++) {
                if (x[i_walker]<x_min || x[i_walker]>=x_max) {
                    continue;
                }
                int i_bin = (int)((x[i_walker]-x_min)/x_bin_width);
                bin_counts[i_bin]++;
                total_bin_count++;
            }
        }
        print_progress(i_step,0,n_steps,false);
    }

    printf("Save results of Task 1...\n");
    // save results to files
    FILE* file = fopen("data/task1.csv", "w");
    fprintf(file, "# {\"gamma\":%.5e, \"dtau\":%.5e, \"n_eq_steps\":%i}\n", gamma, dtau, n_eq_steps);
    fprintf(file, "# i_step, tau, N_walkers, E_T\n");
    for (int i_step=0; i_step<n_steps+1; i_step++) {
        double tau = i_step*dtau;
        fprintf(file, "%i, %.10f, %i, %.10f\n", i_step, tau, N_arr[i_step], E_T_arr[i_step]);
    }
    fclose(file);

    // save the histogram
    file = fopen("data/task1_x_hist.csv", "w");
    fprintf(file, "# x_left, x_right, x_center, bin_density\n");
    for (int i_bin=0; i_bin<n_bins; i_bin++) {
        double x_left = x_min + i_bin*x_bin_width;
        double x_right = x_left + x_bin_width;
        double x_center = (x_left+x_right)/2;
        double bin_density = bin_counts[i_bin]/(x_bin_width*total_bin_count);
        fprintf(file, "%.5f, %.5f, %.5f, %.5f\n", x_left, x_right, x_center, bin_density);
    }
    fclose(file);

    free(x);
    free(temp_x);
    free(N_arr);
    free(E_T_arr);
    free(bin_counts);
}

void init_helium_walker(double* pos, gsl_rng* rng) {
    // initialize one walker r as described in task 2
    for (int i=0; i<2; i++) {
        double r = 0.7 + gsl_rng_uniform(rng);
        double theta = acos(2*gsl_rng_uniform(rng) - 1);
        double phi = 2*M_PI*gsl_rng_uniform(rng);
        double x = r*sin(theta)*cos(phi);
        double y = r*sin(theta)*sin(phi);
        double z = r*cos(theta);
        pos[0+i*3] = x;
        pos[1+i*3] = y;
        pos[2+i*3] = z;
    }
}

double V_helium(double* R) {
    double r1 = vector_norm(R,3);
    double r2 = vector_norm(R+3,3);
    double r12 = 0;
    for (int i=0; i<3; i++) {
        r12 += (R[i]-R[i+3])*(R[i]-R[i+3]);
    }
    r12 = sqrt(r12);
    return - 2/r1 - 2/r2 + 1/r12;
}

void calc_v_F(double* R, double* v_F, double alpha) {
    double temp[3];
    double temp2[3];

    double r1[3] = {R[0], R[1], R[2]};
    double r2[3] = {R[3], R[4], R[5]};

    double r12[3];
    constant_multiplication(temp, r1, -1., 3);
    elementwise_addition(r12, r2, temp, 3);

    double r12_len = vector_norm(r12, 3);
    normalize_vector(r12, 3);

    constant_multiplication(temp, r1, -2., 3);
    constant_multiplication(temp2, r12, -1/(2*(1 + alpha*r12_len)*(1 + alpha*r12_len)), 3);
    elementwise_addition(v_F, temp, temp2, 3);

    constant_multiplication(temp, r2, -2., 3);
    constant_multiplication(temp2, r12, +1/(2*(1 + alpha*r12_len)*(1 + alpha*r12_len)), 3);
    elementwise_addition(v_F+3, temp, temp2, 3);
}

void perform_drift_part_1st_order(double* R, double alpha, double dtau) {
    double v_F[6];
    double temp[6];

    calc_v_F(R, v_F, alpha);

    constant_multiplication(temp, v_F, dtau, 6);
    elementwise_addition(R, temp, R, 6);
}

void perform_drift_part_2nd_order(double* R, double alpha, double dtau) {
    double R_half[6];
    double v_F[6];
    double temp[6];

    // calculate R_1/2
    calc_v_F(R, v_F, alpha);
    constant_multiplication(temp, v_F, dtau/2, 6);
    elementwise_addition(R_half, temp, R, 6);

    // calculate R
    calc_v_F(R_half, v_F, alpha);
    constant_multiplication(temp, v_F, dtau/2, 6);
    elementwise_addition(R, temp, R, 6);
}

double E_L(double* R, double alpha) {
    double temp[3];
    double temp2[3];

    double r1[3] = {R[0], R[1], R[2]};
    double r2[3] = {R[3], R[4], R[5]};

    double r12_vec[3];
    constant_multiplication(temp, r1, -1., 3);
    elementwise_addition(r12_vec, r2, temp, 3);

    double r12 = vector_norm(r12_vec, 3);
    normalize_vector(r12_vec, 3);

    constant_multiplication(temp, r2, -1, 3);
    elementwise_addition(temp, r1, temp, 3);

    normalize_vector(r1, 3);
    normalize_vector(r2, 3);
    constant_multiplication(temp2, r2, -1, 3);
    elementwise_addition(temp2, r1, temp2, 3);

    double res = -4;
    res += dot_product(temp,temp2,3)/(r12*(1+alpha*r12)*(1+alpha*r12));
    res -= 1/(r12*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12));
    res -= 1/(4*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12));
    res += 1/r12;

    return res;
}


// perform the diffusion monte carlo for the helium
void perform_diffusion_monte_carlo(gsl_rng* rng, 
                        bool importance_sampling, bool first_order,
                        double** R0, int N0, int N_max,
                        double E_T0, double gamma, double dtau, double alpha,
                        int n_eq_steps, int n_prod_steps,
                        double** R_out, int* N_out, double* E_T_out,
                        int* N_arr, double* E_T_arr) {

    // we need a temporary array of walkers in each step
    // because the number of walkers changes
    double** temp_R = create_2D_array(N_max, 6);
    double** R = create_2D_array(N_max, 6);
    for (int i=0; i<N_max; i++) {
        for (int j=0; j<6; j++) {
            R[i][j] = R0[i][j];
        }
    }

    int N = N0;
    double E_T = E_T0;
    double E_T_sum = E_T; // sum all previous E_T for the cummulative average
    double sqrt_dtau = sqrt(dtau);
    N_arr[0] = N;
    E_T_arr[0] = E_T;

    int n_steps = n_eq_steps + n_prod_steps;
    
    // do diffusion monte carlo, step through imaginary time 
    for (int i_step=1; i_step<n_steps+1; i_step++) {
        int i_temp = 0;
        for (int i_walker=0; i_walker<N; i_walker++) {
            
            // drift part
            if (importance_sampling) {
                if (first_order) {
                    perform_drift_part_1st_order(R[i_walker], alpha, dtau);
                } else {
                    perform_drift_part_2nd_order(R[i_walker], alpha, dtau);
                }
            }

            // diffusive part
            for (int i_coord=0; i_coord<6; i_coord++) {
                R[i_walker][i_coord] = R[i_walker][i_coord] + sqrt_dtau*gsl_ran_gaussian(rng, 1);
            }

            // reactive part
            double W;
            if (importance_sampling) {
                W = exp(-(E_L(R[i_walker], alpha) - E_T)*dtau);
            } else {
                W = exp(-(V_helium(R[i_walker]) - E_T)*dtau);
            }

            int m = (int)(W + gsl_rng_uniform(rng));

            // sometimes N explodes and we want to prevent this
            if (m>N/2.) {
                printf("WARNING: m = %i in step %i for walker %i/%i\n",m,i_step,i_walker, N);
                m = 100;
            }

            // copy the walker m times into temp_x for the next step
            for (int i_copy=0; i_copy<m; i_copy++) {
                for (int i_coord=0; i_coord<6; i_coord++) {
                    temp_R[i_temp][i_coord] = R[i_walker][i_coord];
                }
                i_temp++;

                if (i_temp >= N_max) {
                    printf("ERROR: there are too many walkers!\n");
                    exit(1);
                }
            }
        }
        N = i_temp;
        // copy temp_R to R so that we have the new list of walkers
        for (int i_walker=0; i_walker<N; i_walker++) {
            for (int i_coord=0; i_coord<6; i_coord++) {
                R[i_walker][i_coord] = temp_R[i_walker][i_coord];
            }
        }

        // adjustment of E_T
        if (i_step < n_eq_steps+1) {
            // in equilibration
            E_T = E_T_sum/i_step - gamma*log(((double)N)/N0);
        } else {
            if (i_step == n_eq_steps+1) {
                // first step of production
                E_T_sum = E_T;
            }
            // in production
            E_T = E_T_sum/(i_step-n_eq_steps) - gamma*log(((double)N)/N0);
        }
        E_T_sum += E_T;
        
        // save the results
        N_arr[i_step] = N;
        E_T_arr[i_step] = E_T;
    }

    *N_out = N;
    *E_T_out = E_T;
    for (int i=0; i<N_max; i++) {
        for (int j=0; j<6; j++) {
            R_out[i][j] = R[i][j];
        }
    }

    destroy_2D_array(temp_R);
    destroy_2D_array(R);
}

void save_diffusion_monte_carlo(char* filepath, double gamma, double dtau, 
                                int n_eq_steps, int n_steps,
                                int* N_arr, double* E_T_arr) {
    // save results to files
    FILE* file = fopen(filepath, "w");
    fprintf(file, "# {\"gamma\":%.5e, \"dtau\":%.5e, \"n_eq_steps\":%i}\n", gamma, dtau, n_eq_steps);
    fprintf(file, "# i_step, tau, N_walkers, E_T\n");
    for (int i_step=0; i_step<n_steps+1; i_step++) {
        double tau = i_step*dtau;
        fprintf(file, "%i, %.10f, %i, %.10f\n", i_step, tau, N_arr[i_step], E_T_arr[i_step]);
    }
    fclose(file);
}

int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code
    gsl_rng* rng = init_rng(42);

    task1(rng);

    ///////////////////////////////////////////////
    // task 2
    ///////////////////////////////////////////////
    printf("Task 2...\n");
    int N0 = 1000; // initial number of walkers
    int N_max = 10*N0; // maximum number of walkers
    double E_T0 = -3.0;
    double gamma = 0.5;
    double alpha = 0;
    double dtau = 0.01;
    int n_eq_steps = 10000;
    int n_prod_steps = 40000;
    int n_steps = n_eq_steps+n_prod_steps;
    bool importance_sampling = false;
    bool first_order = true;

    char* filepath = "data/task2.csv";

    // initialize N walkers
    // walker array for each walker: (x1,y1,z1,x2,y2,z2)
    double** R0 = create_2D_array(N_max, 6);
    for (int i_walker=0; i_walker<N0; i_walker++) {
        init_helium_walker(R0[i_walker], rng);
    }

    int* N_arr = (int*)malloc((n_steps+1)*sizeof(int));
    double* E_T_arr = (double*)malloc((n_steps+1)*sizeof(double));

    perform_diffusion_monte_carlo(rng, importance_sampling, first_order, 
                                R0, N0, N_max, E_T0, gamma, dtau, alpha,
                                n_eq_steps, n_prod_steps, 
                                R0, &N0, &E_T0,
                                N_arr, E_T_arr);
    save_diffusion_monte_carlo(filepath, gamma, dtau, 
                            n_eq_steps, n_steps,
                            N_arr, E_T_arr);

    free(N_arr);
    free(E_T_arr);

    ///////////////////////////////////////////////
    // task 3
    ///////////////////////////////////////////////
    printf("Task 3...\n");
    // take N0, N_max, R0, E_T0 from after task 2
    // first order
    printf("-First Order...\n");
    gamma = 0.5;
    alpha = 0.15;
    dtau = 0.1;
    n_eq_steps = 0;
    n_prod_steps = 5000;
    n_steps = n_eq_steps+n_prod_steps;
    importance_sampling = true;
    first_order = true;

    filepath = "data/task3_1st_order.csv";

    N_arr = (int*)malloc((n_steps+1)*sizeof(int));
    E_T_arr = (double*)malloc((n_steps+1)*sizeof(double));

    perform_diffusion_monte_carlo(rng, importance_sampling, first_order, 
                                R0, N0, N_max, E_T0, gamma, dtau, alpha,
                                n_eq_steps, n_prod_steps, 
                                R0, &N0, &E_T0,
                                N_arr, E_T_arr);
    save_diffusion_monte_carlo(filepath, gamma, dtau, 
                            n_eq_steps, n_steps,
                            N_arr, E_T_arr);

    // second order
    printf("-Second Order...\n");
    first_order = false;
    filepath = "data/task3_2nd_order.csv";
    
    perform_diffusion_monte_carlo(rng, importance_sampling, first_order, 
                                R0, N0, N_max, E_T0, gamma, dtau, alpha,
                                n_eq_steps, n_prod_steps, 
                                R0, &N0, &E_T0,
                                N_arr, E_T_arr);
    save_diffusion_monte_carlo(filepath, gamma, dtau, 
                            n_eq_steps, n_steps,
                            N_arr, E_T_arr);


    free(N_arr);
    free(E_T_arr);


    ///////////////////////////////////////////////
    // task 4
    ///////////////////////////////////////////////
    // take N0, N_max, R0, E_T0 from after task 3
    printf("Task 4...\n");
    gamma = 0.5;
    alpha = 0.15;
    importance_sampling = true;

    double dtau_arr[] = {0.05, 0.1, 0.15, 0.2, 0.3, 0.4}; 
    int n_dtau = sizeof(dtau_arr)/sizeof(double);

    double* E_T_1st_order = (double*)malloc(n_dtau*sizeof(double));
    double* E_T_2nd_order = (double*)malloc(n_dtau*sizeof(double));
    double* E_T_std_1st_order = (double*)malloc(n_dtau*sizeof(double));
    double* E_T_std_2nd_order = (double*)malloc(n_dtau*sizeof(double));

    print_progress(0, 0, n_dtau, true);
    for (int i=0; i<n_dtau; i++) {
        dtau = dtau_arr[i];
        n_eq_steps = 100/dtau;
        n_prod_steps = 400/dtau;
        n_steps = n_eq_steps+n_prod_steps;

        N_arr = (int*)malloc((n_steps+1)*sizeof(int));
        E_T_arr = (double*)malloc((n_steps+1)*sizeof(double));

        first_order = true;
        perform_diffusion_monte_carlo(rng, importance_sampling, first_order, 
                                R0, N0, N_max, E_T0, gamma, dtau, alpha,
                                n_eq_steps, n_prod_steps, 
                                R0, &N0, &E_T0,
                                N_arr, E_T_arr);
        E_T_1st_order[i] = average(E_T_arr+n_eq_steps, n_prod_steps);
        E_T_std_1st_order[i] = standard_deviation(E_T_arr+n_eq_steps, n_prod_steps);
        first_order = false;
        perform_diffusion_monte_carlo(rng, importance_sampling, first_order, 
                                R0, N0, N_max, E_T0, gamma, dtau, alpha,
                                n_eq_steps, n_prod_steps, 
                                R0, &N0, &E_T0,
                                N_arr, E_T_arr);
        E_T_2nd_order[i] = average(E_T_arr+n_eq_steps, n_prod_steps);
        E_T_std_2nd_order[i] = standard_deviation(E_T_arr+n_eq_steps, n_prod_steps);

        free(N_arr);
        free(E_T_arr);
        print_progress(i+1,0,n_dtau, false);
    }

    // save task 4 results
    FILE* file = fopen("data/task4.csv", "w");
    fprintf(file, "# dtau, E_T_1st_order, E_T_std_1st_order, E_T_2nd_order, E_T_std_2nd_order\n");
    for (int i=0; i<n_dtau; i++) {
        fprintf(file, "%.5f, %.5f, %.5f, %.5f, %.5f\n", dtau_arr[i], E_T_1st_order[i], E_T_std_1st_order[i], E_T_2nd_order[i], E_T_std_2nd_order[i]);
    }
    fclose(file);
    
    free(E_T_1st_order);
    free(E_T_std_1st_order);
    free(E_T_2nd_order);
    free(E_T_std_2nd_order);
    gsl_rng_free(rng);
    destroy_2D_array(R0);
    return 0;
}
