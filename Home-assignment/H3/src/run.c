#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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
    int N_max = 3*N0; // maximum number of walkers
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
    int n_eq_steps = 750;
    int n_prod_steps = 3250;
    int n_steps = n_eq_steps + n_prod_steps; // plus the step 0

    int* N_arr = (int*)malloc((n_steps+1)*sizeof(int));
    double* E_T_arr = (double*)malloc((n_steps+1)*sizeof(double));
    double** x_arr = create_2D_array(n_steps+1, N_max);

    N_arr[0] = N;
    E_T_arr[0] = E_T;
    constant_multiplication(x_arr[0], x, 1., N_max);

    double E_T_sum = E_T; // sum all previous E_T for the cummulative average
    // we need a temporary array of walkers in each step
    // because the number of walkers changes
    double* temp_x = (double*)calloc(N_max, sizeof(double));

    // do diffusion monte carlo without importance sampling
    // step through imaginary time 
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
        constant_multiplication(x_arr[i_step], x, 1., N_max);
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
    file = fopen("data/task1_x.csv", "w");
    fprintf(file, "# {\"gamma\":%.5e, \"dtau\":%.5e, \"n_eq_steps\":%i}\n", gamma, dtau, n_eq_steps);
    fprintf(file, "# i_step");
    for (int i_walker=0; i_walker<N_max; i_walker++) {
        fprintf(file, ", x_%06i", i_walker);
    }
    fprintf(file, "\n");
    for (int i_step=0; i_step<n_steps+1; i_step++) {
        fprintf(file, "%8i", i_step);
        for (int i_walker=0; i_walker<N_max; i_walker++) {
            fprintf(file, ", %8.5f", x_arr[i_step][i_walker]);
        }
        fprintf(file, "\n");
    }
    fclose(file);

    free(x);
    free(temp_x);
    free(N_arr);
    free(E_T_arr);
    destroy_2D_array(x_arr);
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

void task2(gsl_rng* rng, double*** R_out, int* N_out, int* N_max_out, double* E_T_out) {
    printf("Perform Task 2...\n");
    // initialize walkers
    int N0 = 1000; // initial number of walkers
    int N_max = 10*N0; // maximum number of walkers
    
    // walker array for each walker: (x1,y1,z1,x2,y2,z2)
    double** R = create_2D_array(N_max, 6);
    int N = N0;
    // initialize N walkers
    for (int i_walker=0; i_walker<N; i_walker++) {
        init_helium_walker(R[i_walker], rng);
    }

    // simulation constants
    double E_T = -3.0;
    double gamma = 0.5;
    double dtau = 0.01;
    double sqrt_dtau = sqrt(dtau);
    int n_eq_steps = 1000;
    int n_prod_steps = 4000;
    int n_steps = n_eq_steps + n_prod_steps; // plus the step 0

    int* N_arr = (int*)malloc((n_steps+1)*sizeof(int));
    double* E_T_arr = (double*)malloc((n_steps+1)*sizeof(double));

    N_arr[0] = N;
    E_T_arr[0] = E_T;

    double E_T_sum = E_T; // sum all previous E_T for the cummulative average
    // we need a temporary array of walkers in each step
    // because the number of walkers changes
    double** temp_R = create_2D_array(N_max, 6);

    // do diffusion monte carlo without importance sampling
    // step through imaginary time 
    for (int i_step=1; i_step<n_steps+1; i_step++) {
        int i_temp = 0;
        for (int i_walker=0; i_walker<N; i_walker++) {
            // diffusive part
            for (int i_coord=0; i_coord<6; i_coord++) {
                R[i_walker][i_coord] = R[i_walker][i_coord] + sqrt_dtau*gsl_ran_gaussian(rng, 1);
            }

            // reactive part
            double W = exp(-(V_helium(R[i_walker]) - E_T)*dtau);
            int m = (int)(W + gsl_rng_uniform(rng));
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

    printf("Save results of Task 2...\n");
    // save results to files
    FILE* file = fopen("data/task2.csv", "w");
    fprintf(file, "# {\"gamma\":%.5e, \"dtau\":%.5e, \"n_eq_steps\":%i}\n", gamma, dtau, n_eq_steps);
    fprintf(file, "# i_step, tau, N_walkers, E_T\n");
    for (int i_step=0; i_step<n_steps+1; i_step++) {
        double tau = i_step*dtau;
        fprintf(file, "%i, %.10f, %i, %.10f\n", i_step, tau, N_arr[i_step], E_T_arr[i_step]);
    }
    fclose(file);

    destroy_2D_array(temp_R);
    free(N_arr);
    free(E_T_arr);
    
    *R_out = R;
    *N_out = N;
    *N_max_out = N_max;
    *E_T_out = E_T;
}

void perform_drift_part(double* R, double alpha, double dtau) {
    double v_F[6] = {0};
    double temp[3];
    double temp2[3];
    double temp3[6];

    double r1[3] = {R[0], R[1], R[2]};
    double r2[3] = {R[3], R[4], R[5]};

    double r12[3];
    constant_multiplication(temp, r1, -1., 3);
    elementwise_addition(r12, r2, temp, 3);

    double r12_len = vector_norm(r12, 3);
    normalize_vector(r12, 3);

    constant_multiplication(temp, r1, -2., 3);
    constant_multiplication(temp2, r12, -1/(2*(1 + 1*alpha*r12_len)*(1 + 1*alpha*r12_len)), 3);
    elementwise_addition(v_F, temp, temp2, 3);

    constant_multiplication(temp, r2, -2., 3);
    constant_multiplication(temp2, r12, +1/(2*(1 + 1*alpha*r12_len)*(1 + 1*alpha*r12_len)), 3);
    elementwise_addition(v_F+3, temp, temp2, 3);

    constant_multiplication(temp3, v_F, dtau, 6);
    elementwise_addition(R, temp3, R, 6);
}

double E_L(double* R, double alpha) {
    double temp[3];
    double temp2[3];

    double r1[3] = {R[0], R[1], R[2]};
    double r2[3] = {R[3], R[4], R[5]};

    double r12_vec[3];
    constant_multiplication(temp, r2, -1., 3);
    elementwise_addition(r12_vec, r1, temp, 3);

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

void task3(gsl_rng* rng, double** R, int N, int N_max, double E_T) {
    printf("Perform Task 3...\n");
    int N0 = 1000;
    // simulation constants
    double alpha = 0.15;
    double gamma = 0.5;
    double dtau = 0.1;
    double sqrt_dtau = sqrt(dtau);
    int n_eq_steps = 0;
    int n_prod_steps = 5000;
    int n_steps = n_eq_steps + n_prod_steps; // plus the step 0

    int* N_arr = (int*)malloc((n_steps+1)*sizeof(int));
    double* E_T_arr = (double*)malloc((n_steps+1)*sizeof(double));

    N_arr[0] = N;
    E_T_arr[0] = E_T;

    double E_T_sum = E_T; // sum all previous E_T for the cummulative average
    // we need a temporary array of walkers in each step
    // because the number of walkers changes
    double** temp_R = create_2D_array(N_max, 6);

    // do diffusion monte carlo without importance sampling
    // step through imaginary time 
    for (int i_step=1; i_step<n_steps+1; i_step++) {
        int i_temp = 0;
        for (int i_walker=0; i_walker<N; i_walker++) {
            // diffusive part
            for (int i_coord=0; i_coord<6; i_coord++) {
                R[i_walker][i_coord] = R[i_walker][i_coord] + sqrt_dtau*gsl_ran_gaussian(rng, 1);
            }

            //drift part
            perform_drift_part(R[i_walker], alpha, dtau);

            // reactive part
            double W = exp(-(E_L(R[i_walker], alpha) - E_T)*dtau);
            int m = (int)(W + gsl_rng_uniform(rng));
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

    printf("Save results of Task 3...\n");
    // save results to files
    FILE* file = fopen("data/task3.csv", "w");
    fprintf(file, "# {\"gamma\":%.5e, \"dtau\":%.5e, \"n_eq_steps\":%i}\n", gamma, dtau, n_eq_steps);
    fprintf(file, "# i_step, tau, N_walkers, E_T\n");
    for (int i_step=0; i_step<n_steps+1; i_step++) {
        double tau = i_step*dtau;
        fprintf(file, "%i, %.10f, %i, %.10f\n", i_step, tau, N_arr[i_step], E_T_arr[i_step]);
    }
    fclose(file);

    destroy_2D_array(temp_R);
    free(N_arr);
    free(E_T_arr);
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

    double** R;
    int N, N_max;
    double E_T;
    task2(rng, &R, &N, &N_max, &E_T);
    task3(rng, R, N, N_max, E_T);
    
    gsl_rng_free(rng);
    destroy_2D_array(R);
    return 0;
}
