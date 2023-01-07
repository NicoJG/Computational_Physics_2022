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
    double gamma = 0.5;
    double dtau = 0.02;
    double sqrt_dtau = sqrt(dtau);
    int n_steps = 400; // plus the step 0

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

    print_progress(0,0, n_steps, true);
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
        N = i_temp+1;
        // copy temp_x to x so that we have the new list of walkers
        for (int i_walker=0; i_walker<N; i_walker++) {
            x[i_walker] = temp_x[i_walker];
        }
        // fill the rest with NAN for cleanliness
        for (int i_walker=N; i_walker<N_max; i_walker++) {
            x[i_walker] = NAN;
        }

        // adjustment of E_T
        E_T = E_T_sum/i_step + gamma*log(((double)N)/N0);
        
        // save the results
        N_arr[i_step] = N;
        E_T_arr[i_step] = E_T;
        constant_multiplication(x_arr[i_step], x, 1., N_max);
        
        print_progress(i_step, 0, n_steps, false);
    }

    printf("Save results of Task 1...\n");
    // save results to files
    FILE* file = fopen("data/task1.csv", "w");
    fprintf(file, "# {\"gamma\":%.5e, \"dtau\":%.5e}\n", gamma, dtau);
    fprintf(file, "# i_step, tau, N_walkers, E_T\n");
    for (int i_step=0; i_step<n_steps+1; i_step++) {
        double tau = i_step*dtau;
        fprintf(file, "%i, %.10f, %i, %.10f\n", i_step, tau, N_arr[i_step], E_T_arr[i_step]);
    }
    fclose(file);
    file = fopen("data/task1_x.csv", "w");
    fprintf(file, "# {\"gamma\":%.5e, \"dtau\":%.5e}\n", gamma, dtau);
    fprintf(file, "# i_step");
    for (int i_walker=0; i_walker<N_max; i_walker++) {
        fprintf(file, ", x_%05i", i_walker);
    }
    fprintf(file, "\n");
    for (int i_step=0; i_step<n_steps+1; i_step++) {
        fprintf(file, "%5i", i_step);
        for (int i_walker=0; i_walker<N_max; i_walker++) {
            fprintf(file, ", %.5f", x_arr[i_step][i_walker]);
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
    
    gsl_rng_free(rng);
    return 0;
}
