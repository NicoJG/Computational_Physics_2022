#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "tools.h"

gsl_rng* init_rng(int seed){
    // seed = 0 means the seed is random
    // set up the random number generator from GSL
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);
    if(seed != 0){
        gsl_rng_set(rng, seed); 
    }

    return rng;
}

/*****************************************************
 * Functions for task 1 and 2
 *****************************************************/
struct mc_result{
    double integral, variance;
    double* x;
};
struct mc_result mc_integrate_uniform(double (*func)(double), int N, gsl_rng* rng){
    // first argument is the function inside the integral
    //https://stackoverflow.com/questions/9410/how-do-you-pass-a-function-as-a-parameter-in-c

    struct mc_result res;

    res.x = (double*)malloc(N*sizeof(double));

    // perform the Monte Carlo Integration
    double I_MC = 0; // Monte Carlo Intergral
    double var_MC = 0; // Monte Carlo Integral variance
    double x, f;
    for(int i=0;i<N;i++){
        x = gsl_rng_uniform(rng);
        res.x[i] = x;
        f = (*func)(x);
        I_MC += f;
        var_MC += f*f;
    }
    I_MC /= N;
    var_MC = var_MC/N - I_MC*I_MC;
    var_MC /= N;

    res.integral = I_MC;
    res.variance = var_MC;

    return res;
}

struct mc_result mc_integrate_importance_sampling(double (*func)(double), double (*pdf)(double), double (*cdf_inv)(double), int N, gsl_rng* rng) {
    // first arguments are functions
    // https://stackoverflow.com/questions/9410/how-do-you-pass-a-function-as-a-parameter-in-c
    // func is f(x)
    // pdf is PDF(x)
    // cdf_inv is CDF^-1(y)

    struct mc_result res;
    res.x = (double*)malloc(N*sizeof(double));

    // perform the Monte Carlo Integration
    double I_MC = 0; // Monte Carlo Intergral
    double var_MC = 0; // Monte Carlo Integral variance
    double u, x, g;
    for (int i=0; i<N; i++) {
        u = gsl_rng_uniform(rng);
        x = (*cdf_inv)(u);
        res.x[i] = x;
        g = (*func)(x) / (*pdf)(x);
        I_MC += g;
        var_MC += g*g;
    }
    I_MC /= N;
    var_MC = var_MC/N - I_MC*I_MC;
    var_MC /= N;

    res.integral = I_MC;
    res.variance = var_MC;

    return res;
}

double f(double x) {
    return x*(1-x);
}

double p(double x) {
    if (0<=x && x<=1) {
        return M_PI/2 * sin(M_PI*x);
    } else {
        return 0.;
    }
}

double F_inv(double y) {
    return acos(1-2*y) / M_PI;
}


/*****************************************************
 * Functions for task 3
 *****************************************************/
void metropolis_algorithm(int n_dim, int N_steps, int N_burn, double init_pos[n_dim],
                        double** samples, double** samples_burn,
                        double (*pdf)(int, double*), void (*new_trial_pos)(int, double*, double*, gsl_rng*),
                        gsl_rng* rng) {
    // Generate N_steps samples (n_dim numbers per sample) using the Metropolis algorithm
    // after a burn in period of N_burn samples
    // given are also initial coordinates
    // the samples are saved in the arrays "samples" and "samples_burn"
    // pdf and new_trial_pos are pointers to functions

    int accept_count = 0;

    // setup
    double trial[n_dim]; // trial position
    double curr[n_dim]; // current position
    for (int i=0; i<n_dim; i++) {
        curr[i] = init_pos[i];
        samples_burn[0][i] = init_pos[i];
    }

    double w_curr = (*pdf)(n_dim, curr);
    double w_trial;

    // sampling
    for (int i=1; i<(N_burn+N_steps); i++) {
        (*new_trial_pos)(n_dim, curr, trial, rng);
        w_trial = pdf(n_dim, trial);
        if (w_trial/w_curr >= gsl_rng_uniform(rng)) {
            // accept new position
            for (int j=0; j<n_dim; j++) {
                curr[j] = trial[j];
            }
            w_curr = w_trial;
            accept_count++;
        }

        // save sample
        for (int j=0; j<n_dim; j++) {
            if (i < N_burn) {
                samples_burn[i][j] = curr[j];
            } else {
                samples[i-N_burn][j] = curr[j];
            }
        }
    }
    printf("acceptance rate: %.5f\n", ((double)accept_count)/(N_burn+N_steps));
}

double weightfunc(int n_dim, double* pos){
    // in this case n_dim must be 3
    if (n_dim != 3) {
        perror("n_dim != 3 in weightfunc!");
        exit(1);
    }

    // pi^(-3/2) https://www.wolframalpha.com/input?key=&i=pi%5E%28-3%2F2%29
    double pi_factor = 0.1795871221251665616890819836276927552821859434348267246213320371;

    double res = 0;
    for (int i=0; i<n_dim; i++) {
        res += pos[i]*pos[i];
    }
    res = pi_factor * exp(-res);
    return res;
}

void get_next_trial(int n_dims, double* curr, double* trial, gsl_rng* rng) {
    double delta = 2;
    for (int i=0; i<n_dims; i++) {
        trial[i] = curr[i] + delta * (gsl_rng_uniform(rng) - 0.5);
    }
}

double f_task3(int n_dim, double* x){
    return x[0]*x[0] * (1 + x[1]*x[1] * (1 + x[2]*x[2]));
}


/*****************************************************
 * Main function
 *****************************************************/
int main() {
    gsl_rng* rng = init_rng(42);

    // Task 1:
    int N[] = {1, 3, 10, 30, 1e2, 3e2, 1e3, 3e3, 1e4};
    int len_N = sizeof(N) / sizeof(int);
    struct mc_result res;
    FILE* file = fopen("data/E3_1.csv", "w");
    fprintf(file, "# N, integral, variance\n");
    for (int i=0; i<len_N; i++) {
        res = mc_integrate_uniform(f, N[i], rng);
        fprintf(file, "%5i, %10.10f, %10.10f \n", N[i], res.integral, res.variance);
        free(res.x);
    }
    fclose(file);

    // Task 2:
    file = fopen("data/E3_2.csv", "w");
    fprintf(file, "# N, integral, variance\n");
    for (int i=0; i<len_N; i++) {
        res = mc_integrate_importance_sampling(f, p, F_inv, N[i], rng);
        fprintf(file, "%5i, %10.10f, %10.10f \n", N[i], res.integral, res.variance);
        if (i < len_N-1) {
            free(res.x);
        }
    }
    fclose(file);
    // save the x for the last one example
    file = fopen("data/E3_2_x.csv", "w");
    fprintf(file, "# x, f(x), p(x)\n");
    for (int i=0; i<N[len_N-1]; i++) {
        fprintf(file, "%.10f, %.10f, %.10f\n", res.x[i], f(res.x[i]), p(res.x[i]));
    }
    fclose(file);
    free(res.x);


    // Task 3:
    int n_dim = 3;
    int N_steps = 10000000;
    int N_burn = 10000;
    int N_steps_save = 100000;
    double** samples = create_2D_array(N_steps,n_dim);
    double** samples_burn = create_2D_array(N_burn,n_dim);
    double init_pos[3] = {2, 2, 2};
    metropolis_algorithm(n_dim, N_steps, N_burn, init_pos, samples, samples_burn,
                        weightfunc, get_next_trial, rng);
    
    // save the samples
    file = fopen("data/E3_3_samples.csv", "w");
    fprintf(file, "# {\"n_dim\": %i, \"N_burn\": %i, \"N_steps\": %i}\n# ", n_dim, N_burn, N_steps);
    fprintf(file, "x_0");
    for (int i=1; i<n_dim; i++) {
        fprintf(file, ", x_%i", i);
    }
    fprintf(file, "\n");
    // burn samples
    for (int i=0; i<N_burn; i++) {
        fprintf(file, "%.10f", samples_burn[i][0]);
        for (int j=1; j<n_dim; j++) {
            fprintf(file, ", %.10f", samples_burn[i][j]);
        }
        fprintf(file, "\n");
    }
    // actual samples (save not all of them)
    fprintf(file, "# actual samples begin:\n");
    for (int i=0; i<N_steps_save; i++) {
        fprintf(file, "%.10f", samples[i][0]);
        for (int j=1; j<n_dim; j++) {
            fprintf(file, ", %.10f", samples[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);

    double I = 0;
    for (int i=0; i<N_steps; i++) {
        I += f_task3(n_dim, samples[i]);
    }
    I /= N_steps;
    printf("Task 3 integral result:\n%.10f\n",I);

    destroy_2D_array(samples);
    destroy_2D_array(samples_burn);
    // free the allocated random number generator
    gsl_rng_free(rng);

    return 0;
}