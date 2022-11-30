#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "tools.h"
#include "lattice_tools.h"

// physical constants
const double k_B = 8.617333262e-5; // eV/K;

// parameters of the simulation
const double E_AA = -436e-3; // eV
const double E_BB = -113e-3; // eV
const double E_AB = -294e-3; // eV
const int n_cells = 10;

double calc_energy(int N_atoms, int* atype, int** nn_idxs) {
    // count the types of neares neighbor pairs
    int N_AA = 0;
    int N_BB = 0;
    int N_AB = 0;

    for (int i=0; i<N_atoms; i++) {
        int atype_a = atype[i];
        for (int j=0; j<8; j++) {
            int atype_b = atype[nn_idxs[i][j]];
            if (atype_a == 0 && atype_b == 0) {
                N_AA++;
            } else if (atype_a == 1 && atype_b == 1) {
                N_BB++;
            } else {
                N_AB++;
            }
        }
    }
    // but now we counted each bond twice
    if (N_AA%2 == 1 || N_BB%2 == 1 || N_AB%2 == 1) {
        perror("ERROR: number of pairs should be even");
        exit(1);
    }
    N_AA /= 2;
    N_BB /= 2;
    N_AB /= 2;

    return N_AA*E_AA + N_BB*E_BB + N_AB*E_AB;
}

double calc_P(int N_atoms, int* atype, int** pos) {
    // long-range order parameter P
    int N_Aa = 0; // number of A atoms on sublattice a
    int N_A = 0; // number of A atoms
    for (int i=0; i<N_atoms; i++) {
        if (atype[i]==0) {
            N_A++;
            if (pos[i][3]==0) {
                N_Aa++;
            }
        }
    }

    return 2*((double)N_Aa)/N_A - 1;
}

double calc_r(int N_atoms, int* atype, int** nn_idxs) {
    // short-range order parameter r
    int N_AB = 0;
    for (int i=0; i<N_atoms; i++) {
        int atype_a = atype[i];
        for (int j=0; j<8; j++) {
            int atype_b = atype[nn_idxs[i][j]];
            if (atype_a != atype_b) {
                N_AB++;
            }
        }
    }
    // we counted each bond twice
    if (N_AB%2 == 1) {
        perror("ERROR: number of pairs should be even");
        exit(1);
    }
    N_AB /= 2;

    return (N_AB-4)/4.;
}

void metropolis_algorithm(int N_atoms, int n_steps, double T, 
                        int* atype, int** pos, int** nn_idxs, gsl_rng* rng, 
                        double* E_out, double* P_out, double* r_out) {
    double E_curr = calc_energy(N_atoms, atype, nn_idxs);
    print_progress(0,0,n_steps-1,true);
    // calc the energy of the current configuration
    for (int i_step=0; i_step<n_steps; i_step++) {
        // make trial change by swapping two atoms
        // choose which atoms to swap
        int swap_idx_a = gsl_rng_uniform_int(rng, N_atoms);
        int swap_idx_b = gsl_rng_uniform_int(rng, N_atoms);
        while (swap_idx_a == swap_idx_b) {
            swap_idx_b = gsl_rng_uniform_int(rng, N_atoms);
        }

        // perform the swap
        int temp_atype = atype[swap_idx_a];
        atype[swap_idx_a] = atype[swap_idx_b];
        atype[swap_idx_b] = temp_atype;

        // check if the swap should be accepted as the new configuration
        double E_trial = calc_energy(N_atoms, atype, nn_idxs);
        double DeltaE = E_trial - E_curr;
        bool accept = (DeltaE <= 0);
        if (!accept) {
            accept = gsl_rng_uniform(rng) <= exp(-DeltaE/(k_B*T));
        }

        if (accept) {
            E_curr = E_trial;
        } else {
            // undo the swap
            temp_atype = atype[swap_idx_a];
            atype[swap_idx_a] = atype[swap_idx_b];
            atype[swap_idx_b] = temp_atype;
        }

        // calculate various instantaneous quantities
        // save the quantities somehow
        E_out[i_step] = E_curr;
        P_out[i_step] = calc_P(N_atoms, atype, pos);
        r_out[i_step] = calc_r(N_atoms, atype, nn_idxs);

        print_progress(i_step,0,n_steps-1,false);
    }
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

    int N_atoms = 2*n_cells*n_cells*n_cells;

    int* atype = (int*)malloc(N_atoms*sizeof(int));
    int** pos = create_2D_int_array(N_atoms, 4);
    int** nn_idxs = create_2D_int_array(N_atoms, 8);
    int idx_by_pos[n_cells][n_cells][n_cells][2];

    construct_bcc_binary_alloy(n_cells, atype, pos, nn_idxs, idx_by_pos);

    // print the constructed arrays for debugging:
    // printf("atype = \n");
    // print_vector_int(atype, N_atoms);
    // printf("pos = \n");
    // print_matrix_int(pos, N_atoms, 4);
    // printf("nn_idxs = \n");
    // print_matrix_int(nn_idxs, N_atoms, 8);

    // Monte Carlo Simulation, Metropolis Algorithm
    double T = 600; // K
    int n_steps = 300000;
    double* E = (double*)malloc(n_steps*sizeof(double));
    double* P = (double*)malloc(n_steps*sizeof(double));
    double* r = (double*)malloc(n_steps*sizeof(double));

    metropolis_algorithm(N_atoms, n_steps, T, atype, pos, nn_idxs, rng, E, P, r);

    // write it to a file
    FILE* file = fopen("data/H2a_sim.csv", "w");
    fprintf(file, "# {\"T[K]\": %.2f, \"n_steps\": %i}\n", T, n_steps);
    fprintf(file, "# i_step, E[eV], P, r\n");
    for (int i_step=0; i_step<n_steps; i_step++) {
        fprintf(file, "%i, %.10e, %.10e, %.10e\n", i_step, E[i_step], P[i_step], r[i_step]);
    }
    fclose(file);    

    // tidy up
    free(atype);
    destroy_2D_int_array(pos);
    destroy_2D_int_array(nn_idxs);
    gsl_rng_free(rng);
    return 0;
}
