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
const double E_AA = -436e-3; // eV E_CuCu
const double E_BB = -113e-3; // eV E_ZnZn
const double E_AB = -294e-3; // eV E_CuZn
const int n_cells = 10;
const int N_atoms = 2*n_cells*n_cells*n_cells;

void count_bonds(int N_atoms, int* atype, int** nn_idxs,
                int* N_AA, int* N_BB, int* N_AB) {
    // count the types of neares neighbor pairs
    (*N_AA) = 0;
    (*N_BB) = 0;
    (*N_AB) = 0;

    for (int i=0; i<N_atoms; i++) {
        int atype_a = atype[i];
        for (int j=0; j<8; j++) {
            int atype_b = atype[nn_idxs[i][j]];
            if (atype_a == 0 && atype_b == 0) {
                (*N_AA)++;
            } else if (atype_a == 1 && atype_b == 1) {
                (*N_BB)++;
            } else {
                (*N_AB)++;
            }
        }
    }
    // but now we counted each bond twice
    if ((*N_AA)%2 == 1 || (*N_BB)%2 == 1 || (*N_AB)%2 == 1) {
        perror("ERROR: number of pairs should be even");
        exit(1);
    }
    (*N_AA) /= 2;
    (*N_BB) /= 2;
    (*N_AB) /= 2;
}

void update_bond_counts_after_swap(int N_atoms, int* atype, 
                int** nn_idxs, int** pos,
                int swap_idx_a, int swap_idx_b,
                int* N_AA, int* N_BB, int* N_AB, int* N_Aa) {
    int atype_a = atype[swap_idx_a]; // old atype at position a
    int atype_b = atype[swap_idx_b]; // old atype at position b
    if (atype_a == atype_b) {
        // the swap introduces no change
        return;
    }

    // update the Number of A atoms in sublattice a
    // if position a is on a different sublattice than position b, N_Aa can change
    if (pos[swap_idx_a][3] != pos[swap_idx_a][3]) {
        if (pos[swap_idx_a][3] == 0) {
            // position a is on sublattice a
            if (atype_a == 0) {
                (*N_Aa)--;
            } else {
                (*N_Aa)++;
            }
        } else {
            // position b is on sublattice a
            if (atype_b == 0) {
                (*N_Aa)--;
            } else {
                (*N_Aa)++;
            }
        }
    }

    // update bonds at position a
    for (int j=0; j<8; j++) {
        int atype_nn = atype[nn_idxs[swap_idx_a][j]];
        // substract old bond and add new bond
        if (atype_a == 0 && atype_nn == 0) {
            (*N_AA)--;
            (*N_AB)++;
        } else if (atype_a == 1 && atype_nn == 1) {
            (*N_BB)--;
            (*N_AB)++;
        } else if (atype_a == 0 && atype_nn == 1) {
            (*N_AB)--;
            (*N_BB)++;
        } else if (atype_a == 1 && atype_nn == 0) {
            (*N_AB)--;
            (*N_AA)++;
        }
    }

    // update bonds at position b
    for (int j=0; j<8; j++) {
        int atype_nn = atype[nn_idxs[swap_idx_b][j]];

        // if the nn of b is a then the change is already acounted for
        if (nn_idxs[swap_idx_b][j] == swap_idx_a) {
            continue;
        }

        // substract old bond and add new bond
        if (atype_b == 0 && atype_nn == 0) {
            (*N_AA)--;
            (*N_AB)++;
        } else if (atype_b == 1 && atype_nn == 1) {
            (*N_BB)--;
            (*N_AB)++;
        } else if (atype_b == 0 && atype_nn == 1) {
            (*N_AB)--;
            (*N_BB)++;
        } else if (atype_b == 1 && atype_nn == 0) {
            (*N_AB)--;
            (*N_AA)++;
        }
    }
}

double calc_energy(int N_AA, int N_BB, int N_AB) {
    return N_AA*E_AA + N_BB*E_BB + N_AB*E_AB;
}

double calc_P(int N_atoms, int N_Aa) {
    // long-range order parameter P
    int N_A = N_atoms/2;
    return 2.*((double)N_Aa)/N_A - 1.;
}

double calc_r(int N_atoms, int N_AB) {
    // short-range order parameter r
    int N = N_atoms/2;
    return (N_AB-4.*N)/(4.*N);
}

void metropolis_algorithm(int N_atoms, int n_steps, double T, 
                        int* atype, int** pos, int** nn_idxs, gsl_rng* rng, 
                        double* E_out, double* P_out, double* r_out) {
    int N_AA, N_BB, N_AB, N_Aa;
    count_bonds(N_atoms,atype,nn_idxs,&N_AA, &N_BB, &N_AB);
    double E_curr = calc_energy(N_AA, N_BB, N_AB);
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

        // save the current bond counts if the trial change is not accepted
        int N_AA_trial = N_AA;
        int N_BB_trial = N_BB;
        int N_AB_trial = N_AB;
        int N_Aa_trial = N_Aa;

        // check if the swap should be accepted as the new configuration
        update_bond_counts_after_swap(N_atoms,atype,nn_idxs, pos,
                            swap_idx_a,swap_idx_b,
                            &N_AA_trial, &N_BB_trial, &N_AB_trial, &N_Aa_trial);
        double E_trial = calc_energy(N_AA_trial, N_BB_trial, N_AB_trial);
        double DeltaE = E_trial - E_curr;
        bool accept = (DeltaE <= 0);
        if (!accept) {
            accept = gsl_rng_uniform(rng) <= exp(-DeltaE/(k_B*T));
        }

        if (accept) {
            E_curr = E_trial;
            // perform the swap
            int temp_atype = atype[swap_idx_a];
            atype[swap_idx_a] = atype[swap_idx_b];
            atype[swap_idx_b] = temp_atype;
            // update the counts
            N_AA = N_AA_trial;
            N_BB = N_BB_trial;
            N_AB = N_AB_trial;
            N_Aa = N_Aa_trial;
        }
        
        if (8*N_atoms != N_AA+N_BB+N_AB) {
            printf("ERROR: bonds calculated wrong: N_AA = %i, N_BB = %i")
        }
        
        // calculate various instantaneous quantities
        // save the quantities somehow
        E_out[i_step] = E_curr;
        P_out[i_step] = calc_P(N_atoms, N_Aa);
        r_out[i_step] = calc_r(N_atoms, N_AB);

        print_progress(i_step,0,n_steps-1,false);
    }
}

void perform_simulation(double T, int n_eq_steps, int n_steps, 
                    bool save_steps, char* save_steps_file_path, gsl_rng* rng,
                    double* E_avg, double* P_avg, double* r_avg, double* C) {
    int* atype = (int*)malloc(N_atoms*sizeof(int));
    int** pos = create_2D_int_array(N_atoms, 4);
    int** nn_idxs = create_2D_int_array(N_atoms, 8);
    int idx_by_pos[n_cells][n_cells][n_cells][2];

    construct_bcc_binary_alloy(n_cells, atype, pos, nn_idxs, idx_by_pos);

    // Monte Carlo Simulation, Metropolis Algorithm
    int n_tot_steps = n_eq_steps+n_steps;
    double* E = (double*)malloc(n_tot_steps*sizeof(double));
    double* P = (double*)malloc(n_tot_steps*sizeof(double));
    double* r = (double*)malloc(n_tot_steps*sizeof(double));

    metropolis_algorithm(N_atoms, n_tot_steps, T, atype, pos, nn_idxs, rng, E, P, r);

    // calculate the average quantities
    *E_avg = average(E+n_eq_steps, n_steps);
    *P_avg = average(P+n_eq_steps, n_steps);
    *r_avg = average(r+n_eq_steps, n_steps);

    // calculate the heat capacity: C = 1/k_BT * (<E^2>-<E>^2)
    double* E_squared = (double*)malloc(n_steps*sizeof(double));
    elementwise_multiplication(E_squared, E+n_eq_steps, E+n_eq_steps, n_steps);
    double E_squared_avg = average(E_squared, n_steps);
    *C = (E_squared_avg - (*E_avg)*(*E_avg))/(k_B*T);

    if (save_steps) {
        printf("Save %i steps for T = %.4f K ...\n", n_tot_steps, T);
        // write it to a file
        FILE* file = fopen(save_steps_file_path, "w");
        fprintf(file, "# {\"T[K]\": %.2f, \"n_eq_steps\":%i, \"n_steps\": %i,\n", T, n_eq_steps, n_steps);
        fprintf(file, "# \"P\": %.10f, \"E[eV]\": %.10f, \"C[eV/K]\": %.10f, \"r\": %.10f}\n", *P_avg, *E_avg, *C, *r_avg);
        fprintf(file, "# i_step, E[eV], P, r\n");
        for (int i_step=0; i_step<n_steps; i_step++) {
            fprintf(file, "%i, %.10e, %.10e, %.10e\n", i_step, E[i_step], P[i_step], r[i_step]);
        }
        fclose(file);    
    }

    // tidy up
    free(E);
    free(E_squared);
    free(P);
    free(r);
    free(atype);
    destroy_2D_int_array(pos);
    destroy_2D_int_array(nn_idxs);
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

    int n_T_saves = 3;
    int n_T = 20;
    int T_min = 100;
    int T_max = 1200;
    int dT = (T_max-T_min)/(n_T-1);
    int dT_save = (T_max-T_min)/(n_T_saves-1);

    double* T = (double*)malloc(n_T*sizeof(double));
    for (int i=0; i<n_T; i++) {
        T[i] = T_min + i*dT;
    }


    double* E = (double*)malloc(n_T*sizeof(double));
    double* P = (double*)malloc(n_T*sizeof(double));
    double* r = (double*)malloc(n_T*sizeof(double));
    double* C = (double*)malloc(n_T*sizeof(double));

    int n_eq_steps_save = 1e6;
    int n_steps_save = 1e6;

    int n_eq_steps = 1e7;
    int n_steps = 1e7;

    // save a few full simulations
    for (int i=0; i<n_T_saves; i++) {
        // save 3 full simulations
        double T_ = T_min + i*dT_save;
        char save_step_file_path[100];
        sprintf(save_step_file_path, "data/full_simulations/H2a_simsteps_T%.2fK.csv", T_);

        perform_simulation(T_, n_eq_steps_save, n_steps_save, 
                        true, save_step_file_path, rng,
                        E, P, r, C);
    }

    // perform alot of simulations without saving
    for (int i=0; i<n_T; i++) {
        perform_simulation(T[i], n_eq_steps, n_steps, 
                        false, "", rng,
                        E+i, P+i, r+i, C+i);
        printf("%i/%i done. (T = %.2f K)\n", i+1,n_T,T[i]);
    }

    // save the results
    FILE* file = fopen("data/H2a.csv", "w");
    fprintf(file, "# {\"n_eq_steps\": %i, \"n_steps\": %i}\n", n_eq_steps, n_steps);
    fprintf(file, "# T, E, P, r, C\n");
    for (int i=0; i<n_T; i++) {
        fprintf(file, "%.5f, %.10f, %.10f, %.10f, %.10f\n", T[i], E[i], P[i], r[i], C[i]);
    }
    fclose(file);


    // tidy up
    free(T);
    free(E);
    free(P);
    free(r);
    free(C);
    gsl_rng_free(rng);
    return 0;
}
