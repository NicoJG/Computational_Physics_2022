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
const int n_cells = 10;
const int N_atoms = 2*n_cells*n_cells*n_cells;

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
    if (pos[swap_idx_a][3] != pos[swap_idx_b][3]) {
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

        // if the nn of a is b then the change is already acounted for
        if (nn_idxs[swap_idx_a][j] == swap_idx_b) {
            continue;
        }

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

void metropolis_algorithm(int N_atoms, int n_steps, double T, 
                        int* atype, int** pos, int** nn_idxs, gsl_rng* rng, 
                        double* E_out, double* P_out, double* r_out) {
    int N_AA, N_BB, N_AB, N_Aa=0;
    count_bonds(N_atoms,atype,nn_idxs,&N_AA, &N_BB, &N_AB);
    // count A atoms on a sublattice
    for (int i=0; i<N_atoms; i++) {
        if (atype[i] == 0 && pos[i][3] == 0) {
            N_Aa++;
        }
    }

    double E_curr = calc_energy(N_AA, N_BB, N_AB);
    print_progress(0,0,n_steps-1,true);
    // calc the energy of the current configuration
    for (int i_step=0; i_step<n_steps; i_step++) {
        // make trial change by swapping two atoms
        // choose which atoms to swap
        int swap_idx_a = gsl_rng_uniform_int(rng, N_atoms);
        int swap_idx_b = gsl_rng_uniform_int(rng, N_atoms);
        while (swap_idx_a == swap_idx_b || atype[swap_idx_a] == atype[swap_idx_b]) {
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
        //DeltaE = E_curr - E_trial;
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

        //if (T>600 && i_step>8e6) {
        //    printf("N_AA = %i, N_BB = %i, N_AB = %i, N_Aa = %i, E = %f\n",N_AA,N_BB,N_AB,N_Aa,E_curr);
        //}
        
        if (4*N_atoms != N_AA+N_BB+N_AB || N_AA<0 || N_BB<0 || N_AB<0) {
            printf("\nERROR: bonds calculated wrong: N_AA = %i, N_BB = %i, N_AB = %i\n", N_AA, N_BB, N_AB);
            exit(1);
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
                    int* atype, int** pos, int** nn_idxs, int idx_by_pos[n_cells][n_cells][n_cells][2],
                    bool save_steps, char* save_steps_file_path, gsl_rng* rng, int n_skip_saves,
                    double* E_avg, double* P_avg, double* r_avg,
                    double* E_std, double* P_std, double* r_std, double* C) {

    // Monte Carlo Simulation, Metropolis Algorithm
    int n_tot_steps = n_eq_steps+n_steps;
    double* E = (double*)malloc(n_tot_steps*sizeof(double));
    double* P = (double*)malloc(n_tot_steps*sizeof(double));
    double* r = (double*)malloc(n_tot_steps*sizeof(double));

    metropolis_algorithm(N_atoms, n_tot_steps, T, atype, pos, nn_idxs, rng, E, P, r);

    // calculate the average quantities (without the equilibration steps)
    *E_avg = average(E+n_eq_steps, n_steps);
    *P_avg = average(P+n_eq_steps, n_steps);
    *r_avg = average(r+n_eq_steps, n_steps);
    *E_std = standard_deviation(E+n_eq_steps, n_steps);
    *P_std = standard_deviation(P+n_eq_steps, n_steps);
    *r_std = standard_deviation(r+n_eq_steps, n_steps);

    // calculate the heat capacity: C = 1/k_BT * (<E^2>-<E>^2)
    double* E_squared = (double*)malloc(n_steps*sizeof(double));
    elementwise_multiplication(E_squared, E+n_eq_steps, E+n_eq_steps, n_steps);
    double E_squared_avg = average(E_squared, n_steps);
    *C = (E_squared_avg - (*E_avg)*(*E_avg))/(k_B*T*T);

    if (save_steps) {
        int n_save_steps = n_tot_steps / n_skip_saves;
        printf("Save %i steps from %i for T = %.4f K ...\n", n_save_steps, n_tot_steps, T);
        // write it to a file
        FILE* file = fopen(save_steps_file_path, "w");
        fprintf(file, "# {\"T[K]\": %.2f, \"n_eq_steps\":%i, \"n_steps\": %i, \"n_skip_saves\": %i,\n", T, n_eq_steps, n_steps, n_skip_saves);
        fprintf(file, "# \"P\": %.10f, \"E[eV]\": %.10f, \"C[eV/K]\": %.10f, \"r\": %.10f}\n", *P_avg, *E_avg, *C, *r_avg);
        fprintf(file, "# i_step, E[eV], P, r\n");
        for (int i_step=0; i_step<n_eq_steps+n_steps; i_step+=n_skip_saves) {
            fprintf(file, "%i, %.10e, %.10e, %.10e\n", i_step, E[i_step], P[i_step], r[i_step]);
        }
        fclose(file);    
    }

    // tidy up
    free(E);
    free(E_squared);
    free(P);
    free(r);
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

    // make 3 sections of temperature, where around the critical temperature it more dense
    int T_0 = 300;
    int T_1 = 600;
    int T_2 = 800;
    int T_3 = 1000;
    // delta T steps in the sections
    int dT01 = 10;
    int dT12 = 2;
    int dT23 = 10;
    // number of T steps in each section
    int n_T01 = (T_1-T_0)/dT01;
    int n_T12 = (T_2-T_1)/dT12;
    int n_T23 = (T_3-T_2)/dT23;
    int n_T = n_T01 + n_T12 + n_T23;

    int n_T_saves = 10; // number of T's to fully save

    double* T = (double*)malloc(n_T*sizeof(double));
    double temp_T = T_0;
    for (int i=0; i<n_T; i++) {
        T[i] = temp_T;
        if (T_0 <= temp_T && temp_T < T_1) {
            temp_T += dT01;
        } else if (T_1 <= temp_T && temp_T < T_2) {
            temp_T += dT12;
        } else if (T_2 <= temp_T && temp_T < T_3) {
            temp_T += dT23;
        }
    }

    double* E = (double*)malloc(n_T*sizeof(double));
    double* P = (double*)malloc(n_T*sizeof(double));
    double* r = (double*)malloc(n_T*sizeof(double));
    double* E_std = (double*)malloc(n_T*sizeof(double));
    double* P_std = (double*)malloc(n_T*sizeof(double));
    double* r_std = (double*)malloc(n_T*sizeof(double));
    double* C = (double*)malloc(n_T*sizeof(double));

    int n_eq_steps = 0.5e6;
    int n_steps = 4.5e6;
    int n_skip_saves = 10;

    // init the lattice completely ordered
    int* atype = (int*)malloc(N_atoms*sizeof(int));
    int** pos = create_2D_int_array(N_atoms, 4);
    int** nn_idxs = create_2D_int_array(N_atoms, 8);
    int idx_by_pos[n_cells][n_cells][n_cells][2];
    construct_bcc_binary_alloy(n_cells, atype, pos, nn_idxs, idx_by_pos);

    // perform the simulations
    for (int i=0; i<n_T; i++) {
        bool save = false;
        char save_step_file_path[100];
        // check if this step should be saved
        if (i%((n_T-1)/(n_T_saves-1)) == 0) {
            save = true;
            sprintf(save_step_file_path, "data/full_simulations/H2a_simsteps_%0*i_T%.0fK.csv", ((int)log10(n_T))+1, i, T[i]);
        }
        perform_simulation(T[i], n_eq_steps, n_steps, 
                        atype, pos, nn_idxs, idx_by_pos,
                        save, save_step_file_path, rng, n_skip_saves,
                        E+i, P+i, r+i,E_std+i, P_std+i, r_std+i, C+i);
        printf("%i/%i done. (T = %.2f K)\n", i+1,n_T,T[i]);
    }

    // save the results
    FILE* file = fopen("data/H2a.csv", "w");
    fprintf(file, "# {\"n_eq_steps\": %i, \"n_steps\": %i}\n", n_eq_steps, n_steps);
    fprintf(file, "# T[K], E[eV], E_std[eV], P, P_std, r, r_std, C[eV/K]\n");
    for (int i=0; i<n_T; i++) {
        fprintf(file, "%.5f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n", 
                T[i], E[i], E_std[i], P[i], P_std[i], r[i], r_std[i], C[i]);
    }
    fclose(file);

    // tidy up
    free(T);
    free(E);
    free(P);
    free(r);
    free(E_std);
    free(P_std);
    free(r_std);
    free(C);
    gsl_rng_free(rng);
    free(atype);
    destroy_2D_int_array(pos);
    destroy_2D_int_array(nn_idxs);
    return 0;
}
