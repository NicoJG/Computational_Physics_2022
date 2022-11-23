#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools_.h"

// Constants
const double pi = 3.14159265359;

/*
 * Calculate the acceleration (for one timestep)
 * @a - vector that is filled with the accelerations
 * @u - vector with the current positions
 * @m - mass of the particles (all have the same mass)
 * @kappa - Spring constant
 * @alpha - anharmonic coupling strength
 * @n_particles - number of particles and length of a and u
 */
void calc_acc(double *a, double *u, double m, double kappa, double alpha, int n_particles)
{
    double u_left_wall = 0;
    double u_right_wall = 0;
    // Calculate the acceleration for the boundaries because u_0 and u_(N+1) are the walls
    a[0] = kappa * (u[1] - 2*u[0] + u_left_wall) * (1 + alpha * (u[1] - u_left_wall)) / m;
    a[n_particles-1] = kappa * (u_right_wall - 2*u[n_particles-1] + u[n_particles-2]) * (1 + alpha * (u_right_wall - u[n_particles-2])) / m;
    
    // Calculate the accelerations
    for (int i = 1; i < n_particles - 1; i++){
        a[i] = kappa * (u[i+1] - 2*u[i] + u[i-1]) * (1 + alpha * (u[i+1] - u[i-1])) / m;
    }
}

/* generate the transformation matrix to and from normal coordinates
 * @trans_matrix - shape: n_particles x n_particles
 * @n_particles - number of particles / normal modes
 */
void transformation_matrix(double** trans_matrix, int n_particles){
    double factor = 1 / ((double)n_particles + 1);
    for(int k = 0; k < n_particles; k++){
        for(int i = 0; i < n_particles; i++){
            trans_matrix[k][i] = sqrt(2 * factor) * sin((i + 1) * (k + 1) * pi * factor);
        }
    }
}

/*
 * Calculate the energy in each mode for one timestep
 * @E_out -  array to be filled with the modal energies: sizeof(E_out)=n_particles
 * @Q - array of the normal positions: sizeof(Q)=n_particles
 * @P - array of the normal momenta: sizeof(P)=n_particles
 * @n_particles - number of particles(normal modes)
 * @m - mass of all the particles
 * @kappa - spring constant
 */
void calc_modal_energies(double* E_out, 
                        double* Q, double* P,
                        int n_particles,
                        double m, double kappa){
    for(int k=0;k<n_particles;k++){
        double omega_k = 2*sqrt(kappa/m)*sin((k+1)*pi/(2*(n_particles+1)));
        E_out[k] = 0.5*(P[k]*P[k] + omega_k*omega_k*Q[k]*Q[k]);
    }
}

/*
 * Perform the velocity verlet alogrithm and calculate the modal energies
 * @t_out - array to be filled with the timesteps that are saved sizeof(t_out)=(n_timesteps+1)
 * @E_out - 2D array to be filled with the modal energies : sizeof(E_out)=(n_timesteps+1) x n_particles
 * @u0 : array with the initial positions :sizeof(u0)=n_particles
 * @p0 : array with the initial momenta :sizeof(p0)=n_particles
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @dt - timestep
 * @n_save_skips - how often the energy calculation and saving should be skipped in a row
 * @m - mass of all the particles
 * @kappa - Spring constant
 * @alpha - anharmonic coupling constant
 */
void velocity_verlet(double* t_out, double** E_out, double* u0, double* p0,
                    int n_timesteps, int n_particles, 
                    double dt, int n_save_skips, double m, 
                    double kappa, double alpha)
{
    int skip_counter = 0; // counts how often the saving was skipped
    int save_step = 0;
    if(m != 1.){
        perror("m is not 1! this code does not account for this case!");
        exit(1);
    } // we just say q=u and v=p because m=1

    // initialize arrays for the current timestep
    double q[n_particles];
    double v[n_particles];
    double a[n_particles];
    double Q[n_particles];
    double P[n_particles];
    double E[n_particles];
    double t = 0;

    // apply initial conditions
    for(int i=0;i<n_particles;i++){
        q[i] = u0[i];
        v[i] = p0[i];
    }

    // generate the transformation matrix
    double** trans_matrix = create_2D_array(n_particles,n_particles);
    transformation_matrix(trans_matrix, n_particles);

    // Transform into the normal coordinates
    matrix_vector_product(Q, trans_matrix, q, n_particles, n_particles);
    matrix_vector_product(P, trans_matrix, v, n_particles, n_particles);

    // Calculate the modal energies
    calc_modal_energies(E, Q, P, n_particles, m, kappa);
    
    /* Save the Energies of the modes*/
    for(int k=0;k<n_particles;k++){
        E_out[0][k] = E[k];
    }

    // Save the time
    t_out[0] = t;

    // a(t=0)
    calc_acc(a, q, m, kappa, alpha, n_particles);

    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];
        }
        
        /* a(t+dt) */
        calc_acc(a, q, m, kappa, alpha, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }

        //increase the time
        t += dt;

        // check if this timestep should be skipped (E2_4)
        // (save the last timestep regardless of the skip counter)
        if(skip_counter<n_save_skips && i<n_timesteps) {
            skip_counter++;
            continue;
        }

        skip_counter = 0;
        save_step++;

        printf("\33[2K\rCurrent savestep: %i at timestep: %i ; t = %.2f", save_step+1, i+1, t);

        // Transform into the normal coordinates
        matrix_vector_product(Q, trans_matrix, q, n_particles, n_particles);
        matrix_vector_product(P, trans_matrix, v, n_particles, n_particles);

        // Calculate the modal energies
        calc_modal_energies(E, Q, P, n_particles, m, kappa);
        
        /* Save the Energies of the modes*/
        for(int k=0;k<n_particles;k++){
            E_out[save_step][k] = E[k];
        }

        // save the time
        t_out[save_step] = t;
    }
    printf("\n");

    destroy_2D_array(trans_matrix);
}

void simulate_and_save(char* filename, int n_particles, 
                    double t_max, double dt, int n_save_skips,
                    double m, double kappa, double alpha){
    int n_timesteps = t_max/dt;
    int n_save_steps = (n_timesteps / (n_save_skips+1)) +1;
    if(n_save_skips>0){
        if((n_timesteps % (n_save_skips+1)) > 0){
            n_save_steps++;
        }
    }
    printf("n_timesteps = %i ; n_save_steps = %i\n",n_timesteps, n_save_steps);
    // prepare output arrays (they are initialized to 0)
    double* t = (double*)calloc(n_save_steps,sizeof(double));
    double** E = create_2D_array(n_save_steps,n_particles);

    // arrays for the initial conditions:
    double u0[n_particles];
    double p0[n_particles];
    double Q0[n_particles];
    double P0[n_particles];

    // initial conditions so that only mode k=0 is active
    t[0] = 0;
    E[0][0] = n_particles;
    P0[0] = sqrt(2*E[0][0]);
    // transform it to the atom coordinates (is the same as the other way)
    double** trans_matrix = create_2D_array(n_particles,n_particles);
    transformation_matrix(trans_matrix, n_particles);
    matrix_vector_product(u0,trans_matrix,Q0, n_particles, n_particles);
    matrix_vector_product(p0,trans_matrix,P0, n_particles, n_particles);

    // perform the simulation
    velocity_verlet(t,E,u0,p0,n_timesteps,n_particles,dt,n_save_skips,m,kappa,alpha);

    // ******************************
    // write to the datafile
    FILE* file = fopen(filename, "w");
    // header
    fprintf(file, "# dt = %.5f ps ; t_max = %.5e ; n_timesteps = %i ; n_save_steps = %i ; \n", dt, t_max, n_timesteps, n_save_steps);
    fprintf(file, "# n_particles = %i ; mass = %.2e ; kappa = %.4e ; alpha = %.4e\n", n_particles, m, kappa, alpha);
    // table header
    fprintf(file, "# t");
    for(int j=0;j<n_particles;j++){
        fprintf(file,", E_%i", j);
    }
    fprintf(file, "\n");
    // write the data
    for(int i=0;i<(n_save_steps);i++){
        fprintf(file, "%.5e", t[i]);
        for(int k=0;k<n_particles;k++){
            fprintf(file, ", %.5e", E[i][k]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    // **********************************

    // free stuff
    free(t);
    destroy_2D_array(E);
    destroy_2D_array(trans_matrix);
}

int main(){
    double m = 1;
    double kappa = 1;
    double alpha = 0;

    int n_particles = 32;
    double t_max = 25000;
    double dt = 0.1;

    int n_save_skips = 0;

    printf("Simulate E2_2 ...\n");
    simulate_and_save("data/E2_2.csv", n_particles, t_max, dt, n_save_skips, m, kappa, alpha);

    alpha = 0.01;
    printf("Simulate E2_3_alpha0.01 ...\n");
    simulate_and_save("data/E2_3_alpha0.01.csv", n_particles, t_max, dt, n_save_skips, m, kappa, alpha);

    alpha = 0.1;
    printf("Simulate E2_3_alpha0.1 ...\n");
    simulate_and_save("data/E2_3_alpha0.1.csv", n_particles, t_max, dt, n_save_skips, m, kappa, alpha);
    
    t_max = 1e6;
    n_save_skips = 999;
    alpha = 0.01;
    printf("Simulate E2_4_alpha0.01 ...\n");
    simulate_and_save("data/E2_4_alpha0.01.csv", n_particles, t_max, dt, n_save_skips, m, kappa, alpha);

    alpha = 0.1;
    printf("Simulate E2_4_alpha0.1 ...\n");
    simulate_and_save("data/E2_4_alpha0.1.csv", n_particles, t_max, dt, n_save_skips, m, kappa, alpha);
}