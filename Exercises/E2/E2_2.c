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


/*
 * Perform the velocity verlet alogrithm
 * @u_out - 2D array to be filled with the positions (holds the initial conditions) : sizeof(u_out)=(n_timesteps+1) x n_particles
 * @p_out - 2D array to be filled with the momenta (holds the initial conditions) : sizeof(p_out)=(n_timesteps+1) x n_particles
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @dt - timestep
 * @m - mass of all the particles
 * @kappa - Spring constant
 * @alpha - anharmonic coupling constant
 */
void velocity_verlet(double** u_out, double** p_out,
                    int n_timesteps, int n_particles, 
                    double dt, double m, 
                    double kappa, double alpha)
{
    // initialize arrays for the current timestep
    double q[n_particles];
    double v[n_particles];
    double a[n_particles];

    // apply initial conditions
    for(int i=0;i<n_particles;i++){
        q[i] = u_out[0][i];
        v[i] = p_out[0][i]/m;
    }

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
        
        /* Save the displacement and momentum of the atoms */
        for(int j=0;j<n_particles;j++){
            u_out[i][j] = q[j];
            p_out[i][j] = m*v[j];
        }
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
 * Transform the coordinates u and p to the normal coordinates Q and P
 * @q_in - 2D array of the positions: sizeof(u)=(n_timesteps+1) x n_particles
 * @p_in - 2D array of the momenta: sizeof(p)=(n_timesteps+1) x n_particles
 * @Q_out - 2D array to be filled with the normal positions: sizeof(Q)=(n_timesteps+1) x n_particles
 * @P_out - 2D array to be filled with the normal momenta: sizeof(P)=(n_timesteps+1) x n_particles
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @m - mass of all the particles
 */
void transform_coords(double** q_in, double** p_in,
                      double** Q_out, double** P_out,
                      int n_timesteps, int n_particles,
                      double m){
    if(m != 1.){
        perror("m is not 1! this code does not account for this case!");
        exit(1);
    } 
    // create the transformation matrix
    double** trans_matrix = create_2D_array(n_particles,n_particles);
    transformation_matrix(trans_matrix, n_particles);

    // apply the transformations u->Q and p->P (matrix multiplication)
    for(int t_i=0; t_i < n_timesteps+1; t_i++){
        matrix_vector_product(Q_out[t_i], trans_matrix, q_in[t_i], n_particles, n_particles);
        matrix_vector_product(P_out[t_i], trans_matrix, p_in[t_i], n_particles, n_particles);
    }
    destroy_2D_array(trans_matrix);
}

/*
 * Calculate the energy in each mode in all timesteps
 * @E_out - 2D array to  be filled with the modal energies: sizeof(E_out)=(n_timesteps+1) x n_particles
 * @Q - 2D array of the normal positions: sizeof(Q)=(n_timesteps+1) x n_particles
 * @P - 2D array of the normal momenta: sizeof(P)=(n_timesteps+1) x n_particles
 * @n_particles - number of particles(normal modes)
 * @n_timesteps - number of timesteps
 * @m - mass of all the particles
 * @kappa - spring constant
 */
void calc_modal_energies(double** E_out, 
                        double** Q, double** P,
                        int n_particles, int n_timesteps,
                        double m, double kappa){
    for(int t_i=0;t_i<n_timesteps+1;t_i++){
        for(int k=0;k<n_particles;k++){
            double omega_k = 2*sqrt(kappa/m)*sin((k+1)*pi/(2*(n_particles+1)));
            printf("k=%i; omega=%.5e\n",k,omega_k);
            E_out[t_i][k] = 0.5*(P[t_i][k]*P[t_i][k] + omega_k*omega_k*Q[t_i][k]*Q[t_i][k]);
        }
    }
}


int main(){
    // simulation constants (dimensionless units)
    int n_particles = 32;
    double t_max = 250;
    double dt = 0.1;
    int n_timesteps = t_max/dt;
    printf("n_timesteps = %i\n",n_timesteps);

    double m = 1;
    double kappa = 1;
    double alpha = 0;

    // prepare output arrays (they are initialized to 0)
    double** u = create_2D_array(n_timesteps+1,n_particles);
    double** p = create_2D_array(n_timesteps+1,n_particles);
    double** P = create_2D_array(n_timesteps+1,n_particles);
    double** Q = create_2D_array(n_timesteps+1,n_particles);
    double** E = create_2D_array(n_timesteps+1,n_particles);

    // initial conditions so that only mode k=0 is active
    double E_0 = n_particles;
    P[0][0] = sqrt(2*E_0);
    transform_coords(Q,P,u,p,0,n_particles, m); // the transformation is the same as the inverse if m=1

    // perform the simulation
    velocity_verlet(u,p,n_timesteps,n_particles,dt,m,kappa,alpha);

    // transform into normal coordinates
    transform_coords(u,p,Q,P,n_timesteps,n_particles,m);

    // calculate the energies
    calc_modal_energies(E,Q,P,n_particles,n_timesteps,m,kappa);

    // ******************************
    // write to the datafile
    FILE* file = fopen("E2_2.csv", "w");
    // header
    fprintf(file, "# dt = %.5f ps ; n_timesteps = %i ; n_particles = %i\n", dt, n_timesteps, n_particles);
    fprintf(file, "# mass = %.2e ; kappa = %.4e ; alpha = %.4e\n", m, kappa, alpha);
    // table header
    fprintf(file, "# t");
    char* vars = "upQPE";
    int n_vars = 5;
    for(int i=0;i<n_vars;i++){
        for(int j=0;j<n_particles;j++){
            fprintf(file,", %c_%i", vars[i], j);
        }
    }
    fprintf(file, "\n");
    // write the data
    double** arrs[n_vars];
    arrs[0] = u;
    arrs[1] = p;
    arrs[2] = Q;
    arrs[3] = P;
    arrs[4] = E;
    for(int i=0;i<(n_timesteps+1);i++){
        fprintf(file, "%.5e", i*dt);
        for(int j=0;j<n_vars;j++){
            for(int k=0;k<n_particles;k++){
                fprintf(file, ", %.5e", arrs[j][i][k]);
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
    // **********************************

    // free stuff
    destroy_2D_array(u);
    destroy_2D_array(p);
    destroy_2D_array(Q);
    destroy_2D_array(P);
    destroy_2D_array(E);
}