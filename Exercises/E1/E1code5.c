#include <stdio.h>
#include "fft.h"

/******************************************************************************
 * E1code4
 ******************************************************************************
 * Routine that runs the velocity verlet algorithm
 * Use as template to construct your program!
 */

/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @size_of_u - the size of the position, acceleration and mass array
 */
void calc_acc(double *a, double *u, double *m, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(- 2*u[0] + u[1])/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1];
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
    }
}

/*
 * Calculate the kinetic energy of the system
 * @v - vector with the current velocities
 * @m - vector with masses of atoms sizeof(n_particles)
 * @n_particles - the size of the position and mass array
 */
double calc_E_kin(double *v, double *m, int n_particles)
{
    // implement the sum for the kinetic energy
    double E_kin = 0;
    for(int i=0;i<n_particles;i++){
        E_kin += 0.5 * m[i] * v[i] * v[i];
    }

    return E_kin;
}

/*
 * Calculate the potential energy of the system
 * @u - vector with the current positions
 * @kappa - Spring constant
 * @n_particles - the size of the position array
 */
double calc_E_pot(double *u, double kappa, int n_particles)
{
    // Extend the position array to include also the two walls
    double q[n_particles+2];
    q[0] = 0;
    q[n_particles+1] = 0;
    for(int i=0;i<n_particles;i++){
        q[i+1] = u[i];
    }

    // implement the sum for the potential energy
    double E_pot = 0;
    for(int i=0;i<(n_particles+1);i++){
        E_pot += (q[i+1] - q[i])*(q[i+1] - q[i]);
    }
    E_pot *= kappa*0.5;

    return E_pot;
}

/*
 * Perform the velocity verlet alogrithm and calculate the energies
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @q_n - position of the n'th atom : sizeof(q_n) = n_timesteps+1
 * @v_n - velocity of the n'th atom : sizeof(v_n) = n_timesteps+1
 * @E_kin - kinetic energy per timestep : sizeof(E_kin) = n_timesteps+1
 * @E_pot - potential energy per timestep : sizeof(E_kin) = n_timesteps+1
 * @dt - timestep
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant
 */
void velocity_verlet(int n_timesteps, int n_particles, 
                    double *q_1, double *q_2, double *q_3, 
                    double *v_1, double *v_2, double *v_3, 
                    double *E_kin, double *E_pot,
                    double dt, double *m, double kappa)
{
    double q[n_particles];
    double v[n_particles];
    double a[n_particles];

    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];

    v[0] = v_1[0];
    v[1] = v_2[0];
    v[2] = v_3[0];

    // Calculate the initial kinetic energy
    E_kin[0] = calc_E_kin(v, m, n_particles);

    // Calculate the initial potential energy
    E_pot[0] = calc_E_pot(q, kappa, n_particles);

    calc_acc(a, q, m, kappa, n_particles);
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
        calc_acc(a, q, m, kappa, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }

        // Calculate the current kinetic energy
        E_kin[i] = calc_E_kin(v, m, n_particles);

        // Calculate the current potential energy
        E_pot[i] = calc_E_pot(q, kappa, n_particles);
        
        /* Save the displacement of the three atoms */
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];

        /* Save the velocity of the three atoms */
        v_1[i] = v[0];
        v_2[i] = v[1];
        v_3[i] = v[2];
    }
}

/* *********************************************
 * Write all data to a file
 * @t - timesteps : size
 * @q_n - position of the n'th atom : sizeof(q_n) = n_timesteps+1
 * @v_n - velocity of the n'th atom : sizeof(v_n) = n_timesteps+1
 * @E_kin - kinetic energy per timestep : sizeof(E_kin) = n_timesteps+1
 * @E_pot - potential energy per timestep : sizeof(E_kin) = n_timesteps+1
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant

*/

/* ******************************************************
 * Main function
 */
int main(){
    double t_max = 1.; // ps
    int n_timesteps = 1000;
    double dt = t_max/n_timesteps;
    int n_particles = 3;

    // units (unit conversion)
    double m_asu = 9649; // u
    double m_asu_kg = 1.6028e-23; // kg

    double ps = 1e-12; // s

    // initial conditions
    double m[3] = {12.011/m_asu, 12.011/m_asu, 12.011/m_asu}; // m_asu with m(carbon-12)=12.011u 
    double q0[3] = {0.01, 0.005, -0.005}; // Angstroem (task)
    //double q0[3] = {0.01, 0.0, 0.0}; // Angstroem (reference plots)
    double v0[3] = {0.0, 0.0, 0.0}; // Angstroem/ps

    double kappa = 1.0; // kN/m (right units?)
    // 1 N = 1 kg*m*s^(-2)
    // => 1 kN/m = 10^3 kg*s^(-2) = 10^3 kg/m_asu * s^(-2)/ps^(-2)) m_asu ps^(-2)
    kappa *= 1e3 / m_asu_kg * ps * ps; // m_asu ps^(-2) 
    printf("1 kN/m = %.5e m_asu ps^(-2)\n",1e3 / m_asu_kg * ps * ps);

    // initialize all needed arrays
    double t[n_timesteps+1];
    double q_1[n_timesteps+1];
    double q_2[n_timesteps+1];
    double q_3[n_timesteps+1];
    double v_1[n_timesteps+1];
    double v_2[n_timesteps+1];
    double v_3[n_timesteps+1];
    double E_kin[n_timesteps+1];
    double E_pot[n_timesteps+1];

    // insert the initial conditions
    q_1[0] = q0[0];
    q_2[0] = q0[1];
    q_3[0] = q0[2];
    v_1[0] = v0[0];
    v_2[0] = v0[1];
    v_3[0] = v0[2];

    // create timestep array
    for(int i=0;i<(n_timesteps+1);i++){
        t[i] = dt*i;
    }

    // Perform the velocity verlet algorithm and calculate the energies
    velocity_verlet(n_timesteps, n_particles, 
                    q_1, q_2, q_3, 
                    v_1, v_2, v_3, 
                    E_kin, E_pot,
                    dt, m, kappa);

    // write all data to a file
    FILE* file = fopen("E1_5.csv", "w+");
    fprintf(file, "# dt = %.5f ps ; n_timesteps = %i\n", dt, n_timesteps);
    fprintf(file, "# masses = [%.2e, %.2e, %.2e] m_asu ; kappa = %.4e m_asu ps^(-2)\n", m[0], m[1], m[2], kappa);
    fprintf(file, "# t[ps], q_1[Å], q_2[Å], q_3[Å], v_1[Å/ps], v_2[Å/ps], v_3[Å/ps], E_kin[eV], E_pot[eV]\n");
    for(int i=0;i<n_timesteps+1;i++){
        fprintf(file, "%.5e, %.5e, %.5e, %.5e, %.5e, %.5e, %.5e, %.5e, %.5e \n", t[i], q_1[i], q_2[i], q_3[i], v_1[i], v_2[i], v_3[i], E_kin[i], E_pot[i]);
    }
    fclose(file);

    // calculate the powerspectra
    double frequencies[n_timesteps+1];
    fft_freq_shift(frequencies, dt, n_timesteps+1);

    double powerspec_1[n_timesteps+1];
    powerspectrum(q_1, powerspec_1, n_timesteps+1);
    powerspectrum_shift(powerspec_1, n_timesteps+1);

    double powerspec_2[n_timesteps+1];
    powerspectrum(q_2, powerspec_2, n_timesteps+1);
    powerspectrum_shift(powerspec_2, n_timesteps+1);

    double powerspec_3[n_timesteps+1];
    powerspectrum(q_3, powerspec_3, n_timesteps+1);
    powerspectrum_shift(powerspec_3, n_timesteps+1);

    // write the powerspectra to a file
    file = fopen("E1_5_powerspectrum.csv", "w+");
    fprintf(file, "# dt = %.5f ps ; n_timesteps = %i\n", dt, n_timesteps);
    fprintf(file, "# masses = [%.2e, %.2e, %.2e] m_asu ; kappa = %.4e m_asu ps^(-2)\n", m[0], m[1], m[2], kappa);
    fprintf(file, "# f[1/ps], p_1[Å^-2], p_2[Å^-2], p_3[Å^-2]\n");
    for(int i=0;i<(n_timesteps+1);i++){
        fprintf(file, "%.5e, %.5e, %.5e, %.5e \n", frequencies[i], powerspec_1[i], powerspec_2[i], powerspec_3[i]);
    }
    fclose(file);


    printf("Finished.\n");

}
