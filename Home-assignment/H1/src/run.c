#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"

// unit conversion factors
const double eV = 1.602176634e-19; // J
const double m_asu = 9649; // u
const double m_asu_kg = 1.6028e-23; // kg
const double GPa = 1e-21 / eV; // eV/Å^3
const double bar = 1e-4 / GPa; // eV/Å^3
const double celsius = 273.15; // K

// Simulation/Physics constants
const double m = 26.98153853 / m_asu; // m_asu (https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Al)
const double k_B = 8.617333262e-5; // eV/K (https://en.wikipedia.org/wiki/Boltzmann_constant#Value_in_different_units)
const double kappa_T = 0.01385 * GPa; // Å^3/eV https://www.knowledgedoor.com/2/elements_handbook/isothermal_compressibility.html
const int n_cells = 4; // number of unit cells
const int n_atoms = 4*n_cells*n_cells*n_cells;


void task1(){
    double pos[n_atoms][3];

    double a0_min = 3.98; // Å
    double a0_max = 4.08; // Å
    int n_a0 = 1000; // how many different a0 should be used?
    double da0 = (a0_max-a0_min)/(n_a0-1); // spacing of the different a0

    // prepare the output file
    FILE* file = fopen("data/H1_1.csv","w");
    fprintf(file, "# a0[Å], E_pot[eV]\n");

    for(int i=0;i<n_a0;i++){
        // calculate the cell lenghths
        double a0 = a0_min + i*da0; // Å
        double L_box = n_cells*a0;

        init_fcc(pos, n_cells, a0);
        double E_pot = get_energy_AL(pos, L_box, n_atoms);

        fprintf(file, "%.8f, %.10f\n", a0, E_pot);
    }
    fclose(file);
}

void initialize_lattice(double pos[][3], double vel[][3], double a0){
    init_matrix_stack(n_atoms,3,pos,0);
    init_matrix_stack(n_atoms,3,vel,0);

    init_fcc(pos, n_cells, a0);

    // initialize random number generator
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);

    double eps_min = -a0*0.065;
    double eps_max = a0*0.065;

    // introduce deviations
    for(int i=0;i<n_atoms;i++){
        pos[i][0] += gsl_rng_uniform(rng)*(eps_max-eps_min) + eps_min;
        pos[i][1] += gsl_rng_uniform(rng)*(eps_max-eps_min) + eps_min;
        pos[i][2] += gsl_rng_uniform(rng)*(eps_max-eps_min) + eps_min;
    }
    gsl_rng_free(rng);
}

void simulate_and_save_data(double pos[][3], double vel[][3], double* t, 
                            char* filename, 
                            int n_timesteps, double dt, double* a0,
                            double temp_scaling_time, double T_desired, 
                            double pressure_scaling_time, double P_desired){
    // if the scaling times are <= 0 then no corresponding scaling happens (no equilibration)

    double t_max = (n_timesteps+1)*dt;

    // prepare output file
    FILE* file = fopen(filename, "w");
    fprintf(file, "# {\"dt\": %.5e, \"n_timesteps\": %i, \"t_max\": %.5e, \n",dt, n_timesteps, t_max);
    fprintf(file, "# \"temp_scaling_time\": %.5e, \"T_desired\": %.5e, \n", temp_scaling_time, T_desired);
    fprintf(file, "# \"pressure_scaling_time\": %.5e, \"P_desired\": %.5e}\n", pressure_scaling_time, P_desired);
    fprintf(file, "# t[ps], E_pot[eV], E_kin[eV], T[K], P[eV/Å^3], a0[Å]\n");

    // velocity verlet:
    printf("Simulate \"%s\"...\n", filename);
    printf("dt = %.5e ; n_timesteps = %i ; t_max = %.5f\n", dt, n_timesteps, t_max);

    double a0_ = *a0;

    double L_box = n_cells*a0_;

    // prepare all needed variables for inside the algorithm
    int percent = -1; // for displaying progress
    double E_pot = 0; // eV
    double E_kin = 0; // eV
    double T = 0; // K
    double P = 0; // eV/Å^3
    double t_= *t; // ps
    double F[n_atoms][3]; //forces eV/Å
    get_forces_AL(F,pos,L_box,n_atoms);

    // evolve the system over time
    for(int i=1;i<n_timesteps+1;i++){
        // show progress
        if((i+1)*100/(n_timesteps+1) != percent) {
            percent = (i+1)*100/(n_timesteps+1);
            printf("\33[2K\r%i %% ; timestep = %i", percent, i);
            fflush(stdout);
        }

        /* v(t+dt/2) */
        for (int j = 0; j < n_atoms; j++) {
            vel[j][0] += dt * 0.5 * F[j][0]/m;
            vel[j][1] += dt * 0.5 * F[j][1]/m;
            vel[j][2] += dt * 0.5 * F[j][2]/m;
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_atoms; j++) {
            pos[j][0] += dt * vel[j][0];
            pos[j][1] += dt * vel[j][1];
            pos[j][2] += dt * vel[j][2];
        }
        
        /* F(t+dt) */
        get_forces_AL(F,pos,L_box,n_atoms);
        
        /* v(t+dt) */
        for (int j = 0; j < n_atoms; j++) {
            vel[j][0] += dt * 0.5 * F[j][0]/m;
            vel[j][1] += dt * 0.5 * F[j][1]/m;
            vel[j][2] += dt * 0.5 * F[j][2]/m;
        }

        //increase the time
        t_ += dt;

        // calculate the energies
        E_pot = get_energy_AL(pos, L_box, n_atoms);
        E_kin = 0;
        for(int j=0;j<n_atoms;j++){
            E_kin += 0.5*m*vel[j][0]*vel[j][0];
            E_kin += 0.5*m*vel[j][1]*vel[j][1];
            E_kin += 0.5*m*vel[j][2]*vel[j][2];
        }

        // calculate the instantaneaous properties (temperature, pressure)
        T = (1/(3*n_atoms*k_B)) * E_kin;
        P = E_kin;
        for(int j=0;j<n_atoms;j++){
            P += pos[j][0]*F[j][0]/2;
            P += pos[j][1]*F[j][1]/2;
            P += pos[j][2]*F[j][2]/2;
        }
        P *= 2/(3*L_box*L_box*L_box);

        if(temp_scaling_time > 0){
            // equilibration scaling(temperature)
            double alpha_T_sqrt = sqrt(1 + 2*dt/temp_scaling_time * (T_desired - T)/T);
            for(int j=0;j<n_atoms;j++){
                vel[j][0] *= alpha_T_sqrt;
                vel[j][1] *= alpha_T_sqrt;
                vel[j][2] *= alpha_T_sqrt;
            }
        }

        if(pressure_scaling_time > 0){
            // equilibration scaling(pressure)
            double alpha_P_cbrt = cbrt(1 - kappa_T*dt/pressure_scaling_time * (P_desired - P));
            for(int j=0;j<n_atoms;j++){
                pos[j][0] *= alpha_P_cbrt;
                pos[j][1] *= alpha_P_cbrt;
                pos[j][2] *= alpha_P_cbrt;
            }
            L_box *= alpha_P_cbrt;
            a0_ *= alpha_P_cbrt;
        }

        // save the data
        fprintf(file, "%.10f, %.10e, %.10e, %.10e, %.10e, %.10e\n", t_, E_pot, E_kin, T, P, a0_);
    }
    printf("\n");
    // update the variables outside this function
    *a0 = a0_;
    *t = t_;

    fclose(file);
}

void task2(double dt, char* filename){
    double pos[n_atoms][3]; //positions
    double vel[n_atoms][3]; //velocities

    double a0 = 4.03075; // Å

    int n_timesteps = 1000;
    double t = 0;

    initialize_lattice(pos, vel, a0);
    
    simulate_and_save_data(pos, vel, &t, filename, n_timesteps, dt, &a0, 0, 0, 0, 0);
}

void task3(){
    double pos[n_atoms][3]; //positions
    double vel[n_atoms][3]; //velocities

    double a0 = 4.03075; // Å
    double dt = 10e-3; //ps
    double t = 0;
    
    double tau_T = 100*dt;
    double T_desired = 500 + celsius; // K
    int n_timesteps_temp_scaling = (5*tau_T)/dt;

    double tau_P = 10*tau_T;
    double P_desired = 1 / bar; // eV/Å^3
    int n_timesteps_pressure_scaling = (5*tau_P)/dt;

    int n_timesteps_simulation = 10000;


    initialize_lattice(pos, vel, a0);
    
    // temperature equilibriation
    simulate_and_save_data(pos, vel, &t, 
                        "data/H1_3_temp_scaling.csv", 
                        n_timesteps_temp_scaling, dt, &a0, 
                        tau_T, T_desired, 0, 0);
    // pressure equilibriation
    simulate_and_save_data(pos, vel, &t, 
                        "data/H1_3_pressure_scaling.csv", 
                        n_timesteps_pressure_scaling, dt, &a0, 
                        tau_T, T_desired, 
                        tau_P, P_desired);
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

    task1();
    task2(10e-3, "data/H1_2_small_enough.csv");
    task2(50e-3, "data/H1_2_far_too_large.csv");
    task2(17.75e-3, "data/H1_2_little_too_large.csv");
    task3();

    return 0;
}
