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
const double bar = 1e-4 * GPa; // eV/Å^3
const double celsius = 273.15; // K

// Simulation/Physics constants
const double m = 26.98153853 / m_asu; // m_asu (https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Al)
const double k_B = 8.617333262e-5; // eV/K (https://en.wikipedia.org/wiki/Boltzmann_constant#Value_in_different_units)
const double kappa_T = 0.01385 / GPa; // Å^3/eV https://www.knowledgedoor.com/2/elements_handbook/isothermal_compressibility.html
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

void velocity_verlet_timestep(double* t, double dt, double a0,
                            double pos[][3], double vel[][3], double F[][3]) {
    double L_box = n_cells*a0;
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
    (*t) += dt;
}

void simulate_and_save_data(double pos[][3], double vel[][3], double* t, 
                            char* filename, 
                            int n_timesteps, double dt, double* a0,
                            double temp_scaling_time, double T_desired, 
                            double pressure_scaling_time, double P_desired,
                            int n_save_positions){
    // if the scaling times are <= 0 then no corresponding scaling happens (no equilibration)

    double run_time = (n_timesteps+1)*dt;

    // prepare output file
    FILE* file = fopen(filename, "w");
    fprintf(file, "# {\"dt\": %.5e, \"n_timesteps\": %i, \"run_time\": %.5e, \n",dt, n_timesteps, run_time);
    fprintf(file, "# \"temp_scaling_time\": %.5e, \"T_desired\": %.5e, \n", temp_scaling_time, T_desired);
    fprintf(file, "# \"pressure_scaling_time\": %.5e, \"P_desired\": %.5e, \n", pressure_scaling_time, P_desired/bar);
    fprintf(file, "# \"n_atoms\": %i, \"n_save_positions\": %i}\n", n_atoms, n_save_positions);

    // table header
    fprintf(file, "# t[ps], E_pot[eV], E_kin[eV], T[K], P[bar], a0[Å]");
    for(int j=0;j<n_save_positions;j++){
        fprintf(file, ", x_1[Å]");
        fprintf(file, ", y_1[Å]");
        fprintf(file, ", z_1[Å]");
    }
    fprintf(file, "\n");

    // velocity verlet:
    printf("Simulate \"%s\"...\n", filename);
    printf("dt = %.5e ; n_timesteps = %i ; run_time = %.5f\n", dt, n_timesteps, run_time);

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

        velocity_verlet_timestep(&t_, dt, a0_, pos, vel, F);

        // calculate the energies
        E_pot = get_energy_AL(pos, L_box, n_atoms);
        E_kin = 0;
        for(int j=0;j<n_atoms;j++){
            E_kin += 0.5*m*vel[j][0]*vel[j][0];
            E_kin += 0.5*m*vel[j][1]*vel[j][1];
            E_kin += 0.5*m*vel[j][2]*vel[j][2];
        }

        // calculate the instantaneaous properties (temperature, pressure)
        T = (2/(3*n_atoms*k_B)) * E_kin;
        P = (n_atoms * k_B * T + get_virial_AL(pos, L_box, n_atoms)) / (L_box*L_box*L_box);

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

        // Task 6: radial distribution function
        double

        P = P / bar; // GPa 

        // save the data
        fprintf(file, "%.10f, %.10e, %.10e, %.10e, %.10e, %.10e", t_, E_pot, E_kin, T, P, a0_);
        for(int j=0;j<n_save_positions;j++){
            fprintf(file, ", %.5e", pos[j][0]);
            fprintf(file, ", %.5e", pos[j][1]);
            fprintf(file, ", %.5e", pos[j][2]);
        }
        fprintf(file, "\n");
    }
    printf("\n");
    // update the variables outside this function
    *a0 = a0_;
    *t = t_;

    fclose(file);
}

void save_system_state(char* file_path, int n_atoms, double pos[n_atoms][3], double vel[n_atoms][3], double a0, double t, double T, double P) {
    FILE* file = fopen(file_path, "w");
    // Header
    fprintf(file, "# t[ps] = %.10f\n", t);
    fprintf(file, "# a0[Å] = %.10f\n", a0);
    fprintf(file, "# T[K] = %.10f\n", T);
    fprintf(file, "# P[bar] = %.10f\n", P/bar);
    fprintf(file, "# n_atoms = %i\n", n_atoms);
    fprintf(file, "# pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]\n");
    for (int i=0; i<n_atoms; i++) {
        fprintf(file, "%.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n", 
            pos[i][0], pos[i][1], pos[i][2], vel[i][0], vel[i][1], vel[i][2]);
    }
    fclose(file);
}

void load_system_state(char* file_path, int n_atoms, double pos[n_atoms][3], double vel[n_atoms][3], double* a0, double* t, double* T, double* P) {
    FILE* file = fopen(file_path, "r");
    char buf[255];
    int n_atoms_test = 0;
    int res;

    res = fscanf(file, "# t[ps] = %lf\n", t);
    res = fscanf(file, "# a0[Å] = %lf\n", a0);
    res = fscanf(file, "# T[K] = %lf\n", T);
    res = fscanf(file, "# P[bar] = %lf\n", P);
    (*P) *= bar;
    res = fscanf(file, "# n_atoms = %i\n", &n_atoms_test);
    if (n_atoms_test != n_atoms) {
        perror("ERROR: n_atoms_test != n_atoms");
        exit(1);
    }
    res = fscanf(file, "# pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]\n");
    for (int i=0; i<n_atoms; i++) {
        res = fscanf(file, "%lf, %lf, %lf, %lf, %lf, %lf\n", 
            &pos[i][0], &pos[i][1], &pos[i][2], &vel[i][0], &vel[i][1], &vel[i][2]);
    }
    fclose(file);
}

void task2(double dt, char* filename){
    double pos[n_atoms][3]; //positions
    double vel[n_atoms][3]; //velocities

    double a0 = 4.03139; // Å

    double t_max = 15;
    int n_timesteps = t_max/dt;

    double t = 0;

    initialize_lattice(pos, vel, a0);
    
    simulate_and_save_data(pos, vel, &t, filename, n_timesteps, dt, &a0, 0, 0, 0, 0, 0);
}

void task3(){
    double pos[n_atoms][3]; //positions
    double vel[n_atoms][3]; //velocities

    double a0 = 4.03139; // Å
    double dt = 5e-3; //ps
    double t = 0;
    
    double tau_T = 100*dt;
    double T_desired = 500 + celsius; // K
    int n_timesteps_temp_scaling = (6*tau_T)/dt;

    double tau_P = 3*tau_T;
    double P_desired = 1 * bar; // eV/Å^3
    int n_timesteps_pressure_scaling = (16*tau_P)/dt;

    double dt_simulation = 2e-3; //ps
    int n_timesteps_simulation = 10000;

    printf("P_desired = %.5e eV/Å^3\n", P_desired);
    printf("kappa_T = %.5e eV/Å^3\n", kappa_T);

    initialize_lattice(pos, vel, a0);
    
    // temperature equilibriation
    simulate_and_save_data(pos, vel, &t, 
                        "data/H1_3_temp_scaling.csv", 
                        n_timesteps_temp_scaling, dt, &a0, 
                        tau_T, T_desired, 0, 0, 0);
    // pressure equilibriation
    simulate_and_save_data(pos, vel, &t, 
                        "data/H1_3_pressure_scaling.csv", 
                        n_timesteps_pressure_scaling, dt, &a0, 
                        tau_T, T_desired, 
                        tau_P, P_desired, 0);

    // simulate
    simulate_and_save_data(pos, vel, &t, 
                        "data/H1_3_after_scaling.csv", 
                        n_timesteps_simulation, dt_simulation, &a0, 
                        0,0,0,0,
                        5);

    save_system_state("data/H1_solid_state.csv", n_atoms, pos, vel, a0, t, T_desired, P_desired);
}

void task4(){
    double pos[n_atoms][3]; //positions
    double vel[n_atoms][3]; //velocities

    double a0 = 4.03139; // Å
    double dt = 5e-3; //ps
    double t = 0;
    
    double tau_T = 100*dt;
    double tau_P = 3*tau_T;

    double T_desired_melting = 1500 + celsius; // K
    double P_desired = 1 * bar; // eV/Å^3
    double T_desired = 700 + celsius; // K


    int n_timesteps_temp_scaling = (4*tau_T)/dt;
    int n_timesteps_pressure_scaling = (20*tau_P)/dt;
    int n_timesteps_temp_decreasing = (20*tau_P)/dt;

    double dt_simulation = 2e-3; //ps
    int n_timesteps_simulation = 10000;

    printf("P_desired = %.5e eV/Å^3\n", P_desired);
    printf("kappa_T = %.5e eV/Å^3\n", kappa_T);

    initialize_lattice(pos, vel, a0);
    
    // temperature equilibriation
    simulate_and_save_data(pos, vel, &t, 
                        "data/H1_4_temp_scaling.csv", 
                        n_timesteps_temp_scaling, dt, &a0, 
                        tau_T, T_desired_melting, 0, 0, 0);
    // pressure equilibriation while melting
    simulate_and_save_data(pos, vel, &t, 
                        "data/H1_4_pressure_scaling.csv", 
                        n_timesteps_pressure_scaling, dt, &a0, 
                        tau_T, T_desired_melting, 
                        tau_P, P_desired, 0);

    // temperature equilibriation to lower temperature while it is melted
    simulate_and_save_data(pos, vel, &t, 
                        "data/H1_4_temp_decreasing.csv", 
                        n_timesteps_temp_decreasing, dt, &a0, 
                        tau_T, T_desired, 
                        tau_P, P_desired, 0);

    // simulate
    simulate_and_save_data(pos, vel, &t, 
                        "data/H1_4_after_scaling.csv", 
                        n_timesteps_simulation, dt_simulation, &a0, 
                        0,0,0,0,
                        5);

    save_system_state("data/H1_liquid_state.csv", n_atoms, pos, vel, a0, t, T_desired, P_desired);
}

void task6() {
    double pos[n_atoms][3]; //positions
    double vel[n_atoms][3]; //velocities
    double F[n_atoms][3]; //forces
    double a0, t, T, P;

    load_system_state("data/H1_liquid_state.csv", n_atoms, pos, vel, &a0, &t, &T, &P);

    double dt = 5e-3; //ps
    int n_timesteps = 10000;
    double run_time = (n_timesteps+1)*dt;

    double L_box = n_cells*a0;
    int percent = -1; // for displaying progress

    // prepare the histogram
    int n_bins = 300;
    double r_max = L_box/2;
    double dr = r_max/n_bins;
    double N_r[n_bins];
    for (int i_bin=0; i_bin<n_bins; i_bin++) {
        N_r[i_bin] = 0;
    }

    printf("Task 6: Simulate and calculate the pair distance histogram...\n");
    printf("dt = %.5e ; n_timesteps = %i ; run_time = %.5f\n", dt, n_timesteps, run_time);
    
    get_forces_AL(F,pos,L_box,n_atoms);

    for (int i_step=0; i_step<n_timesteps; i_step++) {
        // show progress
        if((i_step+1)*100/(n_timesteps+1) != percent) {
            percent = (i_step+1)*100/(n_timesteps+1);
            printf("\33[2K\r%i %% ; timestep = %i", percent, i_step);
            fflush(stdout);
        }

        velocity_verlet_timestep(&t, dt, a0, pos, vel, F);

        // calculate the pair distances for the <N(r)> histogram
        for (int k=0; k<n_atoms; k++) {
        for (int l=k+1; l<n_atoms; l++) {
            // distance in the minimum image convention
            double dx[3];
            for (int i_x=0; i_x<3; i_x++) {
                dx[i_x] = pos[k][i_x] - pos[l][i_x];
                dx[i_x] = dx[i_x] - L_box*round(dx[i_x]/L_box);
            }
            double dist = vector_norm(dx, 3);

            // add it to the right bin (twice to accommodate for both N(r))
            if (dist < r_max) {
                int i_bin = (int) (dist/dr);
                N_r[i_bin] += 2;
            }
        }
        }
    }
    printf("\n");

    // take the average instead of the absolute count
    for (int i_bin=0; i_bin<n_bins; i_bin++) {
        N_r[i_bin] /= (double) n_atoms*n_timesteps;
    }
    
    // write to a file
    FILE* file = fopen("data/H1_6.csv", "w");
    fprintf(file, "# {\"n_atoms\": %i, \"L_box[Å]\": %.10f}\n", n_atoms, L_box);
    fprintf(file, "# r[Å], <N(r)>\n");
    for (int i_bin=0; i_bin<n_bins; i_bin++) {
        fprintf(file, "%.10f, %.10f\n", (i_bin+0.5)*dr, N_r[i_bin]);
    }
    fclose(file);
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

    //task1();
    //task2(5e-3, "data/H1_2_small_enough.csv");
    //task2(50e-3, "data/H1_2_far_too_large.csv");
    //task2(15e-3, "data/H1_2_little_too_large.csv");
    //task3();
    //task4();
    task6();

    return 0;
}
