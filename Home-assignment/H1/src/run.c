#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"

void task1(){
    // constants
    int n_cells = 4; // number of unit cells
    int n_atoms = 4*n_cells*n_cells*n_cells;
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

void task2(double dt, char* filename){
    int n_cells = 4; // number of unit cells
    int n_atoms = 4*n_cells*n_cells*n_cells;
    double pos[n_atoms][3]; //positions
    double vel[n_atoms][3]; //velocities
    double F[n_atoms][3]; //forces

    init_matrix_stack(n_atoms,3,pos,0);
    init_matrix_stack(n_atoms,3,vel,0);
    init_matrix_stack(n_atoms,3,F,0);

    double m_asu = 9649; // u
    double m = 26.98153853 / m_asu; // m_asu (https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Al)

    double a0 = 4.03075; // Å
    double L_box = n_cells*a0;

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
    // prepare output file
    FILE* file = fopen(filename, "w");
    fprintf(file, "# {\"dt\": %.5e}\n",dt);
    fprintf(file, "# t[ps], E_pot[eV], E_kin[eV]\n");

    // velocity verlet:
    //double dt = 5e-3 // ps
    int n_timesteps = 10000;
    double t_max = (n_timesteps+1)*dt;
    printf("dt = %.5e\n", dt);
    printf("n_timesteps = %i\n", n_timesteps);
    printf("t_max = %.5f\n", t_max);

    int percent = -1; // for displaying progress

    double E_pot = 0; // eV
    double E_kin = 0; // eV

    double t = 0; // ps
    get_forces_AL(F,pos,L_box,n_atoms);

    for(int i=1;i<n_timesteps+1;i++){

        // show progress
        if(i*100/(n_timesteps+1) != percent) {
            percent = i*100/(n_timesteps+1);
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
        t += dt;

        // calculate the energies
        E_pot = get_energy_AL(pos, L_box, n_atoms);
        E_kin = 0;
        for(int j=0;j<n_atoms;j++){
            E_kin += 0.5*m*vel[j][0]*vel[j][0];
            E_kin += 0.5*m*vel[j][1]*vel[j][1];
            E_kin += 0.5*m*vel[j][2]*vel[j][2];
        }

        // save the data
        fprintf(file, "%.10f, %.10e, %.10e\n", t, E_pot, E_kin);
    }
    printf("\n");

    fclose(file);
    gsl_rng_free(rng);
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
    task2(15e-3, "data/H1_2_little_too_large.csv");

    return 0;
}
