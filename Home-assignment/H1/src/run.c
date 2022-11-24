#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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
    fprintf(file, "# n_cells = %i\n",n_cells);
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

int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code

    // Task 1:
    task1();

    return 0;
}
