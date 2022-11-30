#include <stdio.h>
#include <stdlib.h>

/*
 * Constructs the arrays that describe a binary alloy with perfect ordering
 * @n_cells - number of unit cells per direction -> N_atoms = 2*n_cells^3
 * @atype - array of shape (N_atoms,) to be filled with the atom type
 * @pos - array of shape (N_atoms,4) to be filled with the atom coordinates (x,y,z,w) 
 * x,y,z determine the unit cell and w=0,1 determines if its at the corner or the center
 * @nn_idxs - array of shape (N_atoms,8) to be filled with the indices of the nearest neighbors
 * @idx_by_pos - array of shape (n_cells, n_cells, n_cells, 2) to be filled 
 * with the index of a specific position (x,y,z,w)
 */
void construct_bcc_binary_alloy(int n_cells, int* atype, int** pos, 
                                int** nn_idxs, int idx_by_pos[n_cells][n_cells][n_cells][2]) {
    int N_atoms = 2*n_cells*n_cells*n_cells;

    int i = 0;

    // construct the lattice
    for (int x=0; x<n_cells; x++) {
    for (int y=0; y<n_cells; y++) {
    for (int z=0; z<n_cells; z++) {
    for (int w=0; w<2; w++) {
        // w==0 -> corner atom
        // w==1 -> center atom
        atype[i] = w;
        pos[i][0] = x;
        pos[i][1] = y;
        pos[i][2] = z;
        pos[i][3] = w;
        idx_by_pos[x][y][z][w] = i;
        i++;
    }
    }
    }
    }

    // assign the nearest neighbors
    // loop over all lattice positions
    for (int x=0; x<n_cells; x++) {
    for (int y=0; y<n_cells; y++) {
    for (int z=0; z<n_cells; z++) {
    for (int w=0; w<2; w++) {
        i = idx_by_pos[x][y][z][w];
        int j = 0;

        // the nearest neighbors have oposite w
        // the positions differ by (-1 or 0) if w==0 or by (0 or 1) if w==1
        // loop over all 8 nearest neighbors
        for (int dx=w-1; dx<=w; dx++) {
        for (int dy=w-1; dy<=w; dy++) {
        for (int dz=w-1; dz<=w; dz++) {
            // position of this particular nearest neighbor
            int x_nn = x+dx;
            int y_nn = y+dy;
            int z_nn = z+dz;
            int w_nn;
            if (w==0) w_nn = 1;
            if (w==1) w_nn = 0;
            // periodic boundary conditions:
            if (x_nn==-1) x_nn = n_cells-1;
            if (y_nn==-1) y_nn = n_cells-1;
            if (z_nn==-1) z_nn = n_cells-1;
            if (x_nn==n_cells) x_nn = 0;
            if (y_nn==n_cells) y_nn = 0;
            if (z_nn==n_cells) z_nn = 0;
            nn_idxs[i][j] = idx_by_pos[x_nn][y_nn][z_nn][w_nn];
            // sanity check:
            if (nn_idxs[i][j] >= N_atoms || nn_idxs[i][j]<0) {
                printf("ERROR: nn_idxs[%i][%i] = %i even though N_atoms = %i\n",i,j,nn_idxs[i][j],N_atoms);
                exit(1);
            }
            j++;
        }
        }
        }
    }
    }
    }
    }
}