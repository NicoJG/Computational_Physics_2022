#pragma once

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
                                int** nn_idxs, int idx_by_pos[n_cells][n_cells][n_cells][2]);