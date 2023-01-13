#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

void
create_2D_array(
		double ***array, unsigned int row_size,
		unsigned int column_size
	       )
{

	double *asentries = (double*) malloc(row_size*column_size*sizeof(double));

	double **as=(double**) malloc(row_size*sizeof(double*));

	for ( int ix =0; ix<row_size; ++ix){
		as[ix] = asentries + column_size*ix;
	}


	for(int i=0;i<row_size;++i){
		for(int k=0;k<column_size;++k){
			as[i][k]=0;
		}
	}
	*array=as;

}
void
destroy_2D_array(
		 double **array
		)		
{
	free(array[0]);
	free(array);
	array = NULL;
}

void write_to_file(char *fname, double *numberwalkers, double *energy, int n_timesteps)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "walkers, Energy\n");
    for(int i = 0; i < n_timesteps; ++i){
	    fprintf(fp, "%f, %f\n", numberwalkers[i], energy[i]);
    }
    fclose(fp);
}

double
potential_v(double position){
    double pot = 0;
    pot = 0.5 * (1 - exp(-position))*(1 - exp(-position)); 

    return pot;
}

int
run(
    int argc,
    char *argv[]
   )
{
    // Define variables
    int N_0 = 200;
    double delta_t = 0.02; 
    double gaussian_vec[30000] ={0}; // Gaussian vector for random displacements
    double E_t[30000] ={0};
    double weight[30000] ={0};
    double **walker_matrix;
    double potential=0;
    double walker_counter[30000]={0};
    walker_counter[0]=200;
    int m=0;
    create_2D_array(&walker_matrix, 30000, 10000); 

    
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    //Initialize the GSL generator with time:
   gsl_rng_set(rng, time(NULL));


    // Equidistant placing -5 to 5.
    double equidist = -5;
    for(int i = 0; i<200; i++){
        walker_matrix[i][0] = equidist; 
        equidist += 10./199.;
    }

    // Starting energy E_t
    E_t[0]=0.5;
    int lag_counter = 0;
    double energy_sum = E_t[0];
// Create displacments 
    for( int j = 0; j<25000; j++){
    // For each step 
    lag_counter = 0;
        for(int i = 0; i <walker_counter[j]; i++){
            walker_matrix[i][j] += sqrt(delta_t) * gsl_ran_gaussian(rng, 1);
            potential=potential_v(walker_matrix[i][j]);
            weight[i]=(double)exp(-(potential-E_t[j])*delta_t);
            if(weight[i]>0){
                m=(int)(weight[i]+gsl_rng_uniform(rng));
                
                // Add new branches
                for(int n = 0; n < m; n++){
                    walker_matrix[lag_counter][j+1] = walker_matrix[i][j];
                    lag_counter += 1;
                } 

                
                //printf("%d\n", lag_counter);
                walker_counter[j+1]=lag_counter;
                m=0;  
            }
            
        }

        // Equilibration
        printf("%.5f\n",energy_sum/(j+1));
        E_t[j+1]=energy_sum/(j+1) - 0.5*log(((double)(walker_counter[j+1]))/walker_counter[0]);
        energy_sum += E_t[j+1];
        

        
       
    }

    write_to_file("walkers.csv", walker_counter, E_t, 25000 );
    gsl_rng_free(rng);
    destroy_2D_array(walker_matrix);
    
    return 0;

    
    
}
