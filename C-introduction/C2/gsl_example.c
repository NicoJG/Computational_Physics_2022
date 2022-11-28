#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

int
main()
{
    // set up the random number generator
    int seed = 42;
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, seed); 

    // Generate some random numbers
    for (int i = 0; i < 10; i++){
	double u = gsl_rng_uniform(rng);
	double g = gsl_ran_gaussian(rng, 1);
	printf ("%.5f %.5f\n", u, g);
    }

    // free the allocated random number generator
    gsl_rng_free(rng);

    return 0;
}
