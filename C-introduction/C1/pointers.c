#include <stdio.h>

int
main()
{ 
    int m = 2;
    int *pm = &m;
    /* A pointer to a pointer */
    int **ppm = &pm;

    printf("m = %i\n", m + 1);
    printf("(*pm) + 1 = %i\n", (*pm) + 1);
    printf("pm = %p\n", pm);
    printf("*(&m) + 1 = %i\n", *(&m) + 1);
    printf("(&m) = %p\n", (&m));
    printf("ppm = %p\n", ppm);
    printf("(*ppm) = %p\n", *ppm);
    printf("&(*ppm) = %p\n", &(*ppm));
    printf("*(*ppm) + 1 = %i\n", *(*ppm) + 1);

    /* Usefulness of pointers: more control over memory management */
    return 0;
}
