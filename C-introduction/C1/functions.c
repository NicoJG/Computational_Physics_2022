#include <stdio.h>


// Function prototypes
void increment_wrong(int);
void increment_correct(int *);


int
main()
{
    int a = 0;
    increment_wrong(a);
    printf("In main after \"increment wrong\" a = %i\n", a);
    printf("--Address in main after \"increment wrong\" &a = %p\n", &a);
    increment_correct(&a);
    printf("In main after \"increment correct\" a = %i\n", a);
    printf("--Address in main after \"increment correct\" &a = %p\n", &a);
    return 0;
}

// Function definitions
void
increment_wrong(
		int a
	       )
{
    a += 1;
    printf("In \"increment wrong\" a = %i\n", a);
    printf("--Address in \"increment wrong\" &a = %p\n", &a);
}

void
increment_correct(
		  int *a
		 )
{
    (*a) += 1;
    printf("In \"increment correct\" a = %i\n", *a);
    printf("--Address in \"increment correct\" a = %p\n", a);
}
