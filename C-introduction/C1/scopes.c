#include <stdio.h>
    
int a = 2; // global

int
main()
{
    int b; // local
    int c = 3;
    {
    int c = 100; // declares new variable (not accessable outside of this scope)
	b = 1; // overwrites the variable declared before
    }
    printf("a: %i\n", a);
    printf("b: %i\n", b);
    printf("c: %i\n", c);
    return 0;
}
