#include <stdlib.h>
#include <stdio.h>

int
main()
{
    //float *heap_array = (float *)malloc(sizeof(float) * 2097152 * 1000); // works
    float *heap_array = (float *)malloc(sizeof(float) * 2097152 * 2000); // fails
    if(heap_array == NULL){
	perror("malloc failed");
    printf("Program continues.\n");
	exit(1);
    }
    heap_array[99] = 99;
    printf("%f\n", heap_array[99]);
    free(heap_array);
    return 0;
}
