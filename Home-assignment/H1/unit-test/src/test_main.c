#include <stdlib.h>
#include <string.h>

#include "test_main.h"

/* ************************************
 * Here follows the test we want
 * to run
 * ***********************************/
START_TEST(test_what_you_wanna_test)
{
    // empty
}


int
main()
{
    test_setup("testing", "core");
    
    // Tests
    add_test(test_what_you_wanna_test);
    
    test_teardown();
    return 0;
}
