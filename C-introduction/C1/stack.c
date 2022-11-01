int
main()
{
    // float is 4 bytes, stacksize is 8192 kbytes
    // => 8192*1024 = 8388608 bytes 
    // => 8388608/4 = 2097152 floats
    float stack_array[2097152]; //fails
    // float stack_array[2097152 - 2000]; //sometimes fails
    // float stack_array[2097152 - 3000]; //does not fail
    stack_array[99] = 0;
    return 0;
}
