CC = gcc

LIBS = \
       -lm

CFLAGS_ADDRESS = \
	 	 -Werror \
		 -Wall \
	 	 -pedantic \
	 	 -g \
	 	 -fsanitize=address \
		 -fno-omit-frame-pointer

# -g was missing
CFLAGS = \
	-Werror \
	-pedantic \
	-g \
	-O0


address: clean
	$(CC) $(CFLAGS_ADDRESS) fix_me.c -o $@ $(LIBS)

gdb: clean
	$(CC) $(CFLAGS) fix_me.c -o $@ $(LIBS)

clean:
	find . -name '*.o' -exec rm -rf {} \;
	rm -rf gdb address address_example

# Added:
address_example: address_example.c
	gcc -g -fsanitize=address -fno-omit-frame-pointer address_example.c -o address_example
	./address_example
	rm address_example

