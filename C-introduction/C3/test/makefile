CC = gcc

# The -Iinclude directive tells the compiler to look for
# headers files inuder the include folder
CFLAGS = \
	 -pedantic \
	 -Wall \
	 -Werror \
	 -g \
	 -Iinclude 
OBJ = \
      obj/vector.o
MAIN = \
       obj/main.o

CFLAGS_OPT = \
	     -O2 \
	     -march=native

# for math.h
LIB += \
     -lm

ifeq ($(MAKECMDGOALS),test)
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

program: obj _program

_program: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o program $(LIB)

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f program run-test
	rm -rf obj


.PHONY: clean
