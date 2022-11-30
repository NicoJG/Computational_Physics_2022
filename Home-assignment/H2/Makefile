CC = gcc

CFLAGS = \
	 -Iinclude \
	 -Werror

CFLAGS_debug = \
	-g \
	-fsanitize=address \
	-fno-omit-frame-pointer \
	-O0

CFLAGS_OPT = \
	     -O3 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas

OBJ = \
      obj/tools.o \
      obj/lattice_tools.o \
      obj/run.o
MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT) 
endif

all: H2 plotting | data plots

plotting: plots/task1_P.pdf

plots/task1_P.pdf: python/task1.py
	python python/task1.py

H2: obj _H2
	./H2

_H2: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o H2 $(LIB)

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

data:
	mkdir -p data

plots:
	mkdir -p plots

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f program run-test
	rm -rf obj obj_test

.PHONY: clean
