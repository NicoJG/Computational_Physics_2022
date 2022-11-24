CC = gcc

CFLAGS = \
	 -Iinclude

CFLAGS_debug = \
	-g \
	-fsanitize=address \
	-fno-omit-frame-pointer \
	-O0

CFLAGS_OPT = \
	     -O2 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas

OBJ = \
      obj/tools.o \
      obj/lattice.o \
      obj/potential.o \
      obj/run.o
MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT) 
endif

all: H1 plotting | data plots

plotting: plots/H1_1.pdf plots/H1_2_small_enough.pdf plots/H1_3_temp_scaling.pdf

plots/H1_1.pdf: data/H1_1.csv python/task1.py
	python python/task1.py

plots/H1_2_small_enough.pdf: data/H1_2_small_enough.csv python/task2.py
	python python/task2.py

plots/H1_3_temp_scaling.pdf: data/H1_3_temp_scaling.csv python/task3.py
	python python/task3.py

H1: obj _H1
	./H1

_H1: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o H1 $(LIB)

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
