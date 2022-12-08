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

plotting: plots/H1_1.pdf plots/H1_2_small_enough.pdf plots/H1_3_equilibration.pdf plots/H1_4_equilibration.pdf task5 plots/H1_6.pdf

plots/H1_1.pdf: data/H1_1.csv python/task1.py
	python python/task1.py

plots/H1_2_small_enough.pdf: data/H1_2_small_enough.csv python/task2.py
	python python/task2.py

plots/H1_3_equilibration.pdf: data/H1_3_temp_scaling.csv python/task3.py
	python python/task3.py

plots/H1_4_equilibration.pdf: data/H1_4_temp_scaling.csv python/task4.py
	python python/task4.py

task5: python/task5.py
	python python/task5.py

plots/H1_6.pdf: data/H1_6.csv python/task6.py
	python python/task6.py


H1: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o H1 $(LIB)
	./H1

obj/%.o: src/%.c | obj
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
