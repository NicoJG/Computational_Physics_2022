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
      obj/run.o
MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_debug) 
endif

all: data/task1.csv plotting

plotting: plots/task1.pdf plots/task2.pdf plots/task3.pdf

plots/task1.pdf: python/task1.py data/task1.csv | plots
	python python/task1.py

plots/task2.pdf: python/task2.py data/task2.csv | plots
	python python/task2.py

plots/task3.pdf: python/task3.py data/task3.csv | plots
	python python/task3.py

data/task1.csv: obj H3 | data
	./H3

H3: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

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
	rm -rf data
	rm -rf plots
	rm -f ./H2

.PHONY: clean