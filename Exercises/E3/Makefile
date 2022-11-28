CFLAGS = \
	-Werror \
	-pedantic \
	-g \
	-O0 \
	-fsanitize=address \
	-fno-omit-frame-pointer \
	-lm -lgsl -lgslcblas

CFLAGS_FAST= \
	-Werror \
	-pedantic \
	-O3 \
	-lm -lgsl -lgslcblas

all: simulate plot

simulate: E3.c 
	mkdir -p data
	gcc E3.c -o E3 $(CFLAGS)
	./E3
	rm E3

plot: E3.py 
	mkdir -p plots
	python E3.py