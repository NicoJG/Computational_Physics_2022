CFLAGS = \
	-Werror \
	-pedantic \
	-g \
	-O0 \
	-fsanitize=address \
	-fno-omit-frame-pointer \
	-lm # -lgsl -lgslcblas

CFLAGS_FAST= \
	-Werror \
	-pedantic \
	-O3 \
	-lm # -lgsl -lgslcblas

all: simulate plot

simulate: E2.c 
	mkdir -p data
	gcc E2.c tools_.c -o E2 $(CFLAGS_FAST)
	./E2
	rm E2

plot: E2.py 
	python E2.py