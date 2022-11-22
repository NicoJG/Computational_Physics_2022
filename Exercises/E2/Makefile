CFLAGS = \
	-Werror \
	-pedantic \
	-g \
	-O0 \
	-fsanitize=address \
	-fno-omit-frame-pointer \
	-lm # -lgsl -lgslcblas

E2_2: E2_2.c
	gcc E2_2.c tools_.c -o E2_2 $(CFLAGS)
	./E2_2
	rm E2_2