CFLAGS = \
	-Werror \
	-pedantic \
	-g \
	-fsanitize=address \
	-fno-omit-frame-pointer \
	-lm -lgsl -lgslcblas \
	-O0

1: E1code1.c E1code2.py
	gcc E1code1.c -o E1code1 -lm
	./E1code1
	python E1code2.py

2: E1code3.c fft.c signal.csv
	gcc E1code3.c fft.c -o E1code3 -lm -lgsl -lgslcblas
	./E1code3
	python E1code3_plot.py
5: E1code5.c
	gcc E1code5.c fft.c $(CFLAGS) -o E1code5
	./E1code5
	python E1code5_plot.py
6: E1code6.c
	gcc E1code6.c fft.c $(CFLAGS) -o E1code6
	./E1code6
	python E1code6_plot.py