LDFLAGS	= -L/home/yechangqing/.local/gsl/lib -lm -lgsl -lgslcblas
CFLAGS = -Wall -I/home/yechangqing/.local/gsl/include

test : test.o numint.o sfem.o

clean :
	rm -f *.o
