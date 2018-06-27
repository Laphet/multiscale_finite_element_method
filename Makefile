LDFLAGS	= -L/home/yechangqing/gsl/lib -lm -lgsl -lgslcblas
CFLAGS = -Wall -I/home/yechangqing/gsl/include

test : test.o numint.o sfem.o

clean :
	rm -f *.o
