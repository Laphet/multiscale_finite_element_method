CPPFLAGS = -Wall
CFLAGS = -I/home/yechangqing/gsl/include
LDFLAGS	= -L/home/yechangqing/gsl/lib -lm -lgsl -lgslcblas

test : test.o numint.o

clean :
	rm -f *.o
