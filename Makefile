gsl_include = /home/yechangqing/gsl/include
gsl_lib = /home/yechangqing/gsl/lib
objects = numerical_integration_on_rectangle.o

${objects}: %.o: %.c
	gcc -Wall -I${gsl_include} -c $<

numerical_integration_on_rectangle: ${objects}
	gcc -L${gsl_lib} $< -lgsl -lgslcblas -lm -o $@

clean :
	rm -f *.o
