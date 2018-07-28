CFLAGS	         = -I/home/yechangqing/.local/gsl/include/
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
NP               = 1
USER_LIBS= 	-L/home/yechangqing/.local/gsl/lib/ -Wl,-rpath,/home/yechangqing/.local/gsl/lib/ -lgsl -lgslcblas -lm

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

test : test.o numint.o sfem.o
	-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} ${USER_LIBS}
