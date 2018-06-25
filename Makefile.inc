CPP		= /share/soft/mvapich2-2.3b_gcc/bin/mpicc -E
CPPFLAGS	= -I/home/yechangqing/phg-0.9.3/include  -D__PHG__
BUILD_CPPFLAGS	= -DMPI_NO_CPPBIND   -DVTK_DIR="\"\""  -DGZIP_PROG="\"/usr/bin/gzip\"" -DBZIP2_PROG="\"/usr/bin/bzip2\""

CC		= /share/soft/mvapich2-2.3b_gcc/bin/mpicc 
CFLAGS		= -g -O2 -Wall -Wmissing-field-initializers -fopenmp
BUILD_CFLAGS	= -fPIC
USER_CFLAGS	= -I/home/yechangqing/gsl/include

CXX		= /share/soft/mvapich2-2.3b_gcc/bin/mpicxx 
CXXFLAGS	= -g -O2 -Wno-deprecated -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -fopenmp
BUILD_CXXFLAGS	= -fPIC
USER_CXXFLAGS	=

FC		= /share/soft/mvapich2-2.3b_gcc/bin/mpif90 
FCFLAGS		= -g -O2
USER_FCFLAGS	=

F77		= /share/soft/mvapich2-2.3b_gcc/bin/mpif77 
FFLAGS		= -g -O2
USER_FFLAGS	=

LINKER		= /share/soft/mvapich2-2.3b_gcc/bin/mpicxx 
LDFLAGS		= -L/home/yechangqing/phg-0.9.3/src  -fopenmp
USER_LDFLAGS	=

LIBS		= -Wl,-rpath,/home/yechangqing/phg-0.9.3/src -lphg  -lmpfr  -lgmp  -lm  -L/soft/mvapich2-2.3b-gcc/lib -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/share/soft/intel_2018_update1/compilers_and_libraries_2018.1.163/linux/ipp/lib/intel64 -L/share/soft/intel_2018_update1/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin -L/share/soft/intel_2018_update1/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin -L/share/soft/intel_2018_update1/compilers_and_libraries_2018.1.163/linux/tbb/lib/intel64/gcc4.7 -L/share/soft/intel_2018_update1/compilers_and_libraries_2018.1.163/linux/daal/lib/intel64_lin -L/share/soft/intel_2018_update1/compilers_and_libraries_2018.1.163/linux/daal/../tbb/lib/intel64_lin/gcc4.4 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../.. -lgfortran -lm -lmpifort -lmpi -lquadmath  -lgfortran -lm -lmpifort -lmpi -lquadmath -Wl,-rpath,/usr/local/lib -Wl,-rpath,/usr/local/lib64
USER_LIBS	= -L/home/yechangqing/gsl/lib -lgsl -lgslcblas

AR		= /usr/bin/ar 
RANLIB		= /usr/bin/ranlib 
BUILD_SHARED	= gcc -shared
BUILD_SHARED_LIBS = 
LIB_SUFFIX	= .so

.o.c:
	${CC} ${CFLAGS} ${CPPFLAGS} ${USER_CFLAGS} -c $*.c

.cxx.c:
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${USER_CXXFLAGS} -c $*.cxx

.f.o:
	${FC} ${FCFLAGS} ${USER_FCFLAGS} -c $*.f

% : %.o
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

% : %.c
	${CC} ${CFLAGS} ${CPPFLAGS} ${USER_CFLAGS} -c $*.c
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

% : %.cxx
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${USER_CXXFLAGS} -c $*.cxx
	${CXX} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

% : %.f
	${FC} ${FCFLAGS} ${USER_FCFLAGS} -c $*.f
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

CPPFLAGS:
	@echo $(CPPFLAGS)

CFLAGS:
	@echo $(CFLAGS)

CXXFLAGS:
	@echo $(CXXFLAGS)

LDFLAGS:
	@echo $(LDFLAGS)

LIBS:
	@echo $(LIBS)
