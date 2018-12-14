# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= mpic++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG= var.exe
${PROG} :	   	main.o wavefunc.o vectormatrixclass.o
			${CPPflags} main.o wavefunc.o vectormatrixclass.o ${LIB} -o ${PROG}

	
main.o :			main.cpp 
		        	${CPPflags} -c main.cpp


	
wavefunc.o :		wavefunc.cpp 
		        	${CPPflags} -c wavefunc.cpp	

	
vectormatrixclass.o :		vectormatrixclass.cpp 
		        	${CPPflags} -c vectormatrixclass.cpp	