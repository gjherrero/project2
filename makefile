clean:
	rm *.o 
	rm *.x
mainCode: compile execute
compile: functionsTest.o mainTest.o
	c++ -O3 -o test.x mainTest.o functionsTest.o -larmadillo -llapack -lblas 
execute:
	./test.x


PROG= testcode.x
tests: ${PROG} executeTests
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
${PROG} :		unitTestsMain.o  unitTests.o functionsTest.o
			${CPPflags} -o ${PROG} unitTestsMain.o unitTests.o functionsTest.o ${LIB} -o ${PROG}

#This creates the object file for the associated .cpp file
#avoid different intstancies 
%.o:			%.cpp 
					${CPPflags} -c $< -o $@
executeTests: 
	./testcode.x

#unitTestsMain.o:			unitTestsMain.cpp 
#					${CPPflags} -c unitTestsMain.cpp

#unitTests.o:			unitTests.cpp
#					${CPPflags} -c unitTests.cpp

#functionsTest.o:		functionsTest.cpp 
#					${CPPflags} -c functionsTest.cpp
