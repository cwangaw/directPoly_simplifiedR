#TARG = common
CODE = mesh_$(TARG).a

# Set-up information ==========================================================

HOME = ..
include $(HOME)/Config/Make.cf.$(TARG)

# Files =======================================================================

SRC = baseObjects.cpp polyMesh.cpp
TPL = 
OBJ = baseObjects.o polyMesh.o
INC = baseObjects.h polyMesh.h
C++C = g++

# Main Rules ==================================================================

.SUFFIXES:

.SUFFIXES: .o .cpp

.cpp.o: 
	$(C++C) $(C++FLAGS) -c $<

all: $(CODE)

$(CODE): $(INC) $(SRC) $(TPL)
	@echo 'LIBRARY: Mesh'
	@make $(OBJ)
	$(AR) ruv $(CODE) $(OBJ)
	@rm -f $(OBJ)

testBase: $(CODE) testBase.cpp
	$(C++C) $(C++FLAGS) -o testBase testBase.cpp $(CODE)

testMesh: $(CODE) testMesh.cpp
	$(C++C) $(C++FLAGS) -o testMesh testMesh.cpp $(CODE)

# Miscellaneous Rules =========================================================

clean: 
	rm -f *.o *~ $(CODE) testMesh testMesh.m testMesh.txt

allclean: clean
	rm -f *.a

