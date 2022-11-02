#TARG = common
CODE = utilities_$(TARG).a

# Set-up information ==========================================================

HOME = ..
include $(HOME)/Config/Make.cf.$(TARG)

# Files =======================================================================

SRC = monitor.cpp
TPL = 
OBJ = monitor.o
INC = debug.h monitor.h version.h
C++C = g++

# Main Rules ==================================================================

.SUFFIXES:

.SUFFIXES: .o .cpp

.cpp.o: 
	$(C++C) $(C++FLAGS) -c $<

all: $(CODE)

$(CODE): $(INC) $(SRC) $(TPL)
	@echo 'LIBRARY: Utilities'
	@make $(OBJ)
	$(AR) ruv $(CODE) $(OBJ)
	@rm -f $(OBJ)

# Miscellaneous Rules =========================================================

clean: 
	rm -f *.o *~ $(CODE)

allclean: clean
	rm -f *.a

