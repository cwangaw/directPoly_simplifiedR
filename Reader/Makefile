#TARG = common
CODE = reader_$(TARG).a

# Set-up information ==========================================================

HOME = ..
include $(HOME)/Config/Make.cf.$(TARG)

# Files =======================================================================

SRC = expr.cpp reader.cpp
TPL = 
OBJ = expr.o reader.o
INC = dimen.h expr.h reader.h units.h units.inc
C++C = g++

# Main Rules ==================================================================

.SUFFIXES:

.SUFFIXES: .o .cpp

.cpp.o: 
	$(C++C) $(C++FLAGS) -c $<

all: $(CODE)

$(CODE): $(INC) $(SRC) $(TPL)
	@echo 'LIBRARY: Reader'
	@make $(OBJ)
	$(AR) ruv $(CODE) $(OBJ)
	@rm -f $(OBJ)

# Miscellaneous Rules =========================================================

clean: 
	rm -f *.o *~ $(CODE)

allclean: clean
	rm -f *.a
# DO NOT DELETE
