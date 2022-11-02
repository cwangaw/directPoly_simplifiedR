HOME = .
#TARG = common
include $(HOME)/Config/Make.cf.$(TARG)

###############################################################################

PROG = directpoly

SRC = main.cpp parameterData.cpp fcns.cpp directSerendipity.cpp \
      directSerendipityFE.cpp ellipticPDE.cpp polyQuadrature.cpp \
			directMixedFE.cpp directMixed.cpp mixedPDEHybrid.cpp mixedPDEConf.cpp
INC = parameterData.h fcns.h directSerendipity.h directMixed.h ellipticPDE.h mixedPDE.h polyQuadrature.h

MAIN_MODULE = $(PROG)_$(TARG).a
OBJ = $(patsubst %.cpp, %.o, $(SRC))
MAIN_OBJECTS = $(foreach src_obj, $(OBJ), $(MAIN_MODULE)($(src_obj)))

MODULES = $(MAIN_MODULE) $(HOME)/Utilities/utilities_$(TARG).a \
              $(HOME)/Mesh/mesh_$(TARG).a $(HOME)/Reader/reader_$(TARG).a
C++C = g++

###############################################################################

.SUFFIXES:
.SUFFIXES: .o .cpp

.cpp.o: 
	$(C++C) $(C++FLAGS) -c $<

$(PROG).$(TARG): $(MAIN_OBJECTS) $(MODULES)
	$(C++C) -o $@ $(LDFLAGS) $(MODULES) $(LIBS)

.force:

all init:
	make $(MAIN_MODULE)
	cd $(HOME)/Utilities; make
	cd $(HOME)/Mesh; make
	cd $(HOME)/Reader; make
	make

$(MAIN_MODULE): $(MAIN_OBJECTS)

$(MAIN_MODULE)(%.o) : %.cpp $(INC)
	$(C++C) -c $(C++FLAGS) $<
	$(AR) $(ARFLAGS) $@  $%
	$(RM) $*.o

$(HOME)/Utilities/utilities_$(TARG).a:
	cd $(HOME)/Utilities; make

$(HOME)/Mesh/mesh_$(TARG).a:
	cd $(HOME)/Mesh; make

$(HOME)/Reader/reader_$(TARG).a:
	cd $(HOME)/Reader; make

###############################################################################

tar: $(SRC) $(INC)
	tar cvf $(PROG).tar Makefile infile_template \
          $(SRC) $(INC) \
          $(HOME)/Utilities $(HOME)/Mesh $(HOME)/Reader $(HOME)/Config \
          infile $(HOME)/Polymesher $(HOME)/test

allclean: clean
	cd $(HOME)/Utilities; make allclean
	cd $(HOME)/Mesh; make allclean
	cd $(HOME)/Reader; make allclean
	cd $(HOME)/Config; rm -f *~

clean:
	rm -f *~ echo logfile *.o *.a $(PROG).$(TARG) $(PROG).tar debugfile
